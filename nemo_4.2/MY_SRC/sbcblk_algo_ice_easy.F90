MODULE sbcblk_algo_ice_easy
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ice_easy  ***
   !!       Computes turbulent components of surface fluxes over sea-ice
   !!
   !!   What is the "EASY" algorithm ?
   !!    => Given a constant value for the neutral coefficients C_D_N, C_E_N and C_H_N over sea-ice
   !!       C_D, C_E and C_H consistent with the near-surface atmospheric stability
   !!    ==> which is already a better option than using a constant value for C_D, C_E and C_H as in
   !!        the default (and simplest) option in NEMO for instance.
   !!    ==> the only room for XXX is the pick of the stability functions: we use those of Andreas 2005 for now...
   !!
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (usually 2m) to zu (usually 10m) if needed
   !!   * the "effective" bulk wind speed at zu: Ub (including gustiness contribution in unstable conditions)
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!       Routine turb_ice_easy maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!            Author: Laurent Brodeau, Summer 2023
   !!
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE par_oce,  ONLY: jpi, jpj, Nis0, Nie0, Njs0, Nje0, nn_hls, ntsi, ntsj, ntei, ntej
   USE lib_mpp,  ONLY: ctl_stop         ! distribued memory computing library
   USE phycst          ! physical constants
   USE sbc_phy         ! Catalog of functions for physical/meteorological parameters in the marine boundary layer

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: turb_ice_easy

   INTEGER , PARAMETER ::   nb_iter = 8        ! number of itterations

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ice_easy( zt, zu, Ts_i, t_zt, qs_i, q_zt, U_zu,        &
      &                      CdN, ChN, CeN,                               &
      &                      Cd_i, Ch_i, Ce_i, t_zu_i, q_zu_i,            &
      &                      xz0, xu_star, xL, xUN10 )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ice_easy  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to:
      !!   Andreas, E.L., Jordan, R.E. & Makshtas, A.P. Parameterizing turbulent exchange over sea ice: the ice station weddell results.
      !!   Boundary-Layer Meteorology 114, 439–460 (2005). https://doi.org/10.1007/s10546-004-1414-7
      !!
      !!           If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!           Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  Ts_i  : surface temperature of sea-ice                         [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  qs_i  : saturation specific humidity at temp. Ts_i over ice    [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!    * CdN     : neutral-stability drag coefficient
      !!    * ChN     : neutral-stability sensible heat coefficient
      !!    * CeN     : neutral-stability evaporation coefficient
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd_i   : drag coefficient over sea-ice
      !!    *  Ch_i   : sensible heat coefficient over sea-ice
      !!    *  Ce_i   : sublimation coefficient over sea-ice
      !!    *  t_zu_i : pot. air temp. adjusted at zu over sea-ice             [K]
      !!    *  q_zu_i : spec. hum. of air adjusted at zu over sea-ice          [kg/kg]
      !!
      !! OPTIONAL OUTPUT:
      !! ----------------
      !!    * xz0     : return the aerodynamic roughness length (integration constant for wind stress) [m]
      !!    * xu_star : return u* the friction velocity                    [m/s]
      !!    * xL      : return the Obukhov length                          [m]
      !!    * xUN10   : neutral wind speed at 10m                          [m/s]
      !!
      !! ** Author: L. Brodeau, July 2023 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in )                     :: zt    ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in )                     :: zu    ! height for U_zu                             [m]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: Ts_i  ! ice surface temperature                [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: t_zt  ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: qs_i  ! sat. spec. hum. at ice/air interface    [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: q_zt  ! spec. air humidity at zt               [kg/kg]
      REAL(wp), INTENT(in ), DIMENSION(jpi,jpj) :: U_zu  ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(in )                 :: CdN, ChN, CeN
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Cd_i  ! drag coefficient over sea-ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ch_i  ! transfert coefficient for heat over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: Ce_i  ! transfert coefficient for sublimation over ice
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: t_zu_i ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: q_zu_i ! spec. humidity adjusted at zu           [kg/kg]
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xz0  ! Aerodynamic roughness length   [m]
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xu_star  ! u*, friction velocity
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xL  ! zeta (zu/L)
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj), OPTIONAL :: xUN10  ! Neutral wind at zu
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Ubzu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp0, ztmp1      ! temporary stuff
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dt_zu, dq_zu
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: u_star, t_star, q_star
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zeta_u, zeta_t           ! stability parameter at height zu
      !!
      REAL(wp) :: zsqrtCDN, zlog1, zlog2
      INTEGER :: jit
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !!
      LOGICAL :: lreturn_z0=.FALSE., lreturn_ustar=.FALSE., lreturn_L=.FALSE., lreturn_UN10=.FALSE.
      !!
      CHARACTER(len=40), PARAMETER :: crtnm = 'turb_ice_easy@sbcblk_algo_ice_easy.f90'
      !!----------------------------------------------------------------------------------
      ALLOCATE (  Ubzu(jpi,jpj), u_star(jpi,jpj),  t_star(jpi,jpj),  q_star(jpi,jpj),  &
         &      zeta_u(jpi,jpj),  dt_zu(jpi,jpj),   dq_zu(jpi,jpj),  &
         &      ztmp0(jpi,jpj),   ztmp1(jpi,jpj)  )

      lreturn_z0    = PRESENT(xz0)
      lreturn_ustar = PRESENT(xu_star)
      lreturn_L     = PRESENT(xL)
      lreturn_UN10  = PRESENT(xUN10)

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )
      IF( .NOT. l_zt_equal_zu )  ALLOCATE( zeta_t(jpi,jpj) )

      zsqrtCDN = SQRT(CdN)
      zlog1 = LOG(zt/zu)
      zlog2 = LOG(zu/10._wp)

      !! Scalar wind speed cannot be below 0.2 m/s
      Ubzu(:,:) = MAX( U_zu, wspd_thrshld_ice )

      !! First guess of temperature and humidity at height zu:
      t_zu_i = MAX( t_zt ,   100._wp )   ! who knows what's given on masked-continental regions...
      q_zu_i = MAX( q_zt , 0.1e-6_wp )   !               "

      !! First guess of transfer coefficients:
      Cd_i(:,:) = CdN
      Ch_i(:,:) = ChN
      Ce_i(:,:) = CeN

      !! ITERATION BLOCK
      DO jit = 1, nb_iter

         !! Air-Ice differences:
         dt_zu(:,:) = t_zu_i - Ts_i ; !  dt_zu = SIGN( MAX(ABS(dt_zu),1.E-6_wp), dt_zu ) ; !RM 2nd part!!!
         dq_zu(:,:) = q_zu_i - qs_i ; !  dq_zu = SIGN( MAX(ABS(dq_zu),1.E-9_wp), dq_zu )

         !! Now we can get the urbulent scales:
         ztmp0(:,:)  = SQRT(Cd_i(:,:)) ; ! == u*/Ubzu
         u_star = ztmp0(:,:) * Ubzu(:,:)
         ztmp0(:,:)  = 1. / MAX(ztmp0(:,:), 1.E-15) ; ! == Ubzu/u*
         t_star =      Ch_i(:,:)  * dt_zu(:,:) * ztmp0(:,:)
         q_star =      Ce_i(:,:)  * dq_zu(:,:) * ztmp0(:,:)

         !!Inverse of Obukov length (1/L) :
         ztmp0(:,:) = One_on_L(t_zu_i, q_zu_i, u_star, t_star, q_star)  ! 1/L == 1/[Obukhov length]
         ztmp0(:,:) = SIGN( MIN(ABS(ztmp0(:,:)),200._wp), ztmp0(:,:) ) ! (prevents FPE from stupid values from masked region later on...)

         !! Stability parameters "zeta" :
         zeta_u(:,:) = zu*ztmp0(:,:)
         zeta_u(:,:) = SIGN( MIN(ABS(zeta_u(:,:)),50.0_wp), zeta_u(:,:) )
         IF( .NOT. l_zt_equal_zu ) THEN
            zeta_t(:,:) = zt*ztmp0(:,:)
            zeta_t(:,:) = SIGN( MIN(ABS(zeta_t(:,:)),50.0_wp), zeta_t(:,:) )
         END IF

         !! Update C_D:
         ztmp0(:,:) = 1._wp + zsqrtCDN/vkarmn*(zlog2 - psi_m_ice(zeta_u(:,:)))
         Cd_i(:,:)  = MIN( MAX( CdN / ( ztmp0(:,:)*ztmp0(:,:) ), Cx_min ) , 1.9E-3_wp )

         !! Update C_H and C_E
         ztmp0(:,:) = ( zlog2 - psi_h_ice(zeta_u) ) / vkarmn / zsqrtCDN
         ztmp1(:,:) = SQRT(Cd_i(:,:)) / zsqrtCDN
         Ch_i(:,:)  = MIN( MAX( ChN*ztmp1(:,:) / ( 1._wp + ChN*ztmp0(:,:) ) , Cx_min ) , 1.9E-3_wp )
         Ce_i(:,:)  = MIN( MAX( CeN*ztmp1(:,:) / ( 1._wp + CeN*ztmp0(:,:) ) , Cx_min ) , 1.9E-3_wp )

         IF( .NOT. l_zt_equal_zu ) THEN
            !! Re-updating temperature and humidity at zu if zt /= zu :
            ztmp0(:,:)  = psi_h_ice(zeta_u(:,:)) - psi_h_ice(zeta_t(:,:)) + zlog1
            t_zu_i(:,:) =            t_zt(:,:) - t_star(:,:)/vkarmn*ztmp0(:,:)
            q_zu_i(:,:) = MAX(0._wp, q_zt(:,:) - q_star(:,:)/vkarmn*ztmp0(:,:) )
         END IF

      END DO !DO jit = 1, nb_iter


      IF( lreturn_z0 )    xz0 = z0_from_Cd( zu, Cd_i(:,:), ppsi=psi_m_ice(zeta_u(:,:)) )
      IF( lreturn_ustar ) xu_star = u_star
      IF( lreturn_L )     xL      = 1./One_on_L(t_zu_i, q_zu_i, u_star, t_star, q_star)
      IF( lreturn_UN10 )  xUN10   = UN10_from_CD(zu, Ubzu(:,:), Cd_i(:,:), ppsi=psi_m_ice(zeta_u(:,:)) )

      DEALLOCATE ( Ubzu, u_star, t_star, q_star, zeta_u, dt_zu, dq_zu, ztmp0, ztmp1 )
      IF( .NOT. l_zt_equal_zu ) DEALLOCATE( zeta_t )

   END SUBROUTINE turb_ice_easy




   FUNCTION psi_m_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for momentum
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ice
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )            !
            zta = pzeta(ji,jj)
            !
            ! Unstable stratification:
            zx = ABS(1._wp - 16._wp*zta)**.25              !  (16 here, not 15!)

            zpsi_u =      LOG( (1._wp + zx*zx)/2. ) &  ! Eq.(30) Jordan et al. 1999
               &     + 2.*LOG( (1._wp + zx       )/2. ) &
               &    - 2.*ATAN( zx ) + 0.5*rpi

            ! Stable stratification:
            zpsi_s = - ( 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35*zta) + 10.7_wp )  ! Eq.(33) Jordan et al. 1999

            !! Combine:
            zstab = 0.5 + SIGN(0.5_wp, zta)
            psi_m_ice(ji,jj) = (1._wp - zstab) * zpsi_u   & ! Unstable (zta < 0)
               &                   + zstab     * zpsi_s     ! Stable (zta > 0)
            !
      END_2D
   END FUNCTION psi_m_ice


   FUNCTION psi_h_ice( pzeta )
      !!----------------------------------------------------------------------------------
      !! ** Purpose: compute the universal profile stability function for
      !!             temperature and humidity
      !!
      !!
      !!     Andreas et al 2005 == Jordan et al. 1999
      !!
      !!     Psi:
      !!     Unstable => Paulson 1970
      !!     Stable   => Holtslag & De Bruin 1988
      !!
      !!             pzeta : stability paramenter, z/L where z is altitude
      !!                     measurement and L is M-O length
      !!
      !! ** Author: L. Brodeau, 2020 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ice
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zta, zx, zpsi_u, zpsi_s, zstab
      !!----------------------------------------------------------------------------------
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )            !
            zta = pzeta(ji,jj)
            !
            ! Unstable stratification:
            zx = ABS(1._wp - 16._wp*zta)**.25              !  (16 here, not 15!)

            zpsi_u =   2.*LOG( (1._wp + zx*zx)/2. )  ! Eq.(31) Jordan et al. 1999

            ! Stable stratification (identical to Psi_m!):
            zpsi_s = - ( 0.7_wp*zta + 0.75_wp*(zta - 14.3_wp)*EXP( -0.35*zta) + 10.7_wp )  ! Eq.(33) Jordan et al. 1999

            !! Combine:
            zstab = 0.5 + SIGN(0.5_wp, zta)
            psi_h_ice(ji,jj) = (1._wp - zstab) * zpsi_u   & ! Unstable (zta < 0)
               &                   + zstab     * zpsi_s     ! Stable (zta > 0)
            !
      END_2D
   END FUNCTION psi_h_ice

   !!======================================================================
END MODULE sbcblk_algo_ice_easy
