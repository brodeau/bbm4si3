MODULE icedyn_rhg_bbm
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_bbm  ***
   !!   Sea-Ice dynamics : rheology Britle Maxwell X
   !!======================================================================
   !! History :
   !!            4.2  !  2022     (L. Brodeau) `BBM` [starting from `icedyn_rhg_evp`]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_rhg_bbm : computes ice velocities from BBM rheology
   !!   rhg_bbm_rst     : read/write BBM fields in ice restart
   !!----------------------------------------------------------------------
   USE phycst         ! Physical constant
   USE dom_oce        ! Ocean domain
   USE sbc_oce , ONLY : ln_ice_embd, nn_fsbc, ssh_m
   USE sbc_ice , ONLY : utau_ice, vtau_ice, snwice_mass, snwice_mass_b
   USE ice            ! sea-ice: ice variables
   USE icevar         ! ice_var_sshdyn
   USE bdy_oce , ONLY : ln_bdy
   USE bdyice
#if defined key_agrif
   USE agrif_ice_interp
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE prtctl         ! Print control

   USE icedyn_rhg_util

   USE iceistate , ONLY : ln_iceini

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg_bbm_init  ! called by icedyn_rhg.F90
   PUBLIC   ice_dyn_rhg_bbm       ! called by icedyn_rhg.F90
   PUBLIC   rhg_bbm_rst           ! called by icedyn_rhg.F90

   !! Parameters too confidential to be in the namelist:
   REAL(wp), PARAMETER :: &
      &                   rz_nup   = 1._wp/3._wp, &  !: Poisson's ratio
      &                   rz_muMC  = 0.7_wp,      &  !: slope of Mohr-Coulomb enveloppe
      &                   reps6    = 1.e-6_wp,    &
      &                   reps12   = 1.e-12_wp,   &
      &                   reps24   = 1.e-24_wp

   REAL(wp),  SAVE :: rk0  ! factor to stiffness matrix => 1._wp / ( 1._wp - rz_nup*rz_nup)
   REAL(wp),  SAVE :: rk11, rk22, rk12, rk33 ! elements of stiffness matrix
   REAL(wp),  SAVE :: rlambda0, rsqrt_nu_rhoi, rsqrt_E0 ! Constant part of Eq.28

   !! Arrays to be allocated into `ice_dyn_rhg_bbm_init()`:
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Uu_sub, Uv_sub, Vv_sub, Vu_sub !: ice velocities that evolve at sub-time-step
   !!    those that remain constant:
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Xmsk00, Xdxt, Xdxf, Xcohst, Xcohsf, XNlimt, XNlimf
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Xe1t2, Xe2t2, Xe1f2, Xe2f2   !optimization, will avoid these array multiplications countless times...
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xtcoast, xfcoast ! to prevent doing cross-nudging at the coast!
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xCNt, xCNf   ! cross nudging coefficients (time-dependant)

   LOGICAL, PARAMETER :: l_use_v_for_h = .FALSE.

   LOGICAL, SAVE :: l_CN        !: whether cross nudging is used ?
   LOGICAL, SAVE :: l_CN_is_2d  !: whether cross nudging coefficient is a 2D array, not a scalar

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_rhg_bbm.F90 13646 2020-10-20 15:33:01Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg_bbm_init( )
      !!-------------------------------------------------------------------
      !! Called into `ice_dyn_rhg_init()@icedyn_rhg.F90`
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   ierror
      REAL(wp) :: zr
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zt1, zt2, zt3, zt4
      INTEGER :: jm
      !!-------------------------------------------------------------------
      l_CN       = ( rn_crndg > 0._wp )
      l_CN_is_2d = ( ln_boost_CN_coast .OR. ln_boost_CN_high_dmg ) ! cross nudging coefficient will be a 2D array, not a scalar

      xmskt(:,:) =     tmask(:,:,1)
      xmskf(:,:) = MIN(fmask(:,:,1), 1._wp)

      !! Fill the stiffness matrix:
      rk0  = 1._wp / ( 1._wp - rz_nup*rz_nup)
      rk11 = rk0
      rk12 = rk0 * rz_nup
      rk22 = rk0
      rk33 = rk0 * (1._wp - rz_nup)

      rlambda0 = rn_mu0 / rn_E0     !: Viscosity / Elasticity of undamaged ice (aka relaxation time) [s]
      rsqrt_nu_rhoi = SQRT( 2._wp*(1._wp + rz_nup)*rhoi )
      rsqrt_E0      = SQRT( rn_E0 )

      ALLOCATE( zt1(jpi,jpj), zt2(jpi,jpj), utauVice(jpi,jpj), vtauUice(jpi,jpj) ,  STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate `utauVice,vtauUice`' )

      ALLOCATE(    Xdxt(jpi,jpj),  Xdxf(jpi,jpj), &
         &      Xcohst(jpi,jpj), Xcohsf(jpi,jpj), XNlimt(jpi,jpj), XNlimf(jpi,jpj),               &
         &      Xe1t2(jpi,jpj), Xe2t2(jpi,jpj), Xe1f2(jpi,jpj), Xe2f2(jpi,jpj),                   &
         &      Xmsk00(jpi,jpj),     STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate arrays' )

      Xe1t2(:,:) = e1t(:,:) * e1t(:,:) * xmskt(:,:)
      Xe2t2(:,:) = e2t(:,:) * e2t(:,:) * xmskt(:,:)
      Xe1f2(:,:) = e1f(:,:) * e1f(:,:) * xmskf(:,:)
      Xe2f2(:,:) = e2f(:,:) * e2f(:,:) * xmskf(:,:)

      ALLOCATE( Uu_sub(jpi,jpj), Uv_sub(jpi,jpj), Vv_sub(jpi,jpj), Vu_sub(jpi,jpj),       STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate ice velocity arrays' )

      ALLOCATE( uVice(jpi,jpj) , vUice(jpi,jpj) , a_f(jpi,jpj) , h_f(jpi,jpj) , STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate PUBLIC arrays' )

      IF( l_CN ) THEN
         !
         IF( ln_boost_CN_high_dmg ) THEN
            IF(rn_max_CN_dmg<=rn_crndg) CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: `rn_max_CN_dmg` must be > `rn_crndg`' )
         ENDIF
         !
         IF( ln_boost_CN_coast ) THEN
            IF(rn_max_CN_coast<=rn_crndg) CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: `rn_max_CN_coast` must be > `rn_crndg`' )
            ALLOCATE( zt3(jpi,jpj), zt4(jpi,jpj), xtcoast(jpi,jpj),  xfcoast(jpi,jpj) , STAT=ierror )
            IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate `xtcoast` & `xfcoast` arrays' )
            !
            xtcoast(:,:) = xmskt(:,:)
            xfcoast(:,:) = xmskf(:,:)
            zt1(:,:)   = xmskt(:,:)
            zt3(:,:)   = xmskf(:,:)
            DO jm=1, 3
               zr = rn_max_CN_coast - REAL(jm-1) * (rn_max_CN_coast-rn_crndg)/3._wp
               !
               zt2(:,:) = 0._wp
               zt2(2:jpi-1,2:jpj-1) =   zt1(2:jpi-1,3:jpj) + zt1(1:jpi-2,2:jpj-1) + zt1(2:jpi-1,1:jpj-2) + zt1(3:jpi,2:jpj-1) &
                  &                       + zt1(3:jpi,3:jpj)   + zt1(1:jpi-2,3:jpj)   + zt1(1:jpi-2,1:jpj-2) + zt1(3:jpi,1:jpj-2)
               zt2(:,:) = zt2(:,:)/8._wp*xmskt(:,:)
               WHERE( (zt2 < 1._wp).AND.(zt1 > 0._wp) ) xtcoast = zr
               xtcoast(:,:) = xtcoast(:,:)*xmskt(:,:)
               WHERE( (xtcoast > 1.01_wp).OR.(xtcoast < 0.09_wp) ) zt1 = 0.
               !
               zt4(:,:) = 0._wp
               zt4(2:jpi-1,2:jpj-1) =   zt3(2:jpi-1,3:jpj) + zt3(1:jpi-2,2:jpj-1) + zt3(2:jpi-1,1:jpj-2) + zt3(3:jpi,2:jpj-1) &
                  &                       + zt3(3:jpi,3:jpj)   + zt3(1:jpi-2,3:jpj)   + zt3(1:jpi-2,1:jpj-2) + zt3(3:jpi,1:jpj-2)
               zt4(:,:) = zt4(:,:)/8._wp*xmskf(:,:)
               WHERE( (zt4 < 1._wp).AND.(zt3 > 0._wp) ) xfcoast = zr
               xfcoast(:,:) = xfcoast(:,:)*xmskf(:,:)
               WHERE( (xfcoast > 1.01_wp).OR.(xfcoast < 0.09_wp) ) zt3 = 0.
               !
               CALL lbc_lnk( 'icedyn_rhg_bbm', zt1,'T',1._wp, xtcoast,'T',1._wp, zt3,'F',1._wp, xfcoast,'F',1._wp )
            END DO
            DEALLOCATE( zt3, zt4 )

         ENDIF !IF( ln_boost_CN_coast )

         IF( l_CN_is_2d ) THEN
            ALLOCATE( xCNt(jpi,jpj),  xCNf(jpi,jpj) , STAT=ierror )
            IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate `xCNt` & `xCNf` arrays' )
         ENDIF

      ENDIF !IF( l_CN )

      !! Going to use a constant `dx`:
      Xdxt(:,:) = 0._wp ; Xdxf(:,:) = 0._wp
      Xdxt(:,:) = SQRT( e1t(:,:)*e2t(:,:) )  ! Local `dx` of grid cell [m]
      Xdxf(:,:) = SQRT( e1f(:,:)*e2f(:,:) )  ! Local `dx` of grid cell [m]
      zr = SUM(Xdxt(:,:)*xmskt(:,:)) / MAX( SUM(xmskt(:,:)) , reps6 )
      IF( lwp ) WRITE(numout,*) '-- ice_dyn_rhg_bbm_init: average `dx` for BBM => ', REAL(zr/1000._wp,4), ' km'

      !! Cohesion and upper limit for compressive stress:
      zt1(:,:) = SQRT( rn_l_ref/Xdxt(:,:) )
      zt2(:,:) = SQRT( rn_l_ref/Xdxf(:,:) )
      !
      Xcohst(:,:) = rn_c_ref*zt1(:,:)
      Xcohsf(:,:) = rn_c_ref*zt2(:,:)
      !
      XNlimt(:,:) = rn_Nref*zt1(:,:) ! `N` of [Eq.29]
      XNlimf(:,:) = rn_Nref*zt2(:,:) ! `N` of [Eq.29]

      DEALLOCATE( zt1, zt2 )

      IF( lwp ) WRITE(numout,*) ''

   END SUBROUTINE ice_dyn_rhg_bbm_init


   SUBROUTINE ice_dyn_rhg_bbm( kt, Kmm, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i )
      !!-------------------------------------------------------------------
      !!                 ***  SUBROUTINE ice_dyn_rhg_bbm  ***
      !!                             BBM-C-grid
      !!
      !! ** purpose : determines sea ice drift from wind stress, ice-ocean
      !!  stress and sea-surface slope. Ice-ice interaction is described by
      !!  the BBM rheology of Olason et al., 2022.
      !!
      !! ** Inputs  : - wind forcing (stress), oceanic currents
      !!                ice total volume (vt_i) per unit area
      !!                snow total volume (vt_s) per unit area
      !!
      !! ** Action  : - compute u_ice, v_ice : the components of the
      !!                sea-ice velocity vector
      !!              - compute delta_i, shear_i, divu_i, which are inputs
      !!                of the ice thickness distribution
      !!
      !! ** Steps   : 0) compute mask at F point
      !!              1) Compute ice snow mass, ice strength
      !!              2) Compute wind, oceanic stresses, mass terms and
      !!                 coriolis terms of the momentum equation
      !!              3) Solve the momentum equation (iterative procedure)
      !!              4) Recompute delta, shear and divergence
      !!                 (which are inputs of the ITD) & store stress
      !!                 for the next time step
      !!              5) Diagnostics including charge ellipse
      !!
      !! ** Notes   :
      !!
      !!
      !!
      !!
      !! References : Olason et al., 2022 #fixme
      !!              Hunke and Dukowicz, JPO97
      !!              Bouillon et al., Ocean Modelling 2009
      !!              Bouillon et al., Ocean Modelling 2013
      !!              Kimmritz et al., Ocean Modelling 2016 & 2017
      !!-------------------------------------------------------------------
      INTEGER                 , INTENT(in ) ::   kt                                    ! time step
      INTEGER                 , INTENT(in ) ::   Kmm                                   ! ocean time level index
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pstress1_i, pstress2_i, pstress12_i   !
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pshear_i  , pdivu_i   , pdelta_i      !
      !!
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   jter         ! local integers
      !
      REAL(wp) ::   zrhoco                                              ! rho0 * rn_cio
      REAL(wp) ::   zdtbbm, z1_dtbbm                                    ! time step for subcycling
      REAL(wp) ::   zm1, zm2, zm3, zmassU, zmassV                       ! ice/snow mass and volume
      REAL(wp) ::   zds, zds2, zdt, zdt2              ! temporary scalars
      REAL(wp) ::   zTauO, zRHS
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zAu  , zAv                      ! ice fraction on U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   zmU_t, zmV_t                    ! (ice-snow_mass / dt) on U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   vUoce, uVoce                    ! ocean/ice u/v component on V/U points
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zsshdyn                         ! array used for the calculation of ice surface slope:
      !                                                                 !    ocean surface (ssh_m) if ice is not embedded
      !                                                                 !    ice bottom surface if ice is embedded
      REAL(wp), DIMENSION(jpi,jpj) ::   zfUu  , zfVv                    ! internal stresses
      REAL(wp), DIMENSION(jpi,jpj) ::   zfUv  , zfVu                    ! internal stresses (E-grid, #E)
      REAL(wp), DIMENSION(jpi,jpj) ::   zspgUu, zspgVv, zspgUv, zspgVu  ! surface pressure gradient at U/V points
      REAL(wp), DIMENSION(jpi,jpj) ::   ztaux_ai, ztauy_ai              ! ice-atm. stress at U-V points
      REAL(wp), DIMENSION(jpi,jpj) ::   ztauxVai, ztauyUai              ! ice-atm. stress
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk01x, zmsk01y                ! dummy arrays
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk00x, zmsk00y                ! mask for ice presence

      REAL(wp), PARAMETER          ::   zepsi  = 1.0e-20_wp             ! tolerance parameter
      REAL(wp), PARAMETER          ::   zmmin  = 1._wp                  ! ice mass (kg/m2)  below which ice velocity becomes very small
      REAL(wp), PARAMETER          ::   zamin  = 0.001_wp               ! ice concentration below which ice velocity becomes very small
      !! --- diags
      REAL(wp) :: zfac
      !!
      !! Add-ons Brodeau for bbm
      REAL(wp), DIMENSION(jpi,jpj) :: ztmp1, ztmp2, ztmp3, ztmp4, zAt, zAf, zht, zhf
      REAL(wp) :: zh, zt1, zt2, ztame, zmacc
      REAL(wp) :: zCorio, zU, zV, zMdt, zUo, zVo
      !!-------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_rhg_bbm: BBM sea-ice rheology'

      ! For diagnostics:
      Xmsk00(:,:) = 0._wp
      WHERE( at_i(:,:) > epsi06 ) Xmsk00(:,:) = 1._wp

      !------------------------------------------------------------------------------!
      ! 1) define some variables and initialize arrays
      !------------------------------------------------------------------------------!

      !! Ice concentration and thicknes @T we are going to work with:
      zAt(:,:) =      at_i(:,:)
      IF( l_use_v_for_h ) THEN
         zht(:,:) = MAX( vt_i(:,:) , 0._wp)
      ELSE
         zht(:,:) = MAX( hm_i(:,:) , 0._wp)
         WHERE( zAt <= 1.E-3_wp ) zht = 0._wp
      ENDIF

      zAf(:,:) = MIN( MAX( rmpT2F( zAt,  lconserv=.TRUE. ) , 0._wp ), rn_amax_n ) * xmskf(:,:) ! Ice conc. at F-points #fixme: add south!
      zhf(:,:) =      MAX( rmpT2F( zht,  lconserv=.TRUE. ) , 0._wp )              * xmskf(:,:) ! Ice thickness at F-points
      IF( .NOT. l_use_v_for_h ) THEN
         WHERE( zAf <= 1.E-3_wp ) zhf = 0._wp
      END IF
      CALL lbc_lnk( 'icedyn_rhg_bbm',  zAf,'F',1._wp,  zhf,'F',1._wp )

      !! Compute h_f... Useless I think... Unless for saving it as output
      a_f(:,:) = zAf(:,:) ! => used for advection of `dmgf` !
      h_f(:,:) = zhf(:,:)

      zrhoco   = rho0 * rn_cio
      zdtbbm   = rdt_ice / REAL( nn_nbbm, wp )
      z1_dtbbm = 1._wp / zdtbbm

      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      zsshdyn(:,:) = ice_var_sshdyn( ssh_m, snwice_mass, snwice_mass_b)

      zspgUu(:,:) = 0._wp ; zspgUv(:,:) = 0._wp
      zspgVv(:,:) = 0._wp ; zspgVu(:,:) = 0._wp

      !! Taming factor for wind-stress at initialization:
      ztame = 1._wp
      IF( ln_tame_ini_ws .AND. (.NOT. ln_rstart) ) THEN
         zh = REAL(kt-nit000,wp) * rdt_ice / 3600._wp   ! nb of hours since initialization
         ztame = 1._wp / ( 1._wp + EXP(-0.25_wp*(zh - rn_half_tame)) )
      ENDIF

      ! Ocean currents at U-V points:
      uVoce(:,:) = rmpU2V( u_oce )
      vUoce(:,:) = rmpV2U( v_oce )

      ztmp1(:,:) = rmpT2F( zsshdyn ) ! SSH at F-points

      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )

            ! ice fraction at U-V points
            zAu(ji,jj) = 0.5_wp * ( zAt(ji,jj) * e1e2t(ji,jj) + zAt(ji+1,jj) * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
            zAv(ji,jj) = 0.5_wp * ( zAt(ji,jj) * e1e2t(ji,jj) + zAt(ji,jj+1) * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)

            ! Ice/snow mass at U-V points
            zm1 = ( rhos * vt_s(ji  ,jj  ) + rhoi * vt_i(ji  ,jj  ) )
            zm2 = ( rhos * vt_s(ji+1,jj  ) + rhoi * vt_i(ji+1,jj  ) )
            zm3 = ( rhos * vt_s(ji  ,jj+1) + rhoi * vt_i(ji  ,jj+1) )
            zmassU = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm2 * e1e2t(ji+1,jj) ) * r1_e1e2u(ji,jj) * umask(ji,jj,1)
            zmassV = 0.5_wp * ( zm1 * e1e2t(ji,jj) + zm3 * e1e2t(ji,jj+1) ) * r1_e1e2v(ji,jj) * vmask(ji,jj,1)

            ! m/dt
            zmU_t(ji,jj)    = zmassU * z1_dtbbm
            zmV_t(ji,jj)    = zmassV * z1_dtbbm

            ! Surface pressure gradient (- m*g*GRAD(ssh)) at U-V points
            zspgUu(ji,jj) = - zmassU * grav * ( zsshdyn(ji+1,jj) - zsshdyn(ji,jj) ) * r1_e1u(ji,jj)
            zspgVv(ji,jj) = - zmassV * grav * ( zsshdyn(ji,jj+1) - zsshdyn(ji,jj) ) * r1_e2v(ji,jj)
            zspgUv(ji,jj) = - zmassV * grav * (   ztmp1(ji,jj)   - ztmp1(ji-1,jj) ) * r1_e1v(ji,jj)  ! `ztmp1` is `zsshdyn` interpolated  @F !
            zspgVu(ji,jj) = - zmassU * grav * (   ztmp1(ji,jj)   - ztmp1(ji,jj-1) ) * r1_e2u(ji,jj)  ! `ztmp1` is `zsshdyn` interpolated  @F !

            ! masks
            zmsk00x(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassU ) )  ! 0 if no ice
            zmsk00y(ji,jj) = 1._wp - MAX( 0._wp, SIGN( 1._wp, -zmassV ) )  ! 0 if no ice

            ! switches
            IF( zmassU <= zmmin .AND. zAu(ji,jj) <= zamin ) THEN
               zmsk01x(ji,jj) = 0._wp
            ELSE
               zmsk01x(ji,jj) = 1._wp
            ENDIF
            IF( zmassV <= zmmin .AND. zAv(ji,jj) <= zamin ) THEN
               zmsk01y(ji,jj) = 0._wp
            ELSE
               zmsk01y(ji,jj) = 1._wp
            ENDIF

      END_2D

      ! Drag ice-atm. in both T- and F-centric contextes:
      ztaux_ai(:,:) = ztame * zAu(:,:) * utau_ice(:,:) * umask(:,:,1)
      ztauy_ai(:,:) = ztame * zAv(:,:) * vtau_ice(:,:) * vmask(:,:,1)
      ztauxVai(:,:) = ztame * zAv(:,:) * utauVice(:,:) * vmask(:,:,1)
      ztauyUai(:,:) = ztame * zAu(:,:) * vtauUice(:,:) * umask(:,:,1)
      
      IF( l_CN_is_2d ) THEN
         !! Preparing 2D version of cross-nudging coefficients:
         IF( ln_boost_CN_high_dmg ) THEN
            !
            ztmp1 = rmpF2T( dmgf, lconserv=.TRUE. ) !! Averaged `dmgf` at T-points:
            ztmp2 = rmpT2F( dmgt, lconserv=.TRUE. ) !! Averaged `dmgt` at F-points:
            !Checked it was not needed: CALL lbc_lnk( 'icedyn_rhg_bbm', ztmp1,'T',1._wp,  ztmp2,'F',1._wp )
            !
            ztmp1(:,:) = ( 0.25_wp * dmgt(:,:)  +  0.75_wp * ztmp1(:,:) ) * xmskt(:,:) ! adding a contribution from original T-point
            ztmp2(:,:) = ( 0.25_wp * dmgf(:,:)  +  0.75_wp * ztmp2(:,:) ) * xmskf(:,:) ! adding a contribution from original F-point
            !
            ztmp3(:,:) = MAX( rn_max_CN_dmg*EXP(rn_C0*(1._wp - ztmp1(:,:))) , rn_crndg ) * xmskt(:,:)
            ztmp4(:,:) = MAX( rn_max_CN_dmg*EXP(rn_C0*(1._wp - ztmp2(:,:))) , rn_crndg ) * xmskf(:,:)
            !
         ELSE
            ztmp3(:,:) = rn_crndg * xmskt(:,:)
            ztmp4(:,:) = rn_crndg * xmskf(:,:)
         ENDIF !IF( ln_boost_CN_high_dmg )
         !
         IF(ln_boost_CN_coast) THEN
            !! Apply boost at the coast:
            ztmp3(:,:) = MAX( ztmp3(:,:), xtcoast(:,:) )
            ztmp4(:,:) = MAX( ztmp4(:,:), xfcoast(:,:) )
         ENDIF
         !
         IF( iom_use('cncoeff_t') ) CALL iom_put( 'cncoeff_t' , ztmp3(:,:) )
         IF( iom_use('cncoeff_f') ) CALL iom_put( 'cncoeff_f' , ztmp4(:,:) )
         !
         xCNt(:,:) = zdtbbm / rdt_ice * ztmp3(:,:)
         xCNf(:,:) = zdtbbm / rdt_ice * ztmp4(:,:)
         !
      ENDIF !IF( l_CN_is_2d )

      !! Because there has been advection of `damage` and stress tensors since last time we called `clean_small_a_all`:
      CALL clean_small_a_all( zAt, zAf,  dmgt, dmgf,  sgm11t, sgm22t, sgm12t,  sgm11f, sgm22f, sgm12f )

      ! Going to average (set to 0 before accumulating during the `nn_nbbm` sub time steps):
      u_ice(:,:) = 0._wp
      uVice(:,:) = 0._wp
      v_ice(:,:) = 0._wp
      vUice(:,:) = 0._wp
      zmacc      = 1._wp/REAL(nn_nbbm)

      !                                               ! ==================== !
      DO jter = 1 , nn_nbbm                           !    loop over jter    !
         !                                            ! ==================== !

         ! ---  Updates the components of the internal stress tensor and the damage in both T- & F-centric worlds ---
         CALL update_stress_dmg( kt, jter, zdtbbm, Uu_sub, Vv_sub, Uv_sub, Vu_sub, zAt, zAf, zht, zhf, & !
            &                                      sgm11t, sgm22t, sgm12f, sgm11f, sgm22f, sgm12t, dmgt, dmgf )


         !! Terms of the divergence of the stress tensor !
         !! ==============================================
         ! --- Ice internal stresses (Appendix C of Hunke and Dukowicz, 2002) --- !
         !!     Stresses used in the following must be vertically-integrated stresses (sigma*h): [Pa.m],
         !!     so the gradient calculated here gives something in [Pa] !!!

         ztmp1(:,:) = sgm11t(:,:) * zht(:,:)
         ztmp2(:,:) = sgm22t(:,:) * zht(:,:)
         ztmp3(:,:) = sgm12f(:,:) * zhf(:,:)
         !!
         CALL div_stress_tensor( 'T',  Xe1t2, Xe2t2,  Xe1f2, Xe2f2,  r1_e2u, r1_e1u, r1_e1v, r1_e2v,  r1_e1e2u, r1_e1e2v,  &
            &                          ztmp1, ztmp2, ztmp3,  zfUu, zfVv )

         ztmp1(:,:) = sgm11f(:,:) * zhf(:,:)
         ztmp2(:,:) = sgm22f(:,:) * zhf(:,:)
         ztmp3(:,:) = sgm12t(:,:) * zht(:,:)
         !!
         CALL div_stress_tensor( 'F',  Xe1f2, Xe2f2,  Xe1t2, Xe2t2,  r1_e2v, r1_e1v, r1_e1u, r1_e2u,  r1_e1e2v, r1_e1e2u,  &
            &                          ztmp1, ztmp2, ztmp3,  zfUv, zfVu )


         ! --- Computation of ice velocity --- !
         IF( MOD(jter,2) == 0 ) THEN
            !! Update `Vv_sub` at V-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtVv.h90"
            !!
            !! Update `Vu_sub` at U-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtVu.h90"
            !!
            !! Update `Uu_sub` at U-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtUu.h90"
            !!
            !! Update `Uv_sub` at V-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtUv.h90"
         ELSE
            !! Update `Uu_sub` at U-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtUu.h90"
            !!
            !! Update `Uv_sub` at V-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtUv.h90"
            !!
            !! Update `Vv_sub` at V-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtVv.h90"
            !!
            !! Update `Vu_sub` at U-points
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~
#           include "icedyn_rhg_bbm_updtVu.h90"

         ENDIF

         !#fixme:
         !! => not needed (double-checked!) it is done through mass stuff...
         !!   => but the question about the implied no-slip condition in F-centric world remains !
         !Uu_sub(:,:) = Uu_sub(:,:) * umask(:,:,1)
         !Vv_sub(:,:) = Vv_sub(:,:) * vmask(:,:,1)
         !Uv_sub(:,:) = Uv_sub(:,:) * vmask(:,:,1) ; ! This applies a no-slip condition !!! #fixme!?
         !Vu_sub(:,:) = Vu_sub(:,:) * umask(:,:,1) ; ! This applies a no-slip condition !!! #fixme!?

         CALL lbc_lnk( 'icedyn_rhg_bbm', Uu_sub,'U',-1._wp, Vv_sub,'V',-1._wp, &
            &                            Uv_sub,'V',-1._wp, Vu_sub,'U',-1._wp )

         u_ice(:,:) = u_ice(:,:) + zmacc*Uu_sub(:,:)
         v_ice(:,:) = v_ice(:,:) + zmacc*Vv_sub(:,:)
         uVice(:,:) = uVice(:,:) + zmacc*Uv_sub(:,:)
         vUice(:,:) = vUice(:,:) + zmacc*Vu_sub(:,:)
         !                                                ! ==================== !
      END DO                                              !  end loop over jter  !
      !                                                   ! ==================== !

      !! Closer look at the components of rate-of-strain tensor:
      IF( iom_use('e11t').OR.iom_use('e22t').OR.iom_use('e12t') ) THEN
         CALL strain_rate( 'T', u_ice, v_ice, uVice, vUice, &
            &              r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, Xe1t2, Xe2t2, tmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3, lblnk=.FALSE. )
         IF( iom_use('e11t') ) CALL iom_put( 'e11t' , ztmp1 )
         IF( iom_use('e22t') ) CALL iom_put( 'e22t' , ztmp2 )
         IF( iom_use('e12t') ) CALL iom_put( 'e12t' , ztmp3 )
      ENDIF
      IF( iom_use('e11f').OR.iom_use('e22f').OR.iom_use('e12f') ) THEN
         CALL strain_rate( 'F', uVice, vUice, u_ice, v_ice, &
            &              r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3,  lblnk=.FALSE. )
         IF( iom_use('e11f') ) CALL iom_put( 'e11f' , ztmp1 )
         IF( iom_use('e22f') ) CALL iom_put( 'e22f' , ztmp2 )
         IF( iom_use('e12f') ) CALL iom_put( 'e12f' , ztmp3 )
      END IF

      !! Saving stress tensor components in the proper units! i.e. Pa :
      IF( iom_use('ice_sig11')  ) CALL iom_put( 'ice_sig11' ,  sgm11t )
      IF( iom_use('ice_sig22')  ) CALL iom_put( 'ice_sig22' ,  sgm22t )
      IF( iom_use('ice_sig12')  ) CALL iom_put( 'ice_sig12' ,  sgm12f )
      IF( iom_use('ice_sig11f') ) CALL iom_put( 'ice_sig11f',  sgm11f )
      IF( iom_use('ice_sig22f') ) CALL iom_put( 'ice_sig22f',  sgm22f )
      IF( iom_use('ice_sig12t') ) CALL iom_put( 'ice_sig12t',  sgm12t )
      !
      IF( iom_use('normstr') ) CALL iom_put( 'normstr' , 0.5_wp*(sgm11t+sgm22t) ) ! First invariant of stress tensor @T
      IF( iom_use('normstrf')) CALL iom_put( 'normstrf', 0.5_wp*(sgm11f+sgm22f) ) ! First invariant of stress tensor @F
      IF( iom_use('sheastr')  ) THEN
         ztmp3(:,:) = 0.5_wp * (sgm11t(:,:) - sgm22t(:,:))
         CALL iom_put( 'sheastr',  SQRT(ztmp3(:,:)*ztmp3(:,:) + sgm12t(:,:)*sgm12t(:,:) ) )    ! Second invariant of stress tensor @T
      ENDIF
      IF( iom_use('sheastrf')  ) THEN
         ztmp3(:,:) = 0.5_wp * (sgm11f(:,:) - sgm22f(:,:))
         CALL iom_put( 'sheastrf', SQRT(ztmp3(:,:)*ztmp3(:,:) + sgm12f(:,:)*sgm12f(:,:) ) )    ! Second invariant of stress tensor @F
      ENDIF

      IF(iom_use('zfUu'))  CALL iom_put( 'zfUu' , zfUu )
      IF(iom_use('zfVv'))  CALL iom_put( 'zfVv' , zfVv )
      IF(iom_use('zfUv'))  CALL iom_put( 'zfUv' , zfUv )
      IF(iom_use('zfVu'))  CALL iom_put( 'zfVu' , zfVu )


      !------------------------------------------------------------------------------!
      ! 4) Recompute delta, shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!

      pdelta_i(:,:) = 0._wp ; ! LB> #fixme: is it really not needed in any context?

      CALL strain_rate( 'T', u_ice, v_ice, uVice, vUice, &
         &              r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, Xe1t2, Xe2t2, tmask(:,:,1), &
         &              ztmp1, ztmp2, ztmp3, lblnk=.TRUE., pdiv=pdivu_i, pmaxshr=pshear_i )
      ! --- divergence of velocity field @T:
      IF( iom_use('icedivt') )  CALL iom_put( 'icedivt' , pdivu_i )
      ! --- shear of velocity field @T:
      IF( iom_use('iceshrt') )  CALL iom_put( 'iceshrt' , ztmp3 )
      ! --- MAXIMUM shear of velocity field @T:
      IF( iom_use('iceshet') )  CALL iom_put( 'iceshet' , pshear_i )
      ! --- total deformation of velocity field @T:
      IF( iom_use('icedeft') ) CALL iom_put( 'icedeft', SQRT( pshear_i*pshear_i + pdivu_i*pdivu_i ) )

      IF( iom_use('icedivf') .OR. iom_use('iceshrf') .OR. iom_use('iceshef') .OR. iom_use('icedeff') ) THEN
         CALL strain_rate( 'F', uVice, vUice, u_ice, v_ice, &
            &              r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3,  lblnk=.TRUE., pdiv=zfUu, pmaxshr=zfVv )
         ! --- divergence of velocity field @F:
         IF( iom_use('icedivf') ) CALL iom_put( 'icedivf' , zfUu )
         ! --- shear of velocity field @F:
         IF( iom_use('iceshrf') ) CALL iom_put( 'iceshrf' , ztmp3 )
         ! --- MAXIMUM shear of velocity field @F:
         IF( iom_use('iceshef') ) CALL iom_put( 'iceshef' , zfVv )
         ! --- total deformation of velocity field @F:
         IF( iom_use('icedeff') ) CALL iom_put( 'icedeff', SQRT( zfVv*zfVv + zfUu*zfUu ) )

      ENDIF



      !------------------------------------------------------------------------------!
      ! 5) diagnostics
      !------------------------------------------------------------------------------!

      ! --- vorticity of velocity field @F:
      IF( iom_use('icevorf') ) THEN
         DO_2D( 0, 0, 0, 0 )
               ztmp3(ji,jj) = (   ( v_ice(ji+1,jj)*r1_e2v(ji+1,jj) - v_ice(ji,jj)*r1_e2v(ji,jj) ) * Xe2f2(ji,jj) &
                  &             - ( u_ice(ji,jj+1)*r1_e1u(ji,jj+1) - u_ice(ji,jj)*r1_e1u(ji,jj) ) * Xe1f2(ji,jj) &
                  &           ) * r1_e1e2f(ji,jj) * fmask(ji,jj,1)   !#fixme: sure about `fmask` here ?
         END_2D
         CALL iom_put( 'icevorf' , ztmp3 )
      ENDIF
      ! --- vorticity of velocity field @T:
      IF( iom_use('icevort') ) THEN
         DO_2D( 0, 0, 0, 0 )
               ztmp3(ji,jj) = (   ( vUice(ji,jj)*r1_e2u(ji,jj) - vUice(ji-1,jj)*r1_e2u(ji-1,jj) ) * Xe2t2(ji,jj) &
                  &             - ( uVice(ji,jj)*r1_e1v(ji,jj) - uVice(ji,jj-1)*r1_e1v(ji,jj-1) ) * Xe1t2(ji,jj) &
                  &           ) * r1_e1e2t(ji,jj) * tmask(ji,jj,1)
         END_2D
         CALL iom_put( 'icevort' , ztmp3 )
      ENDIF

      !     => SI3 expects vertically-integrated stresses in [Pa.m] as output of this routine? (not really used anyway...)
      pstress1_i (:,:) = ( sgm11t(:,:) + sgm22t(:,:) ) * zht(:,:)
      pstress2_i (:,:) = ( sgm11t(:,:) - sgm22t(:,:) ) * zht(:,:)
      pstress12_i(:,:) =          sgm12f(:,:)          * zhf(:,:)

      ! --- ice-atm. stress:
      IF( iom_use('utau_ai') .OR. iom_use('vtau_ai') ) THEN
         CALL iom_put( 'utau_ai' , ztaux_ai )
         CALL iom_put( 'vtau_ai' , ztauy_ai )
      ENDIF
      IF( iom_use('taum_ai') ) THEN
         ztmp1(2:jpi,:) = 0.5_wp * ( ztaux_ai(2:jpi,:) + ztaux_ai(1:jpi-1,:) )
         ztmp2(:,2:jpj) = 0.5_wp * ( ztauy_ai(:,2:jpj) + ztauy_ai(:,1:jpj-1) )
         CALL iom_put( 'taum_ai' , SQRT(ztmp1*ztmp1 + ztmp2*ztmp2) * zAt )
      END IF

      IF( iom_use('utauVai') .OR. iom_use('vtauUai') ) THEN
         CALL iom_put( 'utauVai' , ztauxVai )
         CALL iom_put( 'vtauUai' , ztauyUai )
      ENDIF
      IF( iom_use('taumFai') ) CALL iom_put( 'taumFai' , SQRT(ztauxVai*ztauxVai + ztauyUai*ztauyUai) )  !#fixme: ugly!!!

      ! --- ice-ocean stress:
      IF( iom_use('taum_oi') .OR. iom_use('utau_oi') .OR. iom_use('vtau_oi') ) THEN
         ztmp1(:,:) = umask(:,:,1)
         WHERE( zAu(:,:) < 0.01_wp ) ztmp1(:,:) = 0._wp
         ztmp3(:,:) = u_ice(:,:) - u_oce(:,:) ! dU @U
         ztmp4(:,:) = vUice(:,:) - vUoce(:,:) ! dV @U
         ztmp2(:,:) = SQRT( ztmp3(:,:)*ztmp3(:,:) + ztmp4(:,:)*ztmp4(:,:) ) ! module of relative current at U-point
         zfUu (:,:) = zrhoco * ztmp2(:,:) * ( u_oce(:,:) - u_ice(:,:) )
         IF( iom_use('utau_oi') ) CALL iom_put( 'utau_oi' , zfUu(:,:) * ztmp1 )  ! oce-ice stress /x @ U. MIND: x A !!!
         !!
         ztmp1(:,:) = vmask(:,:,1)
         WHERE( zAv(:,:) < 0.01_wp ) ztmp1(:,:) = 0._wp
         ztmp3(:,:) = v_ice(:,:) - v_oce(:,:) ! dV @V
         ztmp4(:,:) = uVice(:,:) - uVoce(:,:) ! dU @V
         ztmp2(:,:) = SQRT( ztmp3(:,:)*ztmp3(:,:) + ztmp4(:,:)*ztmp4(:,:) ) ! module of relative current at V-point
         zfVv (:,:) = zrhoco * ztmp2(:,:) * ( v_oce(:,:) - v_ice(:,:) )
         IF( iom_use('vtau_oi') ) CALL iom_put( 'vtau_oi' , zfVv(:,:) * ztmp1 )  ! oce-ice stress /y @ V MIND: x A !!!
         !
         IF( iom_use('taum_oi') ) THEN
            ztmp3(:,:) = xmskt(:,:)
            WHERE( zAt(:,:) < 0.01_wp ) ztmp3(:,:) = 0._wp
            ztmp1(2:jpi,:) = 0.5_wp * ( zfUu(2:jpi,:) + zfUu(1:jpi-1,:) )
            ztmp2(:,2:jpj) = 0.5_wp * ( zfVv(:,2:jpj) + zfVv(:,1:jpj-1) )
            CALL iom_put( 'taum_oi' , SQRT(ztmp1*ztmp1 + ztmp2*ztmp2) * ztmp3 )
         END IF
         !
      ENDIF
      !
   END SUBROUTINE ice_dyn_rhg_bbm




   SUBROUTINE update_stress_dmg( kt, kts, pdt, pUu, pVv, pUv, pVu, pAt, pAf, pht, phf,  &
      &                                 ps11t, ps22t, ps12f, ps11f, ps22f, ps12t, pdmgt, pdmgf )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE UPDATE_STRESS_DMG  ***
      !! ** Purpose :
      !!
      !! ** Method  :
      !!
      !! ** Note    : Called at the sub-time-stepping level!
      !!
      !! ** Author : L. Brodeau, 2022
      !!----------------------------------------------------------------------
      INTEGER,                  INTENT(in)    :: kt, kts        ! # of current big and small/sub-time step
      REAL(wp),                 INTENT(in)    :: pdt             ! (small) time-step [s]
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pUu, pVv        ! Ice velocity vector @U & @V
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pUv, pVu        ! Ice velocity vector @V & @U (E-grid)
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pAt, pAf        ! Ice concentration @T & @F
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pht, phf        ! Ice thickness @T & @F
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11t, ps22t, ps12f  ! T-centric Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11f, ps22f, ps12t  ! F-centric Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pdmgt, pdmgf    ! ice damage @T & @F
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: ztp0, ztp1, ztp2, ztp3
      REAL(wp), DIMENSION(jpi,jpj) :: zelat, zelaf, zmulx
      REAL(wp)                     :: ztmp, zmlt
      REAL(wp)                     :: ze11, ze22, ze12
      INTEGER  :: ji, jj
      LOGICAL, PARAMETER :: llbclnk=.FALSE.
      !!----------------------------------------------------------------------

      !! --- First (predictor) update of stress tensor terms [Eq.32] ---
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !! Compute `elasticity`, and `update multiplicator` @ T-points (=> zelat, zmulx)
      CALL PHASE_I( pdt, pAt, pht, pdmgt, ps11t, ps22t,  zelat, zmulx )
      !
      !! Compute the 3 components of the strain-rate tensor @ T-points (=> ztp1, ztp2, ztp3):
      CALL strain_rate( 'T', pUu, pVv, pUv, pVu, r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, Xe1t2, Xe2t2, tmask(:,:,1), &
         &              ztp1, ztp2, ztp3, lblnk=.FALSE. ) !: double-checked that lbc_lnk-ing is NOT NEEDED !!!
      !
      !! Predictor update of the 3 stress tensor components @ T-points:
      DO_2D( nn_hls-1, nn_hls, nn_hls-1, nn_hls )
            ze11 = ztp1(ji,jj)
            ze22 = ztp2(ji,jj)
            ze12 = ztp3(ji,jj)
            zmlt = zmulx(ji,jj) * xmskt(ji,jj)
            ztmp = zelat(ji,jj) * pdt
            !
            ps11t(ji,jj) = zmlt * ( ztmp * ( rk11*ze11 + rk12*ze22) + ps11t(ji,jj) )
            ps22t(ji,jj) = zmlt * ( ztmp * ( rk12*ze11 + rk22*ze22) + ps22t(ji,jj) )
            ps12t(ji,jj) = zmlt * ( ztmp *        rk33 * ze12       + ps12t(ji,jj) )
      END_2D

      !! Compute `elasticity`, and `update multiplicator` @ F-points (=> zelat, zmulx)
      CALL PHASE_I( pdt, pAf, phf, pdmgf, ps11f, ps22f,  zelaf, zmulx ) ! compute `elasticity` and `update multiplicator` @ F
      !
      !! Compute the 3 components of the strain-rate tensor @ F-points (=> ztp1, ztp2, ztp3):
      CALL strain_rate( 'F', pUv, pVu, pUu, pVv, r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
         &              ztp1, ztp2, ztp3, lblnk=.FALSE. ) !: double-checked that lbc_lnk-ing is NOT NEEDED !!!

      !! Predictor update of the 3 stress tensor components @ F-points:
      DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
            ze11 = ztp1(ji,jj)
            ze22 = ztp2(ji,jj)
            ze12 = ztp3(ji,jj)
            zmlt = zmulx(ji,jj) * xmskf(ji,jj)
            ztmp = zelaf(ji,jj) * pdt
            !
            ps11f(ji,jj) = zmlt * ( ztmp * ( rk11*ze11 + rk12*ze22) + ps11f(ji,jj) )
            ps22f(ji,jj) = zmlt * ( ztmp * ( rk12*ze11 + rk22*ze22) + ps22f(ji,jj) )
            ps12f(ji,jj) = zmlt * ( ztmp *        rk33 * ze12       + ps12f(ji,jj) )
      END_2D

      IF( l_CN ) THEN
         !! Cross-nudging on stress tensor components
         !!  Alternate / even/odd sub iterations => gives slightly better results than doing
         !!  it everytime (with half zrbal) on all stress comp.
         !!     (=> double checked that no `lbc_lnk` is needed in `rmpY2X()` !!!)
         CALL CROSS_NUDGING( kts, pdt, pAt, pAf,  ps11t, ps22t, ps12t, ps11f, ps22f, ps12f )
      ELSE
         IF( (kt==1).AND.(lwp) ) WRITE(numout,*) ' *** MIND: no cross-nudging applied!'
      ENDIF

      !! --- Mohr-Coulomb test and britle update ---
      CALL MOHR_COULOMB_DMG( pdt, Xcohst, XNlimt, zelat, Xdxt, ps11t, ps22t, ps12t, pdmgt )
      CALL MOHR_COULOMB_DMG( pdt, Xcohsf, XNlimf, zelaf, Xdxf, ps11f, ps22f, ps12f, pdmgf )

      ! --- Healing of damage with time [Eq.30, Olason et al., 2022]
      !      ==> can suposedely be taken out of sub-time-stepping loop...
      ztp3(:,:) = (rcnd_i*vt_s(:,:))/( rn_cnd_s*MAX(pht(:,:),reps6) ) * xmskt(:,:)   ! => `C` of the `dtemp/(1 + C)` in neXtSIM
      IF(ln_icethd) THEN
         !! Thermo is used, normal stuff:
         ztp1(:,:) = (t_bo(:,:) - t_su(:,:,1)) / (1._wp + ztp3(:,:) ) * xmskt(:,:) ! temperature difference between bottom and surface
      ELSE
         !! Thermo is off, yet we want som refreezing!
         ztp1(:,:) = (-1.8_wp + 25._wp)/ (1._wp + ztp3(:,:) ) * xmskt(:,:) ! temperature difference between bottom and surface
      END IF
      ztp2(:,:) = MIN( ztp1(:,:) / rn_kth , 1._wp/rdt_ice )  ! 1/T_relax
      WHERE( ztp1(:,:) > 0._wp ) pdmgt(:,:) = MAX( pdmgt(:,:) - pdt*ztp2(:,:)      , 0._wp )
      ztp0(:,:) = rmpT2F( ztp1 )
      WHERE( ztp0(:,:) > 0._wp ) pdmgf(:,:) = MAX( pdmgf(:,:) - pdt*rmpT2F( ztp2 ) , 0._wp )

      CALL clean_small_a_all( pAt, pAf,  pdmgt, pdmgf,  ps11t, ps22t, ps12t,  ps11f, ps22f, ps12f )

      !! --- Final LBC linking ---
      CALL lbc_lnk( 'UPDATE_STRESS_DMG@icedyn_rhg_bbm', ps11t,'T',1._wp, ps22t,'T',1._wp, ps12f,'F',1._wp, pdmgt,'T',1._wp, &
         &                                              ps11f,'F',1._wp, ps22f,'F',1._wp, ps12t,'T',1._wp, pdmgf,'F',1._wp  )

   END SUBROUTINE UPDATE_STRESS_DMG



   SUBROUTINE rhg_bbm_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rhg_bbm_rst  ***
      !!
      !! ** Purpose :   Read or write RHG file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER  ::   iter            ! local integer
      INTEGER  ::   id01, id02
      INTEGER  ::   id1, id2, id3, id4, id5, id6, id7, id8
      INTEGER  ::   id11, id12, id13, id14
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialize
         !                                   ! ---------------
         IF( ln_rstart ) THEN                   !* Read the restart file
            !
            id01 = iom_varid( numrir, 'dmgt' , ldstop = .FALSE. )
            id02 = iom_varid( numrir, 'dmgf' , ldstop = .FALSE. )
            !
            id1 = iom_varid( numrir, 'sgm11t' , ldstop = .FALSE. )
            id2 = iom_varid( numrir, 'sgm22t' , ldstop = .FALSE. )
            id3 = iom_varid( numrir, 'sgm12f' , ldstop = .FALSE. )
            !
            id4 = iom_varid( numrir, 'sgm11f' , ldstop = .FALSE. )
            id5 = iom_varid( numrir, 'sgm22f' , ldstop = .FALSE. )
            id6 = iom_varid( numrir, 'sgm12t' , ldstop = .FALSE. )
            !
            id7 = iom_varid( numrir, 'uVice' , ldstop = .FALSE. )
            id8 = iom_varid( numrir, 'vUice' , ldstop = .FALSE. )
            !
            id11 = iom_varid( numrir, 'Uu_sub' , ldstop = .FALSE. )
            id12 = iom_varid( numrir, 'Uv_sub' , ldstop = .FALSE. )
            id13 = iom_varid( numrir, 'Vv_sub' , ldstop = .FALSE. )
            id14 = iom_varid( numrir, 'Vu_sub' , ldstop = .FALSE. )

            IF( MIN( id01, id02 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'dmgt' , dmgt , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'dmgf' , dmgf , cd_type = 'F' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without rheology, set damage @T and @F to 0'
               dmgt(:,:) = 0._wp
               dmgf(:,:) = 0._wp
            ENDIF

            IF( MIN( id1, id2, id3 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'sgm11t' , sgm11t , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22t' , sgm22t , cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12f' , sgm12f , cd_type = 'F' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without rheology, set stresses to 0'
               stress1_i (:,:) = 0._wp
               stress2_i (:,:) = 0._wp
               stress12_i(:,:) = 0._wp
               sgm11t(:,:)     =  0._wp
               sgm22t(:,:)     =  0._wp
               sgm12f(:,:)     =  0._wp

            ENDIF

            IF( MIN( id4, id5, id6 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'sgm11f', sgm11f, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22f', sgm22f, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12t', sgm12t, cd_type = 'T' )
            ELSE                                     ! get them from T-centric stresses...
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without BBM rheology, interpolate F-centric stresses'
               sgm11f(:,:) = rmpT2F( sgm11t, lconserv=.TRUE. )
               sgm22f(:,:) = rmpT2F( sgm22t, lconserv=.TRUE. )
               sgm12t(:,:) = rmpF2T( sgm12f, lconserv=.TRUE. )
               CALL lbc_lnk( 'icedyn_rhg_bbm', sgm11f,'F',1._wp,  sgm22f,'F',1._wp,  sgm12t,'T',1._wp )
            ENDIF

            IF( MIN( id7, id8 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'uVice' , uVice , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'vUice' , vUice , cd_type = 'U', psgn = -1._wp )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without BBM rheology, interpolate F-centric velocities'
               uVice(:,:) = rmpU2V( u_ice )
               vUice(:,:) = rmpV2U( v_ice )
               CALL lbc_lnk( 'rhg_bbm_rst',  uVice,'V',-1._wp, vUice,'U',-1._wp )
            ENDIF

            IF( MIN( id11, id12, id13, id14 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'Uu_sub' , Uu_sub , cd_type = 'U', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Uv_sub' , Uv_sub , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Vv_sub' , Vv_sub , cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numrir, jpdom_auto, 'Vu_sub' , Vu_sub , cd_type = 'U', psgn = -1._wp )
            ELSE
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>>   previous run without BBM rheology, fill sub-ts velocities'
               Uu_sub(:,:) = u_ice(:,:)
               Uv_sub(:,:) = uVice(:,:)
               Vv_sub(:,:) = v_ice(:,:)
               Vu_sub(:,:) = vUice(:,:)
            ENDIF
            !
         ELSE                                   !* Start from rest
            !
            IF(lwp) WRITE(numout,*)
            !
            IF(.NOT. ln_iceini) THEN
               IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set T- and F- centric damage to 0'
               dmgt(:,:) = 0._wp
               dmgf(:,:) = 0._wp
            ELSE
               !PRINT *, 'LOLO[icedyn_rhg_bbm.F90] rhg_bbm_rst(): damage@T taken from `sn_dmg` file => damage@F interpolated!'
               IF(lwp) WRITE(numout,*) '   ==>>>   damage@T taken from `sn_dmg@namini` file => damage@F interpolated!'
               dmgt(:,:) = MIN( MAX(         dmgt(:,:)                , 0._wp ) , rn_dmg_max )
               dmgf(:,:) = MIN( MAX( rmpT2F( dmgt,  lconserv=.TRUE. ) , 0._wp ) , rn_dmg_max )
            ENDIF
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set T-centric stresses to 0'
            sgm11t(:,:) = 0._wp
            sgm22t(:,:) = 0._wp
            sgm12f(:,:) = 0._wp
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set F-centric stresses to 0'
            sgm11f(:,:) = 0._wp
            sgm22f(:,:) = 0._wp
            sgm12t(:,:) = 0._wp
            !
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set F-centric velocities to 0'
            uVice(:,:)  = 0._wp
            vUice(:,:)  = 0._wp
            Uu_sub(:,:) = 0._wp
            Uv_sub(:,:) = 0._wp
            Vv_sub(:,:) = 0._wp
            Vu_sub(:,:) = 0._wp
            !!
         ENDIF
         !
         !
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- rhg-rst ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         CALL iom_rstput( iter, nitrst, numriw, 'dmgt' , dmgt  )
         CALL iom_rstput( iter, nitrst, numriw, 'dmgf' , dmgf  )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sgm11t' , sgm11t )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm22t' , sgm22t )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm12f' , sgm12f )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sgm11f' , sgm11f )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm22f' , sgm22f )
         CALL iom_rstput( iter, nitrst, numriw, 'sgm12t' , sgm12t )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'uVice' , uVice )
         CALL iom_rstput( iter, nitrst, numriw, 'vUice' , vUice )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'Uu_sub' , Uu_sub )
         CALL iom_rstput( iter, nitrst, numriw, 'Uv_sub' , Uv_sub )
         CALL iom_rstput( iter, nitrst, numriw, 'Vv_sub' , Vv_sub )
         CALL iom_rstput( iter, nitrst, numriw, 'Vu_sub' , Vu_sub )
         !
      ENDIF
      !
   END SUBROUTINE rhg_bbm_rst


   SUBROUTINE PHASE_I( pdt, pA, ph, pdmg, ps11, ps22, pelast, pmult )
      !!
      REAL(wp),                 INTENT(in)  :: pdt         ! (small) time-step [s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pA          ! Ice concentration
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ph          ! Ice thickness
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pdmg        ! Ice damage
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11, ps22  ! T-centric Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pelast      ! Ice elasticity
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pmult       ! Multiplicator for stress tensor update
      !!
      REAL(wp) :: zxpC, zsigI, zPmax, zc0, zc1, zPtld, z1md, zlamb
      INTEGER  :: ji, jj

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

            zxpC = EXP( rn_C0*(1._wp - pA(ji,jj)) )  ! `expC` [Eq.8]

            z1md = 1._wp - pdmg(ji,jj)

            pelast(ji,jj) = rn_E0 * z1md * zxpC            ! Elasticity [Eq.9]

            !! --- Lambda ( viscosity / elasticity ) ---
            !!
            !! Normal / V. Danserau:
            !! *    E = E0 * (1 - d)    * exp[-C*(1-A)]    ! elasticity
            !! *    V = V0 * (1 - d)**a * exp[-C*(1-A)]    ! viscosity
            !! *    L = V/E
            !! * => L = L0 * (1 - d)**[a -1]     with L0 == V0/E0
            !zlamb = rlambda0 * z1md **(rz_alrlx - 1._wp) ! Viscous relaxation time aka `Lambda` [Eq.10/Eq.9]   #vero
            !!
            !! Tweaked Olason et al. 2022 / G. Boutin:
            !! *    E = E0 * (1 - d)    * exp[  -C*(1-A)]    ! elasticity
            !! *    V = V0 * (1 - d)**a * exp[-a*C*(1-A)]    ! viscosity
            !! * => L = L0 * (1 - d)**[a -1] * ( exp[  -C*(1-A)] )**[a -1]    with L0 == V0/E0
            !! * => L = L0 * [ (1 - d) * exp[  -C*(1-A)] ]**[a -1]    with L0 == V0/E0
            zlamb = rlambda0 * (z1md * zxpC)**(rn_btrlx - 1._wp) ! Viscous relaxation time aka `Lambda` [Eq.10/Eq.9] #einar

            !! --- P~ / Plastic failure [Eq.8]---
            zsigI  = 0.5_wp * ( ps11(ji,jj) + ps22(ji,jj) ) ! sigI: normal stress aka first invariant
            !!
            zPmax = -rn_P0 * ph(ji,jj)**1.5_wp * zxpC       ! `-P_max` (for sigI<0)
            !zc0 = 0.5_wp + SIGN( 0.5_wp, -zsigI )           ! => if sigI>  0  : zc0=0 else: zc0=1
            !zc1 = 0.5_wp + SIGN( 0.5_wp, zPmax - zsigI )    ! => if sigI<-Pmax: zc1=1 else: zc1=0
            !!
            !zPtld = zc0 * (                                      &  ! if sigI>0 => P~=0 !!!
            !   &              1._wp - zc1                        &  ! Default:    P~ = 1.     when `-Pmax < sigI < 0`
            !   &            + zc1 * zPmax / MIN(zsigI,-reps6) )  ! sigI<-Pmax: P~ = -Pmax / sigI
            zPtld = 0._wp                               ! `sigma_n` > 0
            IF( zsigI < -reps24 ) zPtld = MIN( zPmax/zsigI , 1._wp )         ! Attention! it is -P~ w.r.t paper!

            !! --- "multiplicator" ---
            pmult(ji,jj)  = MIN( zlamb/(zlamb + pdt*(1._wp - zPtld)) , 1._wp - reps24 ) ! Multiplicator term [Eq.32]

      END_2D

   END SUBROUTINE PHASE_I


   SUBROUTINE MOHR_COULOMB_DMG( pdt, pcohe, pNlim, pE, pdx,  ps11, ps22, ps12, pdmg )
      !!
      REAL(wp),                 INTENT(in)    :: pdt              ! (small) time-step [s]
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pcohe            ! cohesion
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pNlim            ! N
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pE               ! Elasticity of damaged ice
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pdx              ! Local grid resolution [m]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11, ps22, ps12 ! Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pdmg             ! damage
      !!
      REAL(wp) :: zsigI, ztmp, ztA, ztB, zsigII, zMC
      REAL(wp) :: zdcrit, zmul, zc1, zc2, zsqrtE, zTd
      REAL(wp) :: rEdmgd
      INTEGER  :: ji, jj
      !!
      !! Maybe a superstition, but I want to avoid a IF block inside the upcoming DO loop:
      rEdmgd = 0.
      IF( ln_damaged_E ) rEdmgd = 1.

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

            zsqrtE = rEdmgd * SQRT(MAX(pE(ji,jj),reps6))   +   (1. - rEdmgd) * rsqrt_E0  ! `sqrt(E)` damaged or undamaged...
            zTd    = MAX( pdx(ji,jj) * rsqrt_nu_rhoi / zsqrtE , reps6 )              ! characteristic time for damage [s] |  (we shall divide by it)...

            zsigI = 0.5_wp * (ps11(ji,jj)+ps22(ji,jj))
            !
            ztmp  = 0.5_wp * (ps11(ji,jj)-ps22(ji,jj))
            ztA   = ztmp*ztmp
            ztB   = ps12(ji,jj)*ps12(ji,jj)
            zsigII  = SQRT( ztA + ztB )

            zMC =  zsigII + rz_muMC*zsigI                         ! Mohr-Coulomb  [Eq.29.2]
            zMC  = SIGN(1._wp, zMC) * MAX( ABS(zMC), reps12 ) ! get rid of values too close to 0 for upcoming division

            zc1 = 0.5_wp + SIGN( 0.5_wp, -pNlim(ji,jj) - zsigI ) !  if (zsigI < -pNlim) => zc1=1. else => zc1=0.

            zdcrit =  (1._wp - zc1) *  pcohe(ji,jj) / zMC                      &  ! `c/MC`   default value for `d_crit` [Eq.29.2]
               &     +      zc1     * -pNlim(ji,jj) / MIN( zsigI , -1._wp )       ! `-N/sigI`[Eq.29.1]

            zc1 = 0.5_wp + SIGN( 0.5_wp, zdcrit - reps24          ) ! if (zdcrit > 0) => zc1=1. else zc1=0.
            zc2 = 0.5_wp + SIGN( 0.5_wp,  1._wp - reps24 - zdcrit ) ! if (zdcrit < 1) => zc2=1. else zc2=0.
            zmul = (1._wp - zdcrit) * pdt / zTd * zc1*zc2           ! Multiplicator for updating stress and damage: `(1 - d_crit)*(pdt/t_d)`

            pdmg(ji,jj) =  MIN( pdmg(ji,jj) + (1._wp - pdmg(ji,jj)) * zmul , rn_dmg_max )
            ps11(ji,jj) =       ps11(ji,jj) -          ps11(ji,jj)  * zmul
            ps22(ji,jj) =       ps22(ji,jj) -          ps22(ji,jj)  * zmul
            ps12(ji,jj) =       ps12(ji,jj) -          ps12(ji,jj)  * zmul

      END_2D

   END SUBROUTINE MOHR_COULOMB_DMG





   SUBROUTINE CROSS_NUDGING( kts, pdt, pAt, pAf,  ps11t, ps22t, ps12t, ps11f, ps22f, ps12f )
      !!
      INTEGER,                  INTENT(in)    :: kts                  ! current small sub-time-step
      REAL(wp),                 INTENT(in)    :: pdt                  ! small time-step [s]
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pAt, pAf             ! ice concentration @T and @F
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11t, ps22t, ps12t  ! T-centric Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11f, ps22f, ps12f  ! F-centric Sigmas [Pa]
      !!
      REAL(wp) :: zcnc ! cross-nudging coefficient to use
      !
      zcnc = rn_crndg * pdt / rdt_ice
      !
      IF( MOD(kts,2) == 0 ) THEN
         !! Only T-centric stress comp. are corrected w.r.t F-centric stress comp.
         IF(l_CN_is_2d) THEN
            ps11t(:,:) = ps11t(:,:) - xCNt(:,:) * ( ps11t(:,:) - rmpF2T( ps11f, lconserv=.TRUE. ) )
            ps22t(:,:) = ps22t(:,:) - xCNt(:,:) * ( ps22t(:,:) - rmpF2T( ps22f, lconserv=.TRUE. ) )
            ps12f(:,:) = ps12f(:,:) - xCNf(:,:) * ( ps12f(:,:) - rmpT2F( ps12t, lconserv=.TRUE. ) )
         ELSE
            ps11t(:,:) = ps11t(:,:) - zcnc * ( ps11t(:,:) - rmpF2T( ps11f, lconserv=.TRUE. ) ) * xmskt(:,:)
            ps22t(:,:) = ps22t(:,:) - zcnc * ( ps22t(:,:) - rmpF2T( ps22f, lconserv=.TRUE. ) ) * xmskt(:,:)
            ps12f(:,:) = ps12f(:,:) - zcnc * ( ps12f(:,:) - rmpT2F( ps12t, lconserv=.TRUE. ) ) * xmskf(:,:)
         ENDIF
         !#tested! => ok! CALL lbc_lnk( 'CROSS_NUDGING_T@icedyn_rhg_bbm',  ps11t,'T',1._wp,  ps22t,'T',1._wp,  ps12f,'F',1._wp )
         CALL clean_small_a_sgm( 'T', pAt, pAf,  ps11t, ps22t, ps12f )
         !!
      ELSE
         !! Only F-centric stress comp. are corrected w.r.t T-centric stress comp.
         IF(l_CN_is_2d) THEN
            ps11f(:,:) = ps11f(:,:) - xCNf(:,:) * ( ps11f(:,:) - rmpT2F( ps11t, lconserv=.TRUE. ) )
            ps22f(:,:) = ps22f(:,:) - xCNf(:,:) * ( ps22f(:,:) - rmpT2F( ps22t, lconserv=.TRUE. ) )
            ps12t(:,:) = ps12t(:,:) - xCNt(:,:) * ( ps12t(:,:) - rmpF2T( ps12f, lconserv=.TRUE. ) )
         ELSE
            ps11f(:,:) = ps11f(:,:) - zcnc * ( ps11f(:,:) - rmpT2F( ps11t, lconserv=.TRUE. ) ) * xmskf(:,:)
            ps22f(:,:) = ps22f(:,:) - zcnc * ( ps22f(:,:) - rmpT2F( ps22t, lconserv=.TRUE. ) ) * xmskf(:,:)
            ps12t(:,:) = ps12t(:,:) - zcnc * ( ps12t(:,:) - rmpF2T( ps12f, lconserv=.TRUE. ) ) * xmskt(:,:)
         ENDIF
         !#tested! => ok! CALL lbc_lnk( 'CROSS_NUDGING_F@icedyn_rhg_bbm',  ps11f,'F',1._wp,  ps22f,'F',1._wp,  ps12t,'T',1._wp )
         CALL clean_small_a_sgm( 'F', pAt, pAf,  ps11f, ps22f, ps12t )
         !!
      ENDIF
   END SUBROUTINE CROSS_NUDGING


#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!==============================================================================
END MODULE icedyn_rhg_bbm
