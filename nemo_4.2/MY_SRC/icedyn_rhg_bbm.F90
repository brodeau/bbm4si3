MODULE icedyn_rhg_bbm
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_bbm  ***
   !!   Sea-Ice dynamics : rheology Britle Maxwell X
   !!======================================================================
   !! History :
   !!            4.2  !  2024     (L. Brodeau) original code
   !!                   `BBM` [starting from `icedyn_rhg_evp`]
   !!
   !!                   Please cite: Brodeau et al. 2024,
   !!                           "Implementation of a brittle sea-ice rheology in an Eulerian,
   !!                            finite-difference, C-grid modeling framework: Impact on the
   !!                            simulated deformation of sea-ice in the Arctic"
   !!                            Geoscientific Model Development (GMD)
   !!                     =>    https://doi.org/10.5194/gmd-17-6051-2024
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
   USE sbc_oce , ONLY : nn_fsbc, ssh_m
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
      &                   rz_muMC  = 0.7_wp          !: slope of Mohr-Coulomb enveloppe

   REAL(wp),  SAVE :: ridlzd
   REAL(wp),  SAVE :: rk0  ! factor to stiffness matrix => 1._wp / ( 1._wp - rz_nup*rz_nup)
   REAL(wp),  SAVE :: rk11, rk22, rk12, rk33 ! elements of stiffness matrix
   REAL(wp),  SAVE :: rlambda0, rsqrt_nu_rhoi, rsqrt_E0, rCe0 ! Constant part of Eq.28

   !! Arrays to be allocated into `ice_dyn_rhg_bbm_init()`:
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Uu_sub, Uv_sub, Vv_sub, Vu_sub !: ice velocities that evolve at sub-time-step
   !!    those that remain constant:
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Xdxt, Xdxf, Xcohst, Xcohsf, XNlimt, XNlimf
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: Xe1t2, Xe2t2, Xe1f2, Xe2f2   !optimization, will avoid these array multiplications countless times...
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xtcoast, xfcoast ! to prevent doing cross-nudging at the coast!
   REAL(wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xCNt, xCNf   ! cross nudging coefficients (space-dependant)

   LOGICAL, PARAMETER :: l_use_v_for_h = .TRUE.  !#LB: we should not! But `h` is sloppy by essence in SI3...

   REAL(wp), SAVE :: zrhoco      !: sea-water density x ocean-ice drag coeff. [kg/m3]
   REAL(wp), SAVE :: zdtbbm, z1_dtbbm !: small time step (time splitting) [s] and its inverse [1/s]

   LOGICAL,  SAVE :: l_CN        !: whether cross nudging is used ?
   LOGICAL,  SAVE :: l_CN_is_2d  !: whether cross nudging coefficient is a 2D array, not a scalar
   REAL(wp), SAVE :: rCNC_eff    !: effective cross-nudging coefficient [-]

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
      REAL(wp) :: zr, zdx_m, zce, zdts
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zt1, zt2, zt3, zt4
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE ::   i2d
      INTEGER :: jm
      !!-------------------------------------------------------------------
      IF( lwp ) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) '    *** Initialization of BBM rheology (ice_dyn_rhg_bbm_init) ***'
      ENDIF

      ridlzd = 1._wp
      IF(ln_idealized) ridlzd = 0._wp
      IF( lwp ) WRITE(numout,*) '  * Disregarding Coriolis and SSH terms in momentum eq.:',ln_idealized,'=> ridlzd =',INT(ridlzd,1)

      zrhoco   = rho0 * rn_cio
      zdtbbm   = rdt_ice / REAL( nn_nbbm, wp )
      z1_dtbbm = 1._wp / zdtbbm

      l_CN       = ( rn_crndg > 0._wp )
      l_CN_is_2d = ( ln_boost_CN_coast .OR. ln_boost_CN_high_dmg ) ! cross nudging coefficient will be a 2D array, not a scalar

      rCNC_eff = rn_crndg / REAL( nn_nbbm, wp )


      !! Masks at T- & F-points for the ice ( rn_shlat MUST = 2 ! )
      !! **********************************
      xmskt(:,:) =     tmask(:,:,1)
      !!
      !!  => for F-points it's tricky, because technically, the F-point is an interface between ocean and land, and
      !!     as such, it can be both considered as ocean and land depending on the variable of interest
      !!     => So we define two of them:
      xmskf(:,:) = MIN(fmask(:,:,1), 1._wp) ! ==> the variable may have a value other than 0 at the F-point (i.e. compressive stresse!)
      !!                                    !     => the `2`s of `fmask` are changed in `1`s !
      ALLOCATE( i2d(jpi,jpj) )
      i2d(:,:) = 0
      i2d(Nis0:Nie0,Njs0:Nje0) = INT( fmask(Nis0:Nie0,Njs0:Nje0,1) &
         &  + fmask(Nis0+1:Nie0+1,Njs0:Nje0,1) + fmask(Nis0:Nie0,Njs0+1:Nje0+1,1) &
         &  + fmask(Nis0-1:Nie0-1,Njs0:Nje0,1) + fmask(Nis0:Nie0,Njs0-1:Nje0-1,1) , 1 )
      WHERE( i2d == 4 ) xmskf = 1._wp
      DEALLOCATE ( i2d )
      !!
      xms0f(:,:) = 1._wp-ABS(fmask(:,:,1)-1._wp) ! ==> the variable has to be 0 at the F-point (i.e. a velocity module!)
      !!                                         !     => the `2`s of `fmask` are changed in `0`s !
      CALL lbc_lnk( 'ice_dyn_rhg_bbm_init',  xmskt,'T',1._wp,  xmskf,'F',1._wp,  xms0f,'F',1._wp )
      
      !! Stiffness matrix
      !! ****************
      rk0  = 1._wp / ( 1._wp - rz_nup*rz_nup)
      rk11 = rk0
      rk12 = rk0 * rz_nup
      rk22 = rk0
      rk33 = rk0 * (1._wp - rz_nup)

      rlambda0 = rn_eta0 / rn_E0     !: Viscosity / Elasticity of undamaged ice (aka relaxation time) [s]
      IF( lwp ) WRITE(numout,*) '  * Viscous relaxation time scale => rlambda0 =', REAL(rlambda0,4), ' [s]'

      rsqrt_nu_rhoi = SQRT( 2._wp*(1._wp + rz_nup)*rhoi )
      rsqrt_E0      = SQRT( rn_E0 )

      ALLOCATE( zt1(jpi,jpj), zt2(jpi,jpj), utauVice(jpi,jpj), vtauUice(jpi,jpj) ,  STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate `utauVice,vtauUice`' )

      ALLOCATE(    Xdxt(jpi,jpj),  Xdxf(jpi,jpj), &
         &      Xcohst(jpi,jpj), Xcohsf(jpi,jpj), XNlimt(jpi,jpj), XNlimf(jpi,jpj),               &
         &      Xe1t2(jpi,jpj), Xe2t2(jpi,jpj), Xe1f2(jpi,jpj), Xe2f2(jpi,jpj),     STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate arrays' )

      Xe1t2(:,:) = e1t(:,:) * e1t(:,:) * xmskt(:,:)
      Xe2t2(:,:) = e2t(:,:) * e2t(:,:) * xmskt(:,:)
      Xe1f2(:,:) = e1f(:,:) * e1f(:,:) * xmskf(:,:)
      Xe2f2(:,:) = e2f(:,:) * e2f(:,:) * xmskf(:,:)

      ALLOCATE( Uu_sub(jpi,jpj), Uv_sub(jpi,jpj), Vv_sub(jpi,jpj), Vu_sub(jpi,jpj),       STAT=ierror )
      IF( ierror /= 0 )  CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: failed to allocate ice velocity arrays' )

      ALLOCATE( uVice(jpi,jpj) , vUice(jpi,jpj) , af_i(jpi,jpj) , STAT=ierror )
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

      !! Grid resolution for T- and F-centric cells:
      Xdxt(:,:) = 0._wp ; Xdxf(:,:) = 0._wp
      Xdxt(:,:) = SQRT( e1t(:,:)*e2t(:,:) )  ! Local `dx` of grid cell [m]
      Xdxf(:,:) = SQRT( e1f(:,:)*e2f(:,:) )  ! Local `dx` of grid cell [m]

      !! A typical `dx` for the Arctic:
      zt1(:,:) = 0._wp
      WHERE(gphit(:,:) > 55._wp) zt1(:,:) = 1._wp
      zt1(:,:) = zt1(:,:)*xmskt(:,:)
      zdx_m = SUM(Xdxt(:,:)*zt1(:,:)) / MAX( SUM(zt1(:,:)) , epsi06 ) ; ! => mean `dx` on this proc domain!
      CALL mpp_sum( 'ice_dyn_rhg_bbm_init', zdx_m)
      zdx_m = zdx_m / REAL(jpnij)  ; ! => mean `dx` of full computational domain

      zce = rsqrt_E0 / rsqrt_nu_rhoi ! propagation speed of shearing elastic waves based on the mean `dx`
      zdts = 0.5_wp*(zdx_m/zce)        ! largest possible small time-step to consider...

      IF( lwp ) THEN
         WRITE(numout,*) '  * Big time step (advection & thermo)  => rdt_ice  =', rdt_ice, ' [s]'
         WRITE(numout,*) '  * Average `dx` of computational domain north of 55N => ', REAL(zdx_m/1000._wp,4), ' [km]'
         WRITE(numout,*) '     ==> propagation speed of shearing elastic waves =>',  zce, '[m/s]'
         WRITE(numout,*) '     ==> the small time-step should therefore be about or below:', zdts, ' [s]'
         WRITE(numout,*) '     ==> that would the case with a `nbbm` >',INT(rdt_ice/zdts,1)
         WRITE(numout,*) '     ==> you have chosen `nbbm` =', nn_nbbm
         WRITE(numout,*) '  * Small time step (rheology) => zdtbbm =', zdtbbm,  ' [s]'
         IF(l_CN) THEN
            WRITE(numout,*) '  * About cross-nudging:'
            WRITE(numout,*) '      - CN parameter (gamma) => rn_crndg =', rn_crndg,' [-]'
            !WRITE(numout,*) '      - nudging coefficient  =>    k_cn  =', rn_crndg/rdt_ice,' [s**-1]'
            !WRITE(numout,*) '      - effective CN coeff.  => rCNC_eff (= k_cn*dt) =', rCNC_eff,' [-]'
            !WRITE(numout,*) '      - theoretical half-life of ajustment reached after', &
            !   &                  INT( 2._wp*LOG(2._wp)*REAL(nn_nbbm,wp)/rn_crndg ,2 ), ' sub-time-steps!'
         ENDIF
         WRITE(numout,*) '  * (scaled) Compression threshod => N_lim =',REAL(rn_Nref*SQRT( rn_l_ref/zdx_m ),4),' [Pa]'
         WRITE(numout,*) ''
      ENDIF

      IF(zdtbbm>zdts) CALL ctl_stop( 'STOP', 'ice_dyn_rhg_bbm_init: `nn_nbbm` is probably to small' )

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

      IF( lwp ) THEN
         WRITE(numout,*) '**********************************************************************'
         WRITE(numout,*) ''
      ENDIF
   END SUBROUTINE ice_dyn_rhg_bbm_init


   SUBROUTINE ice_dyn_rhg_bbm( kt, Kmm, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i )
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
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pshear_i  , pdivu_i
      !!
      INTEGER ::   ji, jj       ! dummy loop indices
      INTEGER ::   jter         ! local integers
      !
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
      REAL(wp), DIMENSION(jpi,jpj) ::   zfUu  , zfVv                    ! divergence of stress tensor (vector) in T-centric grid
      REAL(wp), DIMENSION(jpi,jpj) ::   zfUv  , zfVu                    ! divergence of stress tensor (vector) in F-centric grid
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
      REAL(wp), DIMENSION(jpi,jpj) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, zAt, zAf, zht, zhf, z1_h_t, z1_h_f
      REAL(wp) :: zh, zt1, zt2, ztame, zmacc
      REAL(wp) :: zCorio, zUi, zVi, zMdt, zUo, zVo
      !!-------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_dyn_rhg_bbm: BBM sea-ice rheology'

      !------------------------------------------------------------------------------!
      ! 1) define some variables and initialize arrays
      !------------------------------------------------------------------------------!

      !! Ice concentration and thicknes @T we are going to work with:
      zAt(:,:) = at_i(:,:)
      IF( l_use_v_for_h ) THEN
         zht(:,:) = MAX( vt_i(:,:) , 0._wp)
      ELSE
         !! => actual ice thickness (`hm_i` cannot be used here because = 0!)
         ztmp1(:,:) = 0._wp
         WHERE( at_i(:,:) > epsi20 ) ztmp1(:,:) = 1._wp / at_i(:,:) ! `ztmp1 == 1/A`
         zht(:,:) = MAX( vt_i(:,:) * ztmp1(:,:) , 0._wp)
         WHERE( at_i(:,:) <= 1.E-3_wp ) zht(:,:) = 0._wp
      ENDIF

      zAf(:,:) = MIN( MAX( rmpT2F( zAt,  lconserv=.TRUE. ) , 0._wp ), rn_amax_n ) * xmskf(:,:) ! Ice conc. at F-points #fixme: add south!
      zhf(:,:) =      MAX( rmpT2F( zht,  lconserv=.TRUE. ) , 0._wp )              * xmskf(:,:) ! Ice thickness at F-points
      IF( .NOT. l_use_v_for_h ) THEN
         WHERE( zAf <= 1.E-3_wp ) zhf = 0._wp
      END IF
      CALL lbc_lnk( 'icedyn_rhg_bbm',  zAf,'F',1._wp,  zhf,'F',1._wp )
      af_i(:,:) = zAf(:,:) ! => used for advection of `dmgf` !

      ! For diagnostics:
      xmsk_ice_t(:,:) = 0._wp ; xmsk_ice_f(:,:) = 0._wp
      WHERE( zAt(:,:) >= rclean_below_A ) xmsk_ice_t(:,:) = 1._wp
      WHERE( zAf(:,:) >= rclean_below_A ) xmsk_ice_f(:,:) = 1._wp
      xmsk_ice_t(:,:) = xmsk_ice_t(:,:)*xmskt(:,:)
      xmsk_ice_f(:,:) = xmsk_ice_f(:,:)*xmskf(:,:)


      !------------------------------------------------------------------------------!
      ! 2) Wind / ocean stress, mass terms, coriolis terms
      !------------------------------------------------------------------------------!
      ! sea surface height
      !    embedded sea ice: compute representative ice top surface
      !    non-embedded sea ice: use ocean surface for slope calculation
      zsshdyn(:,:) = ice_var_sshdyn( ssh_m, snwice_mass, snwice_mass_b )

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
         IF( iom_use('cncoeff_t') ) CALL iom_put( 'cncoeff_t' , ztmp3(:,:)*xmsk_ice_t )
         IF( iom_use('cncoeff_f') ) CALL iom_put( 'cncoeff_f' , ztmp4(:,:)*xmsk_ice_f )
         !
         xCNt(:,:) = ztmp3(:,:) / REAL( nn_nbbm, wp )
         xCNf(:,:) = ztmp4(:,:) / REAL( nn_nbbm, wp )
         !
      ENDIF !IF( l_CN_is_2d )

      ! --- Healing of damage with time [Eq.30, Olason et al., 2022]
      !     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ztmp1(:,:) = 1._wp / MAX(zAt(:,:),epsi06)
      ztmp3(:,:) = rcnd_i*vt_s(:,:)*ztmp1(:,:) / MAX( rn_cnd_s*vt_i(:,:)*ztmp1(:,:), epsi06 ) * xmskt(:,:)   ! => `C` of the `dtemp/(1 + C)` in neXtSIM
      IF(iom_use('ice_heal_c'))  CALL iom_put( 'ice_heal_c' , ztmp3 )
      IF(ln_icethd) THEN
         !! Thermo is used, normal stuff:
         ztmp1(:,:) = (t_bo(:,:) - tm_su(:,:)) / (1._wp + ztmp3(:,:) ) * xmskt(:,:) ! temperature difference between bottom and surface
      ELSE
         !! Thermo is off, yet we want som refreezing!
         ztmp1(:,:) = (-1.8_wp + 25._wp)/ (1._wp + ztmp3(:,:) ) * xmskt(:,:) ! temperature difference between bottom and surface
      END IF
      ztmp2(:,:) = MIN( ztmp1(:,:) / rn_kth , 1._wp/rdt_ice )  ! 1/T_relax
      IF(iom_use('ice_heal_x'))  CALL iom_put( 'ice_heal_x' , ztmp2 )
      WHERE( ztmp1(:,:) > 0._wp ) dmgt(:,:) = MAX( dmgt(:,:) - rdt_ice*ztmp2(:,:)      , 0._wp )
      ztmp5(:,:) = rmpT2F( ztmp1 )
      WHERE( ztmp5(:,:) > 0._wp ) dmgf(:,:) = MAX( dmgf(:,:) - rdt_ice*rmpT2F( ztmp2 ) , 0._wp )
      CALL cap_damage( 'T', 'update_stress_dmg', dmgt )
      CALL cap_damage( 'F', 'update_stress_dmg', dmgf )


      !! Because there has been advection of `damage` and stress tensors since last time we called `clean_small_a_all`:
      CALL clean_small_a_all( zAt, zAf,  dmgt, dmgf,  sgm11t, sgm22t, sgm12t,  sgm11f, sgm22f, sgm12f )

      ! Going to average (set to 0 before accumulating during the `nn_nbbm` sub time steps):
      u_ice(:,:) = 0._wp
      uVice(:,:) = 0._wp
      v_ice(:,:) = 0._wp
      vUice(:,:) = 0._wp
      zmacc      = 1._wp/REAL(nn_nbbm)
      !
      IF(iom_use('fUu')) ztmp4(:,:) = 0._wp
      IF(iom_use('fVv')) ztmp5(:,:) = 0._wp

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
         CALL div_stress_tensor( 'T', Xe1t2, Xe2t2,  Xe1f2, Xe2f2,  r1_e2u, r1_e1u, r1_e1v, r1_e2v,  r1_e1e2u, r1_e1e2v,  &
            &                         ztmp1, ztmp2, ztmp3,  zfUu, zfVv )

         ztmp1(:,:) = sgm11f(:,:) * zhf(:,:)
         ztmp2(:,:) = sgm22f(:,:) * zhf(:,:)
         ztmp3(:,:) = sgm12t(:,:) * zht(:,:)
         !!
         CALL div_stress_tensor( 'F', Xe1f2, Xe2f2,  Xe1t2, Xe2t2,  r1_e2v, r1_e1v, r1_e1u, r1_e2u,  r1_e1e2v, r1_e1e2u,  &
            &                         ztmp1, ztmp2, ztmp3,  zfUv, zfVu )

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

         CALL lbc_lnk( 'icedyn_rhg_bbm', Uu_sub,'U',-1._wp, Vv_sub,'V',-1._wp, &
            &                            Uv_sub,'V',-1._wp, Vu_sub,'U',-1._wp )

         ! Average the velocity (to be used for advection at the big time step):
         u_ice(:,:) = u_ice(:,:) + zmacc*Uu_sub(:,:)
         v_ice(:,:) = v_ice(:,:) + zmacc*Vv_sub(:,:)
         uVice(:,:) = uVice(:,:) + zmacc*Uv_sub(:,:)
         vUice(:,:) = vUice(:,:) + zmacc*Vu_sub(:,:)

         ! Average the 2 components of the divergence of the stress tensor for T-centric cell (for diagnostics):
         IF(iom_use('fUu')) ztmp4(:,:) = ztmp4(:,:) + zmacc*zfUu(:,:)
         IF(iom_use('fVv')) ztmp5(:,:) = ztmp5(:,:) + zmacc*zfVv(:,:)
         !
         !                                                ! ==================== !
      END DO !DO jter = 1 , nn_nbbm                       !  end loop over jter  !
      !                                                   ! ==================== !

      IF(iom_use('fUu'))  CALL iom_put( 'fUu' , ztmp4 )
      IF(iom_use('fVv'))  CALL iom_put( 'fVv' , ztmp5 )

      !! Closer look at the components of rate-of-strain tensor:
      IF( iom_use('e11t').OR.iom_use('e22t').OR.iom_use('e12t') ) THEN
         CALL strain_rate( 'T', u_ice, v_ice, uVice, vUice, &
            &              r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, Xe1t2, Xe2t2, tmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3, lblnk=.FALSE. )
         IF( iom_use('e11t') ) CALL iom_put( 'e11t' , ztmp1*xmsk_ice_t )
         IF( iom_use('e22t') ) CALL iom_put( 'e22t' , ztmp2*xmsk_ice_t )
         IF( iom_use('e12t') ) CALL iom_put( 'e12t' , ztmp3*xmsk_ice_t )
      ENDIF
      IF( iom_use('e11f').OR.iom_use('e22f').OR.iom_use('e12f') ) THEN
         CALL strain_rate( 'F', uVice, vUice, u_ice, v_ice, &
            &              r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3,  lblnk=.FALSE. )
         IF( iom_use('e11f') ) CALL iom_put( 'e11f' , ztmp1*xmsk_ice_f )
         IF( iom_use('e22f') ) CALL iom_put( 'e22f' , ztmp2*xmsk_ice_f )
         IF( iom_use('e12f') ) CALL iom_put( 'e12f' , ztmp3*xmsk_ice_f )
      END IF

      !! Saving stress tensor components in the proper units! i.e. Pa :
      IF( iom_use('ice_sig11')  ) CALL iom_put( 'ice_sig11' ,  sgm11t*xmsk_ice_t )
      IF( iom_use('ice_sig22')  ) CALL iom_put( 'ice_sig22' ,  sgm22t*xmsk_ice_t )
      IF( iom_use('ice_sig12')  ) CALL iom_put( 'ice_sig12' ,  sgm12f*xmsk_ice_f )
      IF( iom_use('ice_sig11f') ) CALL iom_put( 'ice_sig11f',  sgm11f*xmsk_ice_f )
      IF( iom_use('ice_sig22f') ) CALL iom_put( 'ice_sig22f',  sgm22f*xmsk_ice_f )
      IF( iom_use('ice_sig12t') ) CALL iom_put( 'ice_sig12t',  sgm12t*xmsk_ice_t )
      !
      IF( iom_use('normstr') ) CALL iom_put( 'normstr' , 0.5_wp*(sgm11t+sgm22t)*xmsk_ice_t ) ! First invariant of stress tensor @T
      IF( iom_use('normstrf')) CALL iom_put( 'normstrf', 0.5_wp*(sgm11f+sgm22f)*xmsk_ice_f ) ! First invariant of stress tensor @F
      !
      IF( iom_use('sheastr' )) CALL iom_put( 'sheastr',  sigmaII( sgm11t, sgm22t, sgm12t ) )
      IF( iom_use('sheastrf')) CALL iom_put( 'sheastrf', sigmaII( sgm11f, sgm22f, sgm12f ) )
      !

      !------------------------------------------------------------------------------!
      ! 4) Recompute shear and div (inputs for mechanical redistribution)
      !------------------------------------------------------------------------------!

      CALL strain_rate( 'T', u_ice, v_ice, uVice, vUice, &
         &              r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, Xe1t2, Xe2t2, tmask(:,:,1), &
         &              ztmp1, ztmp2, ztmp3, lblnk=.TRUE., pdiv=pdivu_i, pmaxshr=pshear_i )
      ! --- divergence of velocity field @T:
      IF( iom_use('icediv') )  CALL iom_put( 'icediv' , pdivu_i*xmsk_ice_t )
      ! --- shear of velocity field @T:
      IF( iom_use('iceshr') )  CALL iom_put( 'iceshr' , ztmp3*xmsk_ice_t )
      ! --- MAXIMUM shear of velocity field @T:
      IF( iom_use('iceshe') )  CALL iom_put( 'iceshe' , pshear_i*xmsk_ice_t )
      ! --- total deformation of velocity field @T:
      IF( iom_use('icedef') ) CALL iom_put( 'icedef', SQRT( pshear_i*pshear_i + pdivu_i*pdivu_i )*xmsk_ice_t )

      IF( iom_use('icedivf') .OR. iom_use('iceshrf') .OR. iom_use('iceshef') .OR. iom_use('icedeff') ) THEN
         CALL strain_rate( 'F', uVice, vUice, u_ice, v_ice, &
            &              r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
            &              ztmp1, ztmp2, ztmp3,  lblnk=.TRUE., pdiv=zfUu, pmaxshr=zfVv )
         ! --- divergence of velocity field @F:
         IF( iom_use('icedivf') ) CALL iom_put( 'icedivf' , zfUu*xmsk_ice_f )
         ! --- shear of velocity field @F:
         IF( iom_use('iceshrf') ) CALL iom_put( 'iceshrf' , ztmp3*xmsk_ice_f )
         ! --- MAXIMUM shear of velocity field @F:
         IF( iom_use('iceshef') ) CALL iom_put( 'iceshef' , zfVv*xmsk_ice_f )
         ! --- total deformation of velocity field @F:
         IF( iom_use('icedeff') ) CALL iom_put( 'icedeff', SQRT( zfVv*zfVv + zfUu*zfUu )*xmsk_ice_f )

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
         CALL iom_put( 'icevorf' , ztmp3*xmsk_ice_f )
      ENDIF
      ! --- vorticity of velocity field @T:
      IF( iom_use('icevort') ) THEN
         DO_2D( 0, 0, 0, 0 )
               ztmp3(ji,jj) = (   ( vUice(ji,jj)*r1_e2u(ji,jj) - vUice(ji-1,jj)*r1_e2u(ji-1,jj) ) * Xe2t2(ji,jj) &
                  &             - ( uVice(ji,jj)*r1_e1v(ji,jj) - uVice(ji,jj-1)*r1_e1v(ji,jj-1) ) * Xe1t2(ji,jj) &
                  &           ) * r1_e1e2t(ji,jj) * tmask(ji,jj,1)
         END_2D
         CALL iom_put( 'icevort' , ztmp3*xmsk_ice_t )
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
         CALL iom_put( 'taum_ai' , SQRT(ztmp1*ztmp1 + ztmp2*ztmp2) * zAt * xmsk_ice_t(:,:) )
      END IF

      IF( iom_use('utauVai') .OR. iom_use('vtauUai') ) THEN
         CALL iom_put( 'utauVai' , ztauxVai )
         CALL iom_put( 'vtauUai' , ztauyUai )
      ENDIF
      IF( iom_use('taumFai') ) CALL iom_put( 'taumFai' , SQRT(ztauxVai*ztauxVai + ztauyUai*ztauyUai)*xmsk_ice_f )  !#fixme: ugly!!!

      ! --- ice-ocean stress (multiplied by `A`, just as in momentum equation !):
      IF( iom_use('taum_oi') .OR. iom_use('utau_oi') .OR. iom_use('vtau_oi') ) THEN
         ztmp1(:,:) = umask(:,:,1)
         WHERE( zAu(:,:) < rclean_below_A ) ztmp1(:,:) = 0._wp
         ztmp3(:,:) = u_ice(:,:) - u_oce(:,:) ! dU @U
         ztmp4(:,:) = vUice(:,:) - vUoce(:,:) ! dV @U
         ztmp2(:,:) = SQRT( ztmp3(:,:)*ztmp3(:,:) + ztmp4(:,:)*ztmp4(:,:) ) ! module of relative current at U-point
         zfUu (:,:) = zAu(:,:) * zrhoco * ztmp2(:,:) * ( u_oce(:,:) - u_ice(:,:) ) * ztmp1
         !
         ztmp1(:,:) = vmask(:,:,1)
         WHERE( zAv(:,:) < rclean_below_A ) ztmp1(:,:) = 0._wp
         ztmp3(:,:) = v_ice(:,:) - v_oce(:,:) ! dV @V
         ztmp4(:,:) = uVice(:,:) - uVoce(:,:) ! dU @V
         ztmp2(:,:) = SQRT( ztmp3(:,:)*ztmp3(:,:) + ztmp4(:,:)*ztmp4(:,:) ) ! module of relative current at V-point
         zfVv (:,:) = zAv(:,:) * zrhoco * ztmp2(:,:) * ( v_oce(:,:) - v_ice(:,:) ) * ztmp1
         !
         IF( iom_use('taum_oi') ) THEN
            !! Module of O/I stress at T points
            ztmp1(2:jpi,:) = 0.5_wp * ( zfUu(2:jpi,:) + zfUu(1:jpi-1,:) )
            ztmp2(:,2:jpj) = 0.5_wp * ( zfVv(:,2:jpj) + zfVv(:,1:jpj-1) )
            ztmp3(:,:) = SQRT(ztmp1*ztmp1 + ztmp2*ztmp2)
         END IF
         !
         CALL lbc_lnk( 'icedyn_rhg_bbm', zfUu,'U',-1._wp, zfVv,'V',-1._wp, ztmp3,'T',1._wp )
         IF( iom_use('utau_oi') ) CALL iom_put( 'utau_oi' , zfUu(:,:) )  ! oce-ice stress /x @ U. MIND: x A !!!
         IF( iom_use('vtau_oi') ) CALL iom_put( 'vtau_oi' , zfVv(:,:) )  ! oce-ice stress /y @ V MIND: x A !!!
         IF( iom_use('taum_oi') ) CALL iom_put( 'taum_oi' , ztmp3 * xmsk_ice_t )
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
      REAL(wp), DIMENSION(jpi,jpj) :: zelat, zelaf, zlambt, zlambf, zmltt, zmltf
      REAL(wp)                     :: zA, zmlt
      REAL(wp)                     :: ze11, ze22, ze12
      INTEGER  :: ji, jj
      LOGICAL, PARAMETER :: llbclnk=.FALSE.
      !!----------------------------------------------------------------------

      !! --- First (predictor) update of stress tensor terms [Eq.32] ---
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !! Compute `elasticity`, and `update multiplicator` @ T-points (=> zelat, zmltt)
      CALL UPDATE_E_L_MULT( pdt, pAt, pht, pdmgt, ps11t, ps22t, ps12t,  zelat, zlambt, zmltt )
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
            zmlt = zmltt(ji,jj) * xmskt(ji,jj)
            zA   = zelat(ji,jj) * pdt
            !
            ps11t(ji,jj) = zmlt * ( zA * ( rk11*ze11 + rk12*ze22) + ps11t(ji,jj) )
            ps22t(ji,jj) = zmlt * ( zA * ( rk12*ze11 + rk22*ze22) + ps22t(ji,jj) )
            ps12t(ji,jj) = zmlt * ( zA *        rk33 * ze12       + ps12t(ji,jj) )
      END_2D

      !! Compute `elasticity`, and `update multiplicator` @ F-points (=> zelaf, zmltf)
      CALL UPDATE_E_L_MULT( pdt, pAf, phf, pdmgf, ps11f, ps22f, ps12f,  zelaf, zlambf, zmltf ) ! compute `elasticity` and `update multiplicator` @ F
      !
      !! Compute the 3 components of the strain-rate tensor @ F-points (=> ztp1, ztp2, ztp3):
      CALL strain_rate( 'F', pUv, pVu, pUu, pVv, r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, Xe1f2, Xe2f2, fmask(:,:,1), &
         &              ztp1, ztp2, ztp3, lblnk=.FALSE. ) !: double-checked that lbc_lnk-ing is NOT NEEDED !!!

      !! Predictor update of the 3 stress tensor components @ F-points:
      DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
            ze11 = ztp1(ji,jj)
            ze22 = ztp2(ji,jj)
            ze12 = ztp3(ji,jj)
            zmlt = zmltf(ji,jj) * xmskf(ji,jj)
            zA = zelaf(ji,jj) * pdt
            !
            ps11f(ji,jj) = zmlt * ( zA * ( rk11*ze11 + rk12*ze22) + ps11f(ji,jj) )
            ps22f(ji,jj) = zmlt * ( zA * ( rk12*ze11 + rk22*ze22) + ps22f(ji,jj) )
            ps12f(ji,jj) = zmlt * ( zA *        rk33 * ze12       + ps12f(ji,jj) )
      END_2D

      IF( kts==nn_nbbm ) THEN
         IF( iom_use('zelat') ) CALL iom_put( 'zelat' , zelat(:,:)*xmsk_ice_t(:,:) )
         IF( iom_use('zelaf') ) CALL iom_put( 'zelaf' , zelaf(:,:)*xmsk_ice_f(:,:) )
         IF( iom_use('zetat') ) CALL iom_put( 'zetat' , zlambt(:,:)*zelat(:,:)*xmsk_ice_t(:,:) ) ; ! BBM viscosity `eta`
         IF( iom_use('zetaf') ) CALL iom_put( 'zetaf' , zlambf(:,:)*zelaf(:,:)*xmsk_ice_f(:,:) )
         IF( iom_use('zlambt')) CALL iom_put( 'zlambt', zlambt(:,:)*xmsk_ice_t(:,:) )
         IF( iom_use('zlambf')) CALL iom_put( 'zlambf', zlambf(:,:)*xmsk_ice_f(:,:) )
         IF( iom_use('zmult') ) CALL iom_put( 'zmult' , zmltt(:,:)*xmsk_ice_t(:,:) )
         IF( iom_use('zmulf') ) CALL iom_put( 'zmulf' , zmltf(:,:)*xmsk_ice_f(:,:) )
      ENDIF

      IF( l_CN ) THEN
         !! Cross-nudging on stress tensor components
         !!  Alternate / even/odd sub iterations => gives slightly better results than doing
         !!  it everytime (with half zrbal) on all stress comp.
         !!     (=> double checked that no `lbc_lnk` is needed in `rmpY2X()` !!!)
         IF( MOD(kts,2) == 0 ) THEN
            !! Correcting of T-centric stress tensor components:
            CALL CROSS_NUDGING( 'T', pht, phf, ps11f, ps22f, ps12t,   ps11t, ps22t, ps12f )
            CALL clean_small_a_sgm( 'T', pAt, pAf,  ps11t, ps22t, ps12f )
         ELSE
            !! Correcting of F-centric stress tensor components:
            CALL CROSS_NUDGING( 'F', phf, pht, ps11t, ps22t, ps12f,   ps11f, ps22f, ps12t )
            CALL clean_small_a_sgm( 'F', pAt, pAf,  ps11f, ps22f, ps12t )
         END IF
      ELSE
         IF( (kt==1).AND.(lwp) ) WRITE(numout,*) ' *** MIND: no cross-nudging applied!'
      ENDIF

      !! --- Mohr-Coulomb test and britle update ---
      CALL MOHR_COULOMB_DMG( pdt, Xcohst, XNlimt, zelat, Xdxt, ps11t, ps22t, ps12t, pdmgt )
      CALL MOHR_COULOMB_DMG( pdt, Xcohsf, XNlimf, zelaf, Xdxf, ps11f, ps22f, ps12f, pdmgf )

      CALL clean_small_a_all( pAt, pAf,  pdmgt, pdmgf,  ps11t, ps22t, ps12t,  ps11f, ps22f, ps12f )

      CALL cap_damage( 'T', 'update_stress_dmg', pdmgt )
      CALL cap_damage( 'F', 'update_stress_dmg', pdmgf )

      !! --- Final LBC linking ---
      CALL lbc_lnk( 'UPDATE_STRESS_DMG@icedyn_rhg_bbm', ps11t,'T',1._wp, ps22t,'T',1._wp, ps12f,'F',1._wp, pdmgt,'T',1._wp, &
         &                                              ps11f,'F',1._wp, ps22f,'F',1._wp, ps12t,'T',1._wp, pdmgf,'F',1._wp  )

   END SUBROUTINE update_stress_dmg



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

            IF( MIN( id1, id2, id3, id4, id5, id6 ) > 0 ) THEN      ! fields exist
               CALL iom_get( numrir, jpdom_auto, 'sgm11t', sgm11t, cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22t', sgm22t, cd_type = 'T' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12f', sgm12f, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm11f', sgm11f, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm22f', sgm22f, cd_type = 'F' )
               CALL iom_get( numrir, jpdom_auto, 'sgm12t', sgm12t, cd_type = 'T' )
            ELSE                                     ! start rheology from rest
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) '   ==>>> did not find components of stress tensors in restart file => set to 0'
               sgm11t(:,:) =  0._wp
               sgm22t(:,:) =  0._wp
               sgm12f(:,:) =  0._wp
               sgm11f(:,:) =  0._wp
               sgm22f(:,:) =  0._wp
               sgm12t(:,:) =  0._wp
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


   SUBROUTINE UPDATE_E_L_MULT( pdt, pA, ph, pdmg, ps11, ps22, ps12, pelast, plamb, pmult )
      !!
      REAL(wp),                 INTENT(in)  :: pdt         ! (small) time-step [s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pA          ! Ice concentration
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ph          ! Ice thickness
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pdmg        ! Ice damage
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11, ps22, ps12  ! sigma_11 & sigma_22 [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pelast      ! Ice elasticity
      REAL(wp), DIMENSION(:,:), INTENT(out) :: plamb       ! Ice viscous time scale [s]
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pmult       ! Multiplicator for stress tensor update
      !!
      REAL(wp) :: zxpC, zsigI, zPmax, zc0, zPtld, z1md, zE, zeta, zlamb
      REAL(wp) :: z1_zsigI, zsigII, zang
      INTEGER  :: ji, jj

      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

            zxpC = EXP( rn_C0*(1._wp - pA(ji,jj)) )  ! `expC` [Eq.8]

            z1md = 1._wp - pdmg(ji,jj)

            !! Elasticity [Pa]:
            !!    *** BBM and MEB:
            !!       * E = E0 * (1 - d) * exp[-C*(1-A)]
            zE = rn_E0 * z1md * zxpC

            !! Viscosity [Pa.s]:
            !!    *** MEB (Dansereau et al., 2016):
            !!       * V = V0 * (1 - d)**a * exp[-C*(1-A)]  (viscosity)
            !!    *** BBM (Olason et al. 2022) [Eq.10/Eq.9]:
            !!       *    V = V0 * (1 - d)**a * exp[b*-C*(1-A)]    (with b=a in Olason et al. 2022)
            zeta = rn_eta0 * z1md**rn_alrlx * zxpC**rn_btrlx

            zlamb = zeta / MAX( zE, epsi20 )  ! time-scale for viscous relaxation
            zlamb = MAX( zlamb, pdt )

            zsigI  = 0.5_wp * ( ps11(ji,jj) + ps22(ji,jj) ) ! sigI: normal stress aka first invariant
            zsigII = sigmaII( ps11(ji,jj), ps22(ji,jj), ps12(ji,jj) )

            z1_zsigI = SIGN( 1._wp , zsigI ) / MAX( ABS(zsigI), epsi20 )   ! 1/SigI without the SigI=0 singularity...

            zang   = ATAN( zsigII * z1_zsigI )
            zPmax = -rn_P0 * ph(ji,jj)**1.5_wp * zxpC * COS( zang )       ! `-P_max` (for sigI<0)

            zc0 = 0.5_wp + SIGN( 0.5_wp, -zsigI-epsi20 )           ! => if sigI<-epsi20 => zc0=1 else: zc0=0
            zPtld = zc0 * MIN( zPmax*z1_zsigI , 1._wp )

            pmult(ji,jj)  = MIN( zlamb/(zlamb + pdt*(1._wp - zPtld)) , 1._wp ) ! Multiplicator term [Eq.32]

            pelast(ji,jj) = zE
            plamb(ji,jj) = zlamb

      END_2D

   END SUBROUTINE UPDATE_E_L_MULT


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
      REAL(wp) :: zsigI, zsigII, zMC
      REAL(wp) :: zdcrit, zmlt, zsqrtE, zTd, zc0, zc1, z1_zsigI, z1_zMC, zNlim, zcA, zcB
      INTEGER  :: ji, jj
      !!
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

            zNlim = pNlim(ji,jj)

            zdcrit = 9999._wp ; ! just for initialization

            zsqrtE = SQRT(MAX(pE(ji,jj),epsi06))                               ! `sqrt(E)` (damaged ice)...
            zTd    = MAX( pdx(ji,jj) * rsqrt_nu_rhoi / zsqrtE , epsi06 )       ! characteristic time for damage [s] |  (we shall divide by it)...

            zsigI  = 0.5_wp * (ps11(ji,jj) + ps22(ji,jj))

            zsigII = sigmaII( ps11(ji,jj), ps22(ji,jj), ps12(ji,jj) )

            z1_zsigI = SIGN( 1._wp , zsigI ) / MAX( ABS(zsigI), epsi20 )   ! 1/SigI without the SigI=0 singularity...

            zMC = zsigII + rz_muMC*zsigI                             ! Mohr-Coulomb  [Eq.29.2]
            z1_zMC = SIGN( 1._wp , zMC ) / MAX( ABS(zMC), epsi20 )   ! 1/MC without the MC=0 singularity...

            zc0 = 0.5_wp + SIGN( 0.5_wp , zsigI + zNlim       )   ! if zsigI<-Nlim => zc0=0 ; zc0=1 otherwize
            zdcrit = zc0 * pcohe(ji,jj) * z1_zMC  +  (zc0-1._wp) * zNlim * z1_zsigI   ! `zc0-1` because we need `-Nlim` !

            ! Do only `!IF( (ABS(zsigI)>0.1_wp).AND.(zsigII>0.1_wp) )`:
            zcA = 0.5_wp + SIGN( 0.5_wp , ABS(zsigI)-0.1_wp )   ! if |sigI|>0.1 => zcA=1 ; zcA=0 otherwize
            zcB = 0.5_wp + SIGN( 0.5_wp ,     zsigII-0.1_wp )   ! if  sigII>0.1 => zcB=1 ; zcB=0 otherwize
            zc0 = 0.5_wp + SIGN( 0.5_wp , zdcrit       )   ! if zdcrit>0 => zc0=1 ; zc0=0 otherwize
            zc1 = 0.5_wp + SIGN( 0.5_wp , 1._wp-zdcrit )   ! if zdcrit<1 => zc1=1 ; zc0=0 otherwize

            !IF( (zdcrit > 0._wp).AND.(zdcrit < 1._wp) ) THEN
            !! Damage grows, and stresses are decreased accordingly
            zmlt = zc0*zc1*zcA*zcB * (1._wp - zdcrit) * pdt / zTd          ! Multiplicator for updating stress and damage: `(1 - d_crit)*(pdt/t_d)`
            pdmg(ji,jj) =  MIN( pdmg(ji,jj) + (1._wp - pdmg(ji,jj)) * zmlt , rn_dmg_max )
            ps11(ji,jj) =       ps11(ji,jj) -          ps11(ji,jj)  * zmlt
            ps22(ji,jj) =       ps22(ji,jj) -          ps22(ji,jj)  * zmlt
            ps12(ji,jj) =       ps12(ji,jj) -          ps12(ji,jj)  * zmlt

      END_2D

   END SUBROUTINE MOHR_COULOMB_DMG


   SUBROUTINE CROSS_NUDGING( cgt, ph, phx, ps11x, ps22x, ps12x,   ps11, ps22, ps12 )
      !!
      CHARACTER(len=1),         INTENT(in)    :: cgt                  ! 'T' or 'F' -centric grid ?
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: ph, phx              ! ice thickness at center-point, and corner-point, of relevant grid
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: ps11x, ps22x, ps12x  ! 3 stresses reference on grid [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11,  ps22,  ps12   ! 3 stresses on grid to correct [Pa]
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zs11x, zs22x, zs12x, z1_h, z1_hx
      !!
      !! Nudging must use `sigma*h` rather than `sigma`:
      zs11x(:,:) = ps11x(:,:)*phx(:,:)
      zs22x(:,:) = ps22x(:,:)*phx(:,:)
      zs12x(:,:) = ps12x(:,:)* ph(:,:)
      !!
      !!
      IF( cgt == 'T' ) THEN
         !! Only T-centric stress comp. are corrected w.r.t F-centric stress comp.
         z1_h(:,:)  = 1._wp / MAX( ph(:,:),  epsi06) * xmskt(:,:)
         z1_hx(:,:) = 1._wp / MAX( phx(:,:), epsi06) * xmskf(:,:)
         !!
         IF(l_CN_is_2d) THEN
            ps11(:,:) = ps11(:,:) - xCNt(:,:) * ( ps11(:,:) - rmpF2T( zs11x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps22(:,:) = ps22(:,:) - xCNt(:,:) * ( ps22(:,:) - rmpF2T( zs22x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps12(:,:) = ps12(:,:) - xCNf(:,:) * ( ps12(:,:) - rmpT2F( zs12x, lconserv=.TRUE. ) * z1_hx(:,:))
         ELSE
            ps11(:,:) = ps11(:,:) - rCNC_eff  * ( ps11(:,:) - rmpF2T( zs11x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps22(:,:) = ps22(:,:) - rCNC_eff  * ( ps22(:,:) - rmpF2T( zs22x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps12(:,:) = ps12(:,:) - rCNC_eff  * ( ps12(:,:) - rmpT2F( zs12x, lconserv=.TRUE. ) * z1_hx(:,:))
         ENDIF
         !!
      ELSEIF( cgt == 'F' ) THEN
         !! Only F-centric stress comp. are corrected w.r.t F-centric stress comp.
         z1_h(:,:)  = 1._wp / MAX( ph(:,:),  epsi06) * xmskf(:,:)
         z1_hx(:,:) = 1._wp / MAX( phx(:,:), epsi06) * xmskt(:,:)
         !!
         IF(l_CN_is_2d) THEN
            ps11(:,:) = ps11(:,:) - xCNf(:,:) * ( ps11(:,:) - rmpT2F( zs11x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps22(:,:) = ps22(:,:) - xCNf(:,:) * ( ps22(:,:) - rmpT2F( zs22x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps12(:,:) = ps12(:,:) - xCNt(:,:) * ( ps12(:,:) - rmpF2T( zs12x, lconserv=.TRUE. ) * z1_hx(:,:))
         ELSE
            ps11(:,:) = ps11(:,:) - rCNC_eff  * ( ps11(:,:) - rmpT2F( zs11x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps22(:,:) = ps22(:,:) - rCNC_eff  * ( ps22(:,:) - rmpT2F( zs22x, lconserv=.TRUE. ) * z1_h(:,:) )
            ps12(:,:) = ps12(:,:) - rCNC_eff  * ( ps12(:,:) - rmpF2T( zs12x, lconserv=.TRUE. ) * z1_hx(:,:))
         ENDIF
      ELSE
         CALL ctl_stop( 'STOP', 'CROSS_NUDGING() => wrong type of grid-point: '//cgt )
      ENDIF
      !!
   END SUBROUTINE CROSS_NUDGING


#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!==============================================================================
END MODULE icedyn_rhg_bbm
