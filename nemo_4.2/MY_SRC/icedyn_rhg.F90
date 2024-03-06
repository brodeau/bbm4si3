MODULE icedyn_rhg
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg  ***
   !!   Sea-Ice dynamics : master routine for rheology
   !!======================================================================
   !! history :  4.0  !  2018     (C. Rousset)      Original code
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!    ice_dyn_rhg      : computes ice velocities
   !!    ice_dyn_rhg_init : initialization and namelist read
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain
   USE ice            ! sea-ice: variables
   USE icedyn_rhg_evp ! sea-ice: EVP rheology
   USE icedyn_rhg_eap ! sea-ice: EAP rheology
   USE icedyn_rhg_vp  ! sea-ice: VP  rheology
   USE icedyn_rhg_bbm ! sea-ice: BBM rheology !#bbm
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rhg        ! called by icestp.F90
   PUBLIC   ice_dyn_rhg_init   ! called by icestp.F90



   INTEGER ::              nice_rhg   ! choice of the type of rheology
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_rhgEVP = 1   ! EVP rheology
   INTEGER, PARAMETER ::   np_rhgEAP = 2   ! EAP rheology
   INTEGER, PARAMETER ::   np_rhgVP  = 3   ! VP rheology
   INTEGER, PARAMETER ::   np_rhgBBM = 4   ! BBM rheology

   !
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_rhg.F90 14072 2020-12-04 07:48:38Z laurent $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_rhg( kt, Kmm )
      !!-------------------------------------------------------------------
      !!               ***  ROUTINE ice_dyn_rhg  ***
      !!
      !! ** Purpose :   compute ice velocity
      !!
      !! ** Action  : comupte - ice velocity (u_ice, v_ice)
      !!                      - 3 components of the stress tensor (stress1_i, stress2_i, stress12_i)
      !!                      - shear, divergence and delta (shear_i, divu_i, delta_i)
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ice time step
      INTEGER, INTENT(in) ::   Kmm    ! ocean time level index
      !!--------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_rhg')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'ice_dyn_rhg: sea-ice rheology'
         WRITE(numout,*)'~~~~~~~~~~~'
      ENDIF
      !
      !--------------!
      !== Rheology ==!
      !--------------!
      SELECT CASE( nice_rhg )
      !                                !------------------------!
      CASE( np_rhgEVP )                ! Elasto-Viscous-Plastic !
         !                             !------------------------!
         CALL ice_dyn_rhg_evp( kt, Kmm, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )
         !
         !                             !------------------------!
      CASE( np_rhgVP  )                ! Viscous-Plastic        !
         !                             !------------------------!
         CALL ice_dyn_rhg_vp ( kt, shear_i, divu_i, delta_i )
         !
         !                             !----------------------------!
      CASE( np_rhgEAP )                ! Elasto-Anisotropic-Plastic !
         !                             !----------------------------!
         CALL ice_dyn_rhg_eap( kt, Kmm, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i, aniso_11, aniso_12, rdg_conv )
         !
         !                             !-------------------------! 
      CASE( np_rhgBBM )                ! Brittle Bingham Maxwell !
         !                             !-------------------------!         
         CALL ice_dyn_rhg_bbm( kt, Kmm, stress1_i, stress2_i, stress12_i, shear_i, divu_i, delta_i )
         !
      END SELECT
      !
      IF( lrst_ice ) THEN
         IF( ln_rhg_EVP )   CALL rhg_evp_rst( 'WRITE', kt ) !* write EVP fields in the restart file
         IF( ln_rhg_EAP )   CALL rhg_eap_rst( 'WRITE', kt ) !* write EAP fields in the restart file
         ! MV note: no restart needed for VP as there is no time equation for stress tensor
         IF( ln_rhg_BBM )   CALL rhg_bbm_rst( 'WRITE', kt )
      ENDIF
      !
      ! controls
      IF( sn_cfctl%l_prtctl ) &
         &                 CALL ice_prt3D   ('icedyn_rhg')                                                             ! prints
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_rhg', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icedyn_rhg',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      IF( ln_timing    )   CALL timing_stop ('icedyn_rhg')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_rhg


   SUBROUTINE ice_dyn_rhg_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_rhg_init  ***
      !!
      !! ** Purpose : Physical constants and parameters linked to the ice
      !!      dynamics
      !!
      !! ** Method  :  Read the namdyn_rhg namelist and check the ice-dynamic
      !!       parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_rhg
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_rhg/  ln_rhg_EVP, ln_aEVP, ln_rhg_EAP, rn_creepl, rn_ecc , nn_nevp, rn_relast, nn_rhg_chkcvg, &  !-- evp
         &                  ln_rhg_BBM, ln_MEB, rn_Nref, rn_E0, rn_eta0, rn_P0, rn_kth, nn_nbbm, nn_d_adv,          &  !-- bbm
         &                  rn_crndg, ln_boost_CN_coast, rn_max_CN_coast, ln_boost_CN_high_dmg, rn_max_CN_dmg,      &  !-- bbm
         &                  rn_dmg_max, rn_C0, rn_alrlx, rn_btrlx, rn_c_ref, rn_l_ref, ln_damaged_E,                &  !-- bbm
         &                  ln_tame_ini_ws, rn_half_tame,                                                           &  !-- bbm
         &                  ln_rhg_VP, nn_vp_nout, nn_vp_ninn, nn_vp_chkcvg                                            !-- vp
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, namdyn_rhg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namdyn_rhg in reference namelist' )
      READ  ( numnam_ice_cfg, namdyn_rhg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namdyn_rhg in configuration namelist' )
      IF(lwm) WRITE ( numoni, namdyn_rhg )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_rhg_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist : namdyn_rhg:'
         IF( ln_rhg_EVP .OR. ln_rhg_EAP ) THEN
         WRITE(numout,*) '      rheology EVP (icedyn_rhg_evp)                        ln_rhg_EVP    = ', ln_rhg_EVP
         WRITE(numout,*) '         use adaptive EVP (aEVP)                           ln_aEVP       = ', ln_aEVP
         WRITE(numout,*) '         creep limit                                       rn_creepl     = ', rn_creepl ! also used by vp
         WRITE(numout,*) '         eccentricity of the elliptical yield curve        rn_ecc        = ', rn_ecc    ! also used by vp
         WRITE(numout,*) '         number of iterations for subcycling               nn_nevp       = ', nn_nevp
         WRITE(numout,*) '         ratio of elastic timescale over ice time step     rn_relast     = ', rn_relast
         WRITE(numout,*) '         check convergence of rheology                     nn_rhg_chkcvg = ', nn_rhg_chkcvg
         ENDIF
         IF( ln_rhg_VP ) THEN
         WRITE(numout,*) '      rheology VP   (icedyn_rhg_VP)                        ln_rhg_VP     = ', ln_rhg_VP
         WRITE(numout,*) '         number of outer iterations                        nn_vp_nout    = ', nn_vp_nout
         WRITE(numout,*) '         number of inner iterations                        nn_vp_ninn    = ', nn_vp_ninn
         WRITE(numout,*) '         iteration step for convergence check              nn_vp_chkcvg  = ', nn_vp_chkcvg
         ENDIF
         IF( ln_rhg_BBM ) THEN
            IF(.NOT.ln_damage) CALL ctl_stop( 'ice_dyn_rhg_init: BBM rheology => set `ln_damage=.true` in `nampar`' )
            WRITE(numout,*) '    rheology BBM (icedyn_rhg_bbm)                          ln_rhg_BBM    = ', ln_rhg_BBM !#bbm
            IF(ln_MEB) WRITE(numout,*) '         will use the MEB rheology variant rather than pure BBM!' 
            WRITE(numout,*) '         max. compressive stress at the ref. scale [Pa]    rn_Nref       = ', rn_Nref  !#bbm
            WRITE(numout,*) '         elasticity of undamaged ice [Pa]                  rn_E0         = ', rn_E0  !#bbm
            WRITE(numout,*) '         viscosity of undamaged ice  [Pa.s]                rn_eta0       = ', rn_eta0  !#bbm
            IF(.NOT.ln_MEB ) WRITE(numout,*) '         compression factor "P" at play in "P_max"         rn_P0         = ', rn_P0  !#bbm            
            WRITE(numout,*) '         healing constant for damage                       rn_kth        = ', rn_kth  !#bbm
            WRITE(numout,*) '         number of iterations for subcycling               nn_nbbm       = ', nn_nbbm !#bbm
            WRITE(numout,*) '         advection of damage and stresses @T & @F          nn_d_adv  = ', nn_d_adv !#bbm
            IF( nn_d_adv==0 ) WRITE(numout,*) '           => no advection at all!' !#bbm
            IF( nn_d_adv==1 ) WRITE(numout,*) '           => advection of damage only' !#bbm
            IF( nn_d_adv >1 ) WRITE(numout,*) '           => advection of damage and stresses' !#bbm
            IF( nn_d_adv==3 ) WRITE(numout,*) '             ==> add "upper-convected time" terms in tensor advection' !#bbm
            WRITE(numout,*) '         cross-nudging coeff. for stress tensor            rn_crndg      = ', rn_crndg !#bbm
            IF(rn_crndg>0._wp) THEN
               WRITE(numout,*) '      => boost the CN at the coastline?             ln_boost_CN_coast = ', ln_boost_CN_coast
               IF(ln_boost_CN_coast) WRITE(numout,*) '                   ==>          rn_max_CN_coast = ', rn_max_CN_coast
               WRITE(numout,*) '      => boost the CN where damage is high?      ln_boost_CN_high_dmg = ', ln_boost_CN_high_dmg
               IF(ln_boost_CN_high_dmg) WRITE(numout,*) '                   ==>         rn_max_CN_dmg = ', rn_max_CN_dmg
            ENDIF
            WRITE(numout,*) '         ceiling value to cap damage with                  rn_dmg_max    = ', rn_dmg_max !#bbm
            WRITE(numout,*) '         compaction paramater (coeff. of exponential)      rn_C0         = ', rn_C0 !#bbm
            WRITE(numout,*) '         `alpha` of viscosity (dep. on `A` and `h`)        rn_alrlx      = ', rn_alrlx !#bbm
            WRITE(numout,*) '         `alpha`-like of Olason/Boutin in `lamdda`         rn_btrlx      = ', rn_btrlx !#bbm
            WRITE(numout,*) '         ice cohesion value at the lab scale               rn_c_ref      = ', rn_c_ref !#bbm
            WRITE(numout,*) '         scaling parameter for cohesion                    rn_l_ref      = ', rn_l_ref  !#bbm
            WRITE(numout,*) '         use damaged elasticity in MC test                 ln_damaged_E  = ', ln_damaged_E !#bbm
            WRITE(numout,*) '         gently increase wind stress from 0 if cold start  ln_tame_ini_ws= ', ln_tame_ini_ws !#bbm
            IF( ln_tame_ini_ws ) &
               & WRITE(numout,*) '         => delay (h) at which half of increase done       rn_half_tame  = ', rn_half_tame !#bbm
         END IF

         IF( ln_rhg_EVP ) THEN
            IF    ( nn_rhg_chkcvg == 0 ) THEN   ;   WRITE(numout,*) '         no check cvg'
            ELSEIF( nn_rhg_chkcvg == 1 ) THEN   ;   WRITE(numout,*) '         check cvg at the main time step'
            ELSEIF( nn_rhg_chkcvg == 2 ) THEN   ;   WRITE(numout,*) '         check cvg at both main and rheology time steps'
            ENDIF
         ENDIF
         WRITE(numout,*) '      rheology EAP (icedyn_rhg_eap)                        ln_rhg_EAP = ', ln_rhg_EAP
      ENDIF
      !
      
      !                             !== set the choice of ice advection ==!
      ioptio = 0
      IF( ln_rhg_EVP ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgEVP    ;   ENDIF
      IF( ln_rhg_EAP ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgEAP    ;   ENDIF
      IF( ln_rhg_VP  ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgVP     ;   ENDIF
      IF( ln_rhg_BBM ) THEN   ;   ioptio = ioptio + 1   ;   nice_rhg = np_rhgBBM    ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_dyn_rhg_init: choose one and only one ice rheology' )
      !
      IF( ln_rhg_EVP  )   CALL rhg_evp_rst( 'READ' )  !* read or initialize all required files
      IF( ln_rhg_EAP  )   CALL rhg_eap_rst( 'READ' )  !* read or initialize all required files
      
      IF( ln_rhg_BBM  ) THEN
         CALL ice_dyn_rhg_bbm_init() !* allocation of BBM-specific arrays (LB: needs to be done before reading restarts...)
         CALL rhg_bbm_rst( 'READ' )  !* read or initialize all required files
      END IF
      ! no restart for VP as there is no explicit time dependency in the equation
      !
   END SUBROUTINE ice_dyn_rhg_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_rhg
