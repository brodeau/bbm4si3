MODULE icedyn_adv
   !!======================================================================
   !!                       ***  MODULE icedyn_adv   ***
   !!   sea-ice: advection
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv   : advection of sea ice variables
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE sbc_oce , ONLY : nn_fsbc   ! frequency of sea-ice call
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   !
   USE icedyn_rhg,      ONLY : ln_rhg_EVP, ln_rhg_BBM   !#bbm
   USE icedyn_rhg_util, ONLY : strain_rate, clean_small_a_sgm !#bbm (advection of tensor)
   !
   USE icedyn_adv_pra ! sea-ice: advection scheme (Prather)
   USE icedyn_adv_umx ! sea-ice: advection scheme (ultimate-macho)
   !
   USE icedyn_adv_pra_bbm_t ! sea-ice: advection scheme (Prather) @T for `damage and stresses`  !bbm
   USE icedyn_adv_pra_bbm_f ! sea-ice: advection scheme (Prather) @F for `damage and stresses`  !bbm
   !
   USE icedyn_rhg_util, ONLY : cap_damage
   !
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE timing         ! Timing
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv        ! called by icestp
   PUBLIC   ice_dyn_adv_init   ! called by icedyn

   INTEGER ::              nice_adv   ! choice of the type of advection scheme
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_advPRA = 1   ! Prather scheme
   INTEGER, PARAMETER ::   np_advUMx = 2   ! Ultimate-Macho scheme
   !
   ! ** namelist (namdyn_adv) **
   INTEGER         ::   nn_UMx       ! order of the UMx advection scheme
   !
   LOGICAL, SAVE :: l_advct_stresses, l_advct_oldroyd !#bbm

   ! Transformation of stresses to something positive and not too large!!!
   REAL(wp), PARAMETER :: r_sgm_sf    = 1.E-4_wp        ! scale-factor from Pa to 10xM.Pa
   REAL(wp), PARAMETER :: r_1_sgm_sf  = 1._wp/r_sgm_sf
   REAL(wp), PARAMETER :: r_mltpl_s12 = 5.0_wp
   REAL(wp), PARAMETER :: r_mltpl_skk = 1.0_wp


   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_adv.F90 13472 2020-09-16 13:05:19Z smasson $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_dyn_adv ***
      !!
      !! ** purpose : advection of sea ice
      !!
      !! ** method  : One can choose between
      !!     a) an Ultimate-Macho scheme (with order defined by nn_UMx) => ln_adv_UMx
      !!     b) and a second order Prather scheme => ln_adv_Pra
      !!
      !! ** action :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      REAL(wp), DIMENSION(jpi,jpj) :: zdudx, zdvdy, zdudy, zdvdx, zdiv  !#bbm
      REAL(wp), DIMENSION(jpi,jpj) :: zdmg_pos, zs11_pos, zs22_pos, zs12_pos  !#bbm
      REAL(wp)                     :: zm1t, zm1f,zm2t, zm2f
      REAL(wp)                     :: z_skk_ao, z_s12_ao
      !!---------------------------------------------------------------------
      !
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_adv')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv: sea-ice advection'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      !---------------!
      !== Advection ==!
      !---------------!
      SELECT CASE( nice_adv )
      !                                !-----------------------!
      CASE( np_advUMx )                ! ULTIMATE-MACHO scheme !
         !                             !-----------------------!
         IF(lwp) WRITE(numout,'("  *** advects GENERIC fields @T with UMX, order = ",i1," kt=",i5.5)') nn_UMx, kt
         CALL ice_dyn_adv_umx( nn_UMx, kt, u_ice, v_ice, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i )
         !                             !-----------------------!
      CASE( np_advPRA )                ! PRATHER scheme        !
         !                             !-----------------------!
         IF(lwp) WRITE(numout,'("  *** advects GENERIC fields @T with Prather, kt=",i5.5)') kt
         CALL ice_dyn_adv_pra(         kt, u_ice, v_ice, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i )
      END SELECT


      IF( ln_rhg_BBM ) THEN

         IF( nn_d_adv >= 1 ) THEN

            IF( l_advct_stresses ) THEN
               !! 1st, find the maximum values to pick an `add_offset` !
               !! for `skk` we nultiply it by `-1` so negative values will be smaller (normally: `HUGE_NEG_VAL < skk < SMALL_POS_VAL` )
               zm1t = MAXVAL( sgm11t(:,:)*r_sgm_sf*xmskt(:,:) )
               zm2t = MAXVAL( sgm22t(:,:)*r_sgm_sf*xmskt(:,:) )
               zm1f = MAXVAL( sgm11f(:,:)*r_sgm_sf*xmskf(:,:) )
               zm2f = MAXVAL( sgm22f(:,:)*r_sgm_sf*xmskf(:,:) )
               z_skk_ao = MAX( MAX(zm1t,zm2t), MAX(zm1f,zm2f) )
               z_skk_ao = MAX( z_skk_ao, 1.E-3_wp ) ! safety...
               CALL mpp_max( 'icedyn_adv', z_skk_ao )
               z_skk_ao = r_mltpl_skk * REAL( CEILING( z_skk_ao/r_mltpl_skk) , wp ) ! we want it to be the multiple of `r_mltpl_skk` and just above...
               z_skk_ao = MAX( z_skk_ao, 2._wp )
               !
               zm1t = MINVAL( sgm12t(:,:)*r_sgm_sf*xmskt(:,:) )
               zm1f = MINVAL( sgm12f(:,:)*r_sgm_sf*xmskf(:,:) )
               z_s12_ao = -1._wp * MIN( zm1t, zm1f ) ! => positive
               z_s12_ao = MAX( z_s12_ao, 1.E-3_wp ) ! safety...
               CALL mpp_max( 'icedyn_adv', z_s12_ao )
               z_s12_ao = r_mltpl_s12 * REAL( CEILING( z_s12_ao/r_mltpl_s12) , wp ) ! be a multiple of `r_mltpl_s12` and just above it...
               z_s12_ao = MAX( z_s12_ao, 5._wp )
            ENDIF

            !! Advect at T points with u@U, v@V velocities
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! Transform before advection:
            zdmg_pos =      (1._wp - dmgt)           * xmskt  ! cleaner to advect `1-d` rather than `d`...
            zs11_pos = (-sgm11t*r_sgm_sf + z_skk_ao) * xmskt
            zs22_pos = (-sgm22t*r_sgm_sf + z_skk_ao) * xmskt
            zs12_pos = ( sgm12t*r_sgm_sf + z_s12_ao) * xmskt
            !
            IF( ANY(zs11_pos<0._wp) .OR. ANY(zs22_pos<0._wp) ) THEN
               WRITE(numout,*) ' *** Min val for `zs11_pos`:', MINVAL(zs11_pos)
               WRITE(numout,*) ' *** Min val for `zs22_pos`:', MINVAL(zs22_pos)
               CALL ctl_stop( 'ice_dyn_adv: a "transformed" sigma_ii @T still has negative value(s)!!!')
            ENDIF
            IF( ANY(zs12_pos<0._wp) ) THEN
               WRITE(numout,*) ' *** Min val for `zs12_pos`:', MINVAL(zs12_pos)
               CALL ctl_stop( 'ice_dyn_adv: a "transformed" sigma_12 @T still has negative value(s)!!!')
            ENDIF

            IF( l_advct_stresses ) THEN
               IF(lwp) WRITE(numout,'("  *** advects damage & stresses @T with Prather, kt=",i5.5)') kt
               IF(lwp) WRITE(numout,'("      ==> offset for `skk` and `s12`: ",f6.2,", ",f6.2)') z_skk_ao, z_s12_ao
               CALL ice_dyn_adv_pra_t_d( kt, u_ice, v_ice,  zdmg_pos, pdd1=zs12_pos, pdd2=zs11_pos, pdd3=zs22_pos )
            ELSE
               IF(lwp) WRITE(numout,'("  *** advects only damage @T with Prather, kt=",i5.5)') kt
               CALL ice_dyn_adv_pra_t_d( kt, u_ice, v_ice,  zdmg_pos )
            ENDIF
            !
            ! Fall back after advection:
            dmgt   = (1._wp - zdmg_pos)
            CALL cap_damage( 'T', 'ice_dyn_adv', dmgt )
            sgm11t = (z_skk_ao - zs11_pos) * r_1_sgm_sf * xmskt
            sgm22t = (z_skk_ao - zs22_pos) * r_1_sgm_sf * xmskt
            sgm12t = (zs12_pos - z_s12_ao) * r_1_sgm_sf * xmskt


            !! Advect at F points with u@V, v@U velocities
            !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! Transform before advection:
            zdmg_pos = (1._wp - dmgf)                * xmskf
            zs11_pos = (-sgm11f*r_sgm_sf + z_skk_ao) * xmskf
            zs22_pos = (-sgm22f*r_sgm_sf + z_skk_ao) * xmskf
            zs12_pos = ( sgm12f*r_sgm_sf + z_s12_ao) * xmskf
            !
            IF( ANY(zs11_pos<0._wp) .OR. ANY(zs22_pos<0._wp) ) THEN
               WRITE(numout,*) ' *** Min val for `zs11_pos`:', MINVAL(zs11_pos)
               WRITE(numout,*) ' *** Min val for `zs22_pos`:', MINVAL(zs22_pos)
               CALL ctl_stop( 'ice_dyn_adv: a "transformed" sigma_ii @F still has negative value(s)!!!')
            ENDIF
            IF( ANY(zs12_pos<0._wp) ) THEN
               WRITE(numout,*) ' *** Min val for `zs12_pos`:', MINVAL(zs12_pos)
               CALL ctl_stop( 'ice_dyn_adv: a "transformed" sigma_12 @F still has negative value(s)!!!')
            ENDIF
            !
            IF( l_advct_stresses ) THEN
               IF(lwp) WRITE(numout,'("  *** advects damage & stresses @F with Prather, kt=",i5.5)') kt
               CALL ice_dyn_adv_pra_f_d( kt, uVice, vUice,  zdmg_pos, pdd1=zs12_pos, pdd2=zs11_pos, pdd3=zs22_pos )
            ELSE
               IF(lwp) WRITE(numout,'("  *** advects only damage @F with Prather, kt=",i5.5)') kt
               CALL ice_dyn_adv_pra_f_d( kt, uVice, vUice,  zdmg_pos )
            ENDIF

            ! Fall back after advection:
            dmgf   = (1._wp - zdmg_pos)
            CALL cap_damage( 'F', 'ice_dyn_adv', dmgf )
            sgm11f = (z_skk_ao - zs11_pos) * r_1_sgm_sf * xmskf
            sgm22f = (z_skk_ao - zs22_pos) * r_1_sgm_sf * xmskf
            sgm12f = (zs12_pos - z_s12_ao) * r_1_sgm_sf * xmskf

            IF( l_advct_oldroyd ) THEN
               !!
               CALL strain_rate( 'T', u_ice, v_ice, uVice, vUice, &
                  &              r1_e1e2t, e2u, e1v, r1_e2u, r1_e1v, e1t*e1t, e2t*e2t, tmask(:,:,1), &
                  &              zdudx, zdvdy, zdmg_pos, lblnk=.TRUE., pdudy=zdudy, pdvdx=zdvdx, pdiv=zdiv )
               !!                         =>  `zdmg_pos` used for receiving the shear, which is not used anyway...
               IF(nn_d_adv==3) THEN
                  IF(lwp) WRITE(numout,'("  *** Going for LOWER-ADVECTED advection term, kt=",i5.5)') kt
                  CALL add_lower_advected( zdudx, zdudy, zdvdx, zdvdy, zdiv, xmskt,  sgm11t, sgm22t, sgm12t )
               ENDIF
               IF(nn_d_adv==4) THEN
                  IF(lwp) WRITE(numout,'("  *** Going for UPPER-ADVECTED advection term, kt=",i5.5)') kt
                  CALL add_upper_advected( zdudx, zdudy, zdvdx, zdvdy, zdiv, xmskt,  sgm11t, sgm22t, sgm12t )
               ENDIF
               !!
               CALL strain_rate( 'F', uVice, vUice, u_ice, v_ice, &
                  &              r1_e1e2f, e2v, e1u, r1_e2v, r1_e1u, e1f*e1f, e2f*e2f, fmask(:,:,1), &
                  &              zdudx, zdvdy, zdmg_pos,  lblnk=.TRUE., pdudy=zdudy, pdvdx=zdvdx, pdiv=zdiv )
               !!                         =>  `zdmg_pos` used for receiving the shear, which is not used anyway...
               IF(nn_d_adv==3) THEN
                  CALL add_lower_advected( zdudx, zdudy, zdvdx, zdvdy, zdiv, xmskf,  sgm11f, sgm22f, sgm12f )
               ENDIF
               IF(nn_d_adv==4) THEN
                  CALL add_upper_advected( zdudx, zdudy, zdvdx, zdvdy, zdiv, xmskf,  sgm11f, sgm22f, sgm12f )
               ENDIF

            ELSE
               IF(lwp) WRITE(numout,'("  *** No Upper- or Lower-convected advection term used! kt=",i5.5)') kt
            END IF !IF( l_advct_oldroyd )

            CALL clean_small_a_sgm( 'T', at_i, a_f,  sgm11t, sgm22t, sgm12f )
            CALL clean_small_a_sgm( 'F', at_i, a_f,  sgm11f, sgm22f, sgm12t )

         END IF !IF( nn_d_adv >= 1 )
      END IF !IF( ln_rhg_BBM )

      !------------
      ! diagnostics
      !------------
      diag_trp_ei(:,:) = SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice
      diag_trp_es(:,:) = SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
      diag_trp_sv(:,:) = SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice
      diag_trp_vi(:,:) = SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice
      diag_trp_vs(:,:) = SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice
      IF( iom_use('icemtrp') )   CALL iom_put( 'icemtrp' ,  diag_trp_vi * rhoi          )   ! ice mass transport
      IF( iom_use('snwmtrp') )   CALL iom_put( 'snwmtrp' ,  diag_trp_vs * rhos          )   ! snw mass transport
      IF( iom_use('salmtrp') )   CALL iom_put( 'salmtrp' ,  diag_trp_sv * rhoi * 1.e-03 )   ! salt mass transport (kg/m2/s)
      IF( iom_use('dihctrp') )   CALL iom_put( 'dihctrp' , -diag_trp_ei                 )   ! advected ice heat content (W/m2)
      IF( iom_use('dshctrp') )   CALL iom_put( 'dshctrp' , -diag_trp_es                 )   ! advected snw heat content (W/m2)

      ! controls
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn & trp - ')                           ! prints
      IF( ln_timing    )   CALL timing_stop ('icedyn_adv')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_adv


   SUBROUTINE ice_dyn_adv_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_init  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice
      !!                dynamics
      !!
      !! ** Method  :   Read the namdyn_adv namelist and check the ice-dynamic
      !!                parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_adv
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_adv/ ln_adv_Pra, ln_adv_UMx, nn_UMx
      !!-------------------------------------------------------------------
      !
      READ  ( numnam_ice_ref, namdyn_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_adv in reference namelist' )
      READ  ( numnam_ice_cfg, namdyn_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_adv in configuration namelist' )
      IF(lwm) WRITE( numoni, namdyn_adv )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_adv:'
         WRITE(numout,*) '      type of advection scheme (Prather)             ln_adv_Pra = ', ln_adv_Pra
         WRITE(numout,*) '      type of advection scheme (Ulimate-Macho)       ln_adv_UMx = ', ln_adv_UMx
         WRITE(numout,*) '         order of the Ultimate-Macho scheme          nn_UMx     = ', nn_UMx
      ENDIF
      !
      !                             !== set the choice of ice advection ==!
      ioptio = 0
      IF( ln_adv_Pra ) THEN
         ioptio = ioptio + 1
         nice_adv = np_advPRA
      ENDIF
      IF( ln_adv_UMx ) THEN
         ioptio = ioptio + 1
         nice_adv = np_advUMx
      ENDIF
      IF( ioptio /= 1 ) CALL ctl_stop( 'ice_dyn_adv_init: choose one and only one ice adv. scheme (ln_adv_Pra or ln_adv_UMx)' )
      !
      IF( ln_adv_Pra )  CALL adv_pra_init  !* read or initialize all required files
      !
      IF( ln_rhg_BBM ) THEN
         l_advct_stresses = (nn_d_adv >= 2) ! advect non only damage but also components of stress tensors
         l_advct_oldroyd  = (nn_d_adv >= 3) ! add the terms for advection of tensor ("upper-convected time" aka "Oldroyd" derivatives)
         CALL adv_pra_t_d_init( )
         CALL adv_pra_f_d_init( )
      END IF
      !
   END SUBROUTINE ice_dyn_adv_init


   SUBROUTINE add_lower_advected( pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk, ps11, ps22, ps12 )
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11, ps22, ps12
      !! Lower-convected version (sign is inversed because moved from LHS to RHS):
      ps11(:,:) = ps11(:,:) - 2._wp*rDt_ice*( pdudx(:,:)*ps11(:,:) + pdvdx(:,:)*ps12(:,:) )*zmsk(:,:)
      ps22(:,:) = ps22(:,:) - 2._wp*rDt_ice*( pdvdy(:,:)*ps22(:,:) + pdudy(:,:)*ps12(:,:) )*zmsk(:,:)
      ps12(:,:) = ps12(:,:) - rDt_ice*( pdiv(:,:)*ps12(:,:) + pdudy(:,:)*ps11(:,:) + pdvdx(:,:)*ps22(:,:) )*zmsk(:,:)
   END SUBROUTINE add_lower_advected

   SUBROUTINE add_upper_advected( pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk, ps11, ps22, ps12 )
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pdudx, pdudy, pdvdx, pdvdy, pdiv, zmsk
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11, ps22, ps12
      !! Upper-convected version (sign is inversed because moved from LHS to RHS):
      ps11(:,:) = ps11(:,:) + 2._wp*rDt_ice*( pdudx(:,:)*ps11(:,:) + pdudy(:,:)*ps12(:,:) )*zmsk(:,:)
      ps22(:,:) = ps22(:,:) + 2._wp*rDt_ice*( pdvdy(:,:)*ps22(:,:) + pdvdx(:,:)*ps12(:,:) )*zmsk(:,:)
      ps12(:,:) = ps12(:,:) + rDt_ice*( pdiv(:,:)*ps12(:,:) + pdudy(:,:)*ps22(:,:) + pdvdx(:,:)*ps11(:,:) )*zmsk(:,:)
   END SUBROUTINE add_upper_advected


#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv
