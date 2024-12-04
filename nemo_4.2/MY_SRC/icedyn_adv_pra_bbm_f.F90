MODULE icedyn_adv_pra_bbm_f

   !!BBM => advects damage and stress tensor components @T-point

   !!======================================================================
   !!                       ***  MODULE icedyn_adv_pra_bbm_f   ***
   !!   sea-ice : advection => Prather scheme
   !!======================================================================
   !! History :       !  2008-03  (M. Vancoppenolle) original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!            4.2  !  2024     (L. Brodeau)       BBM
   !!--------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_pra_f_d : advection of sea ice using Prather scheme
   !!   adv_x, adv_y    : Prather scheme applied in i- and j-direction, resp.
   !!   adv_pra_f_d_init    : initialisation of the Prather scheme
   !!   adv_pra_f_d_rst     : read/write Prather field in ice restart file, or initialized to zero
   !!----------------------------------------------------------------------
   USE phycst         ! physical constant
   USE dom_oce        ! ocean domain
   USE ice            ! sea-ice variables
   USE sbc_oce , ONLY : nn_fsbc   ! frequency of sea-ice call
   USE icevar         ! sea-ice: operations
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv_pra_f_d   ! called by icedyn_adv
   PUBLIC   adv_pra_f_d_init      ! called by icedyn_adv

   ! Moments for advection
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sx1md, sy1md, sxx1md, syy1md, sxy1md   ! ice damage
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sxdd1, sydd1, sxxdd1, syydd1, sxydd1   ! ice damage
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sxdd2, sydd2, sxxdd2, syydd2, sxydd2   ! ice damage
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sxdd3, sydd3, sxxdd3, syydd3, sxydd3   ! ice damage

   REAL(wp), PARAMETER :: rr_scl_fct = 1.E-6_wp
   REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_adv_pra_bbm_f.F90 14026 2020-12-03 08:48:10Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_pra_f_d( kt, pu_ice, pv_ice,  p1md, pdd1, pdd2, pdd3 )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_pra_f_d  **
      !!
      !! ** BRODEAU, 2022 => Advect `damage` at F-points only
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!                but also their quadratic forms. The method preserves
      !!                tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!----------------------------------------------------------------------
      INTEGER,                            INTENT(in   ) :: kt         ! current time step
      REAL(wp), DIMENSION(:,:),           INTENT(in   ) :: pu_ice     ! ice i-velocity
      REAL(wp), DIMENSION(:,:),           INTENT(in   ) :: pv_ice     ! ice j-velocity
      REAL(wp), DIMENSION(:,:),           INTENT(inout) :: p1md       ! (1 - damage) @ F
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(inout) :: pdd1, pdd2, pdd3 ! stresses @ F
      !
      LOGICAL  :: l_advect_sigma = .FALSE.
      INTEGER  ::   jt      ! dummy loop indices
      INTEGER  ::   icycle                  ! number of sub-timestep for the advection
      REAL(wp) ::   zdt, z1_dt              !   -      -
      REAL(wp), DIMENSION(1)                  ::   zcflprv, zcflnow   ! for global communication
      REAL(wp), DIMENSION(jpi,jpj)            ::   zudy, zvdx
      REAL(wp), DIMENSION(jpi,jpj)            ::   z2d
      REAL(wp), DIMENSION(jpi,jpj)            ::   zarea
      REAL(wp), DIMENSION(jpi,jpj)            ::   z0di, z0d1, z0d2, z0d3
      INTEGER :: ji, jj
      CHARACTER(len=32) :: cstr
      CHARACTER(len=19), PARAMETER :: crtnm = 'ice_dyn_adv_pra_f_d'
      CHARACTER(len=1),  PARAMETER :: cgrt    = 'F'
      !!----------------------------------------------------------------------
      !
      l_advect_sigma = ( PRESENT(pdd1) .AND. PRESENT(pdd2) .AND. PRESENT(pdd3) )
      !
      cstr = 'Prather advection scheme for BBM (damage only)'
      IF(l_advect_sigma) cstr = 'Prather advection scheme for BBM (damage and stresses )'
      IF( kt == nit000 .AND. lwp ) WRITE(numout,*) '-- '//crtnm//': '//TRIM(cstr)//' at '//cgrt//'-points'
      !
      ! --- If ice drift is too fast, use  subtime steps for advection (CFL test for stability) --- !
      !        Note: the advection split is applied at the next time-step in order to avoid blocking global comm.
      !              this should not affect too much the stability
      zcflnow(1) =                  MAXVAL( ABS( pu_ice(:,:) ) * rDt_ice * r1_e1v(:,:) )
      zcflnow(1) = MAX( zcflnow(1), MAXVAL( ABS( pv_ice(:,:) ) * rDt_ice * r1_e2u(:,:) ) )

      ! non-blocking global communication send zcflnow and receive zcflprv
      CALL mpp_delay_max( crtnm, 'cflice', zcflnow(:), zcflprv(:), kt == nitend - nn_fsbc + 1 )

      IF( zcflprv(1) > .5 ) THEN
         icycle = 2
      ELSE
         icycle = 1
      ENDIF
      zdt = rDt_ice / REAL(icycle)
      z1_dt = 1._wp / zdt

      ! --- transport --- !
      zudy(:,:) = pu_ice(:,:) * e2v(:,:) * rr_scl_fct    ! `rr_scl_fct` because we work with km^2 rather than m^2 !!!
      zvdx(:,:) = pv_ice(:,:) * e1u(:,:) * rr_scl_fct    !       "                     "

      DO jt = 1, icycle

         z2d(:,:) = e1e2f(:,:) * rr_scl_fct * xmskf(:,:)  ! cell area in km^2

         ! --- transported fields --- !
         zarea(:,:) = z2d(:,:)
         z0di (:,:) = p1md (:,:) * z2d(:,:)        ! Damage content !#bbm
         IF(l_advect_sigma) THEN
            z0d1 (:,:) = pdd1 (:,:) * z2d(:,:)        ! Damage content !#bbm
            z0d2 (:,:) = pdd2 (:,:) * z2d(:,:)        ! Damage content !#bbm
            z0d3 (:,:) = pdd3 (:,:) * z2d(:,:)        ! Damage content !#bbm
         END IF
         !
         !
         !                                                                  !--------------------------------------------!
         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !                                                               !--------------------------------------------!
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0di, sx1md, sxx1md, sy1md, syy1md, sxy1md )
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0di, sx1md, sxx1md, sy1md, syy1md, sxy1md )
            IF(l_advect_sigma) THEN
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0d1, sxdd1, sxxdd1, sydd1, syydd1, sxydd1 )
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0d1, sxdd1, sxxdd1, sydd1, syydd1, sxydd1 )
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0d2, sxdd2, sxxdd2, sydd2, syydd2, sxydd2 )
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0d2, sxdd2, sxxdd2, sydd2, syydd2, sxydd2 )
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0d3, sxdd3, sxxdd3, sydd3, syydd3, sxydd3 )
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0d3, sxdd3, sxxdd3, sydd3, syydd3, sxydd3 )
            ENDIF
            !
            !
            !                                                               !--------------------------------------------!
         ELSE                                                               !== even ice time step:  adv_y then adv_x  ==!
            !                                                               !--------------------------------------------!
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0di, sx1md, sxx1md, sy1md, syy1md, sxy1md )
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0di, sx1md, sxx1md, sy1md, syy1md, sxy1md )
            IF(l_advect_sigma) THEN
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0d1, sxdd1, sxxdd1, sydd1, syydd1, sxydd1 )
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0d1, sxdd1, sxxdd1, sydd1, syydd1, sxydd1 )
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0d2, sxdd2, sxxdd2, sydd2, syydd2, sxydd2 )
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0d2, sxdd2, sxxdd2, sydd2, syydd2, sxydd2 )
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0d3, sxdd3, sxxdd3, sydd3, syydd3, sxydd3 )
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0d3, sxdd3, sxxdd3, sydd3, syydd3, sxydd3 )
            ENDIF
            !
         ENDIF

         ! --- Lateral boundary conditions --- !
         !     caution: for gradients (sx and sy) the sign changes
         IF(l_advect_sigma) THEN
            CALL lbc_lnk( crtnm,   z0di, cgrt, 1._wp, sx1md , cgrt, -1._wp, sy1md , cgrt, -1._wp  &
               &                , sxx1md, cgrt, 1._wp, syy1md, cgrt,  1._wp, sxy1md, cgrt,  1._wp  &
               &                ,   z0d1, cgrt, 1._wp, sxdd1 , cgrt, -1._wp, sydd1 , cgrt, -1._wp  &
               &                , sxxdd1, cgrt, 1._wp, syydd1, cgrt,  1._wp, sxydd1, cgrt,  1._wp  &
               &                ,   z0d2, cgrt, 1._wp, sxdd2 , cgrt, -1._wp, sydd2 , cgrt, -1._wp  &
               &                , sxxdd2, cgrt, 1._wp, syydd2, cgrt,  1._wp, sxydd2, cgrt,  1._wp  &
               &                ,   z0d3, cgrt, 1._wp, sxdd3 , cgrt, -1._wp, sydd3 , cgrt, -1._wp  &
               &                , sxxdd3, cgrt, 1._wp, syydd3, cgrt,  1._wp, sxydd3, cgrt,  1._wp  )
         ELSE
            CALL lbc_lnk( crtnm,   z0di, cgrt, 1._wp, sx1md , cgrt, -1._wp, sy1md , cgrt, -1._wp  &
               &                , sxx1md, cgrt, 1._wp, syy1md, cgrt,  1._wp, sxy1md, cgrt,  1._wp  )
         ENDIF

         ! --- Recover the properties from their contents --- !
         z2d(:,:) = r1_scl_fct * r1_e1e2f(:,:) * xmskf(:,:)
         p1md (:,:) = z0di (:,:) * z2d(:,:) !#bbm
         IF(l_advect_sigma) THEN
            pdd1 (:,:) = z0d1 (:,:) * z2d(:,:) !#bbm
            pdd2 (:,:) = z0d2 (:,:) * z2d(:,:) !#bbm
            pdd3 (:,:) = z0d3 (:,:) * z2d(:,:) !#bbm
         END IF

         p1md(:,:) = MIN( MAX( p1md(:,:), 1._wp - rn_dmg_max ) , 1._wp ) !! `p1md` is `1-damage` !!!

      END DO !DO jt = 1, icycle
      !
      IF( lrst_ice )   CALL adv_pra_f_d_rst( 'WRITE', kt )   !* write Prather fields in the restart file
      !
   END SUBROUTINE ice_dyn_adv_pra_f_d


   SUBROUTINE adv_x( pdt, puf, pcrh, pe1e2, psm, ps0, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_x  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on x axis
      !!----------------------------------------------------------------------
      REAL(wp)                , INTENT(in   ) ::   pdt                ! the time step
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   puf                ! i-direction ice velocity at V-point [m/s]
      REAL(wp)                , INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pe1e2
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ps0                ! field to be advected
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!
      INTEGER  ::   ji, jj                               ! dummy loop indices
      INTEGER  ::   kj0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zU
      REAL(wp), DIMENSION(jpi,jpj) ::   zf0, zfx , zfy , zbet   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zfm, zfxx, zfyy, zfxy   !  -      -
      REAL(wp), DIMENSION(jpi,jpj) ::   zalg, zalg1, zalg1q     !  -      -
      !---------------------------------------------------------------------
      ! in order to avoid lbc_lnk (communications):
      !    jj loop must be 1:jpj   if adv_x is called first
      !                and 2:jpj-1 if adv_x is called second
      kj0 = NINT(pcrh)
      !
      ! Limitation of moments.
      DO_2D( 1, 0, kj0, kj0 )
            !! Only `ji+1` needed
            !! => [0:jpi-1] for `ji` !
            !
            zmask = xmskf(ji,jj)
            zU    = puf(ji+1,jj)
            !
            ze1e2 = pe1e2(ji,jj)
            zpsm  = psm (ji,jj) ! optimization
            zps0  = ps0 (ji,jj)
            zpsx  = psx (ji,jj)
            zpsxx = psxx(ji,jj)
            zpsy  = psy (ji,jj)
            zpsyy = psyy(ji,jj)
            zpsxy = psxy(ji,jj)

            !  Initialize volumes of boxes (=area if adv_x first called, =psm otherwise)
            zpsm = MAX( pcrh * ze1e2  +  ( 1.0 - pcrh ) * zpsm , epsi20 )
            !
            zslpmax = MAX( 0._wp, zps0 )
            zs1max  = 1.5 * zslpmax
            zs1new  = MIN( zs1max, MAX( -zs1max, zpsx ) )
            zs2new  = MIN( 2.0 * zslpmax - 0.3334 * ABS( zs1new ), MAX( ABS( zs1new ) - zslpmax, zpsxx ) )
            rswitch = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zslpmax) ) ) * zmask ! Case of empty boxes & Apply mask

            zps0  = zslpmax
            zpsx  = zs1new * rswitch
            zpsxx = zs2new * rswitch
            zpsy  = zpsy   * rswitch
            zpsyy = zpsyy  * rswitch
            zpsxy = MIN( zslpmax, MAX( -zslpmax, zpsxy ) ) * rswitch

            !  Calculate fluxes and moments between boxes i<-->i+1
            !                                !  Flux from i to i+1 WHEN u GT 0
            !! #F-centric:
            !! When on T point:  flux from T_i to T_i+1 is done via u@U_i
            !! But when F point: flux from F_i to F_i+1 is done via u@V_i+1 !!!
            zbet(ji,jj)  =  MAX( 0._wp, SIGN( 1._wp, zU ) )
            zalf         =  MAX( 0._wp, zU ) * pdt / zpsm
            zalfq        =  zalf * zalf
            zalf1        =  1.0 - zalf
            zalf1q       =  zalf1 * zalf1
            !
            zfm (ji,jj)  =  zalf  * zpsm
            zf0 (ji,jj)  =  zalf  * ( zps0  + zalf1 * ( zpsx + (zalf1 - zalf) * zpsxx ) )
            zfx (ji,jj)  =  zalfq * ( zpsx  + 3.0 * zalf1 * zpsxx )
            zfxx(ji,jj)  =  zalf  * zalfq * zpsxx
            zfy (ji,jj)  =  zalf  * ( zpsy  + zalf1 * zpsxy )
            zfxy(ji,jj)  =  zalfq * zpsxy
            zfyy(ji,jj)  =  zalf  * zpsyy
            !                                !  Readjust moments remaining in the box.
            zpsm  =  zpsm  - zfm(ji,jj)
            zps0  =  zps0  - zf0(ji,jj)
            zpsx  =  zalf1q * ( zpsx - 3.0 * zalf * zpsxx )
            zpsxx =  zalf1  * zalf1q * zpsxx
            zpsy  =  zpsy  - zfy (ji,jj)
            zpsyy =  zpsyy - zfyy(ji,jj)
            zpsxy =  zalf1q * zpsxy
            !
            psm (ji,jj) = zpsm ! optimization
            ps0 (ji,jj) = zps0
            psx (ji,jj) = zpsx
            psxx(ji,jj) = zpsxx
            psy (ji,jj) = zpsy
            psyy(ji,jj) = zpsyy
            psxy(ji,jj) = zpsxy
            !
      END_2D
      ! ji+ ji- jj+ jj-

      !LB: because of the F-centric world, need to link here => to fill the i+1
      CALL lbc_lnk( 'adv_x', zfm,'F', 1._wp,  zf0,'F',1._wp, zfxy,'F',1._wp, &
         &                   zfx,'F',-1._wp, zfxx,'F',1._wp, zbet,'F',1._wp, &
         &                   zfy,'F',-1._wp, zfyy,'F',1._wp,                 &
         &                   psm,'F', 1._wp,  ps0,'F',1._wp, psxy,'F',1._wp, &
         &                   psx,'F',-1._wp, psxx,'F',1._wp,                 &
         &                   psy,'F',-1._wp, psyy,'F',1._wp  )
      !! Fix due to the use of lbc_lnk, to avoid division by 0 later on...
      IF( ANY(psm<0._wp) ) CALL ctl_stop('STOP', 'adv_x() [icedyn_adv_pra_bbm_f.F90] : `psm<0` somewhere')
      WHERE( psm < epsi20 ) psm = epsi20

      DO_2D( 1, 0, kj0, kj0 )
            !! Only `ji+1` needed
            !! => [0:jpi-1] for `ji` !
            zU = puf(ji+1,jj)
            !                                !  Flux from i+1 to i when u LT 0.
            zalf          = MAX( 0._wp, -zU ) * pdt / psm(ji+1,jj)
            zalg  (ji,jj) = zalf
            zalfq         = zalf * zalf
            zalf1         = 1.0 - zalf
            zalg1 (ji,jj) = zalf1
            zalf1q        = zalf1 * zalf1
            zalg1q(ji,jj) = zalf1q
            !
            zfm   (ji,jj) = zfm (ji,jj) + zalf  *    psm (ji+1,jj)
            zf0   (ji,jj) = zf0 (ji,jj) + zalf  * (  ps0 (ji+1,jj) &
               &          - zalf1 * ( psx(ji+1,jj) - (zalf1 - zalf ) * psxx(ji+1,jj) ) )
            zfx   (ji,jj) = zfx (ji,jj) + zalfq * (  psx (ji+1,jj) - 3.0 * zalf1 * psxx(ji+1,jj) )
            zfxx  (ji,jj) = zfxx(ji,jj) + zalf  *    psxx(ji+1,jj) * zalfq
            zfy   (ji,jj) = zfy (ji,jj) + zalf  * (  psy (ji+1,jj) - zalf1 * psxy(ji+1,jj) )
            zfxy  (ji,jj) = zfxy(ji,jj) + zalfq *    psxy(ji+1,jj)
            zfyy  (ji,jj) = zfyy(ji,jj) + zalf  *    psyy(ji+1,jj)
      END_2D
      ! ji+ ji- jj+ jj-

      DO_2D( 0, 1, kj0, kj0 )
            !! Only `ji-1` needed
            !! => [2:jpi] for `ji` !
            !                                !  Readjust moments remaining in the box.
            zbt  =       zbet(ji-1,jj)
            zbt1 = 1.0 - zbet(ji-1,jj)
            !
            zpsm  = psm (ji,jj) ! optimization
            zps0  = ps0 (ji,jj)
            zpsx  = psx (ji,jj)
            zpsxx = psxx(ji,jj)
            zpsy  = psy (ji,jj)
            zpsyy = psyy(ji,jj)
            zpsxy = psxy(ji,jj)
            !
            zpsm  = zbt * zpsm + zbt1 * ( zpsm - zfm(ji-1,jj) )
            zps0  = zbt * zps0 + zbt1 * ( zps0 - zf0(ji-1,jj) )
            zpsx  = zalg1q(ji-1,jj) * ( zpsx + 3.0 * zalg(ji-1,jj) * zpsxx )
            zpsxx = zalg1 (ji-1,jj) * zalg1q(ji-1,jj) * zpsxx
            zpsy  = zbt * zpsy  + zbt1 * ( zpsy  - zfy (ji-1,jj) )
            zpsyy = zbt * zpsyy + zbt1 * ( zpsyy - zfyy(ji-1,jj) )
            zpsxy = zalg1q(ji-1,jj) * zpsxy

            !   Put the temporary moments into appropriate neighboring boxes.
            !                                !   Flux from i to i+1 IF u GT 0.
            zbt   =       zbet(ji-1,jj)
            zbt1  = 1.0 - zbet(ji-1,jj)
            zpsm  = zbt * ( zpsm + zfm(ji-1,jj) ) + zbt1 * zpsm
            zalf  = zbt * zfm(ji-1,jj) / zpsm
            zalf1 = 1.0 - zalf
            ztemp = zalf * zps0 - zalf1 * zf0(ji-1,jj)
            !
            zps0  =  zbt  * ( zps0 + zf0(ji-1,jj) ) + zbt1 * zps0
            zpsx  =  zbt  * ( zalf * zfx(ji-1,jj) + zalf1 * zpsx + 3.0 * ztemp ) + zbt1 * zpsx
            zpsxx =  zbt  * ( zalf * zalf * zfxx(ji-1,jj) + zalf1 * zalf1 * zpsxx                            &
               &            + 5.0 * ( zalf * zalf1 * ( zpsx  - zfx(ji-1,jj) ) - ( zalf1 - zalf ) * ztemp ) ) &
               &            + zbt1 * zpsxx
            zpsxy =  zbt  * ( zalf * zfxy(ji-1,jj) + zalf1 * zpsxy            &
               &            + 3.0 * (- zalf1*zfy(ji-1,jj)  + zalf * zpsy ) )  &
               &            + zbt1 * zpsxy
            zpsy  =  zbt * ( zpsy  + zfy (ji-1,jj) ) + zbt1 * zpsy
            zpsyy =  zbt * ( zpsyy + zfyy(ji-1,jj) ) + zbt1 * zpsyy

            !                                !  Flux from i+1 to i IF u LT 0.
            zbt   =       zbet(ji,jj)
            zbt1  = 1.0 - zbet(ji,jj)
            zpsm  = zbt * zpsm + zbt1 * ( zpsm + zfm(ji,jj) )
            zalf  = zbt1 * zfm(ji,jj) / zpsm
            zalf1 = 1.0 - zalf
            ztemp = - zalf * zps0 + zalf1 * zf0(ji,jj)
            !
            zps0  = zbt * zps0  + zbt1 * ( zps0 + zf0(ji,jj) )
            zpsx  = zbt * zpsx  + zbt1 * ( zalf * zfx(ji,jj) + zalf1 * zpsx + 3.0 * ztemp )
            zpsxx = zbt * zpsxx + zbt1 * ( zalf * zalf * zfxx(ji,jj) + zalf1 * zalf1 * zpsxx &
               &                         + 5.0 * ( zalf * zalf1 * ( - zpsx + zfx(ji,jj) )    &
               &                         + ( zalf1 - zalf ) * ztemp ) )
            zpsxy = zbt * zpsxy + zbt1 * ( zalf * zfxy(ji,jj) + zalf1 * zpsxy  &
               &                         + 3.0 * ( zalf1 * zfy(ji,jj) - zalf * zpsy ) )
            zpsy  = zbt * zpsy  + zbt1 * ( zpsy  + zfy (ji,jj) )
            zpsyy = zbt * zpsyy + zbt1 * ( zpsyy + zfyy(ji,jj) )
            !
            psm (ji,jj) = zpsm  ! optimization
            ps0 (ji,jj) = zps0
            psx (ji,jj) = zpsx
            psxx(ji,jj) = zpsxx
            psy (ji,jj) = zpsy
            psyy(ji,jj) = zpsyy
            psxy(ji,jj) = zpsxy
            !
      END_2D
      ! ji+ ji- jj+ jj-
      !
   END SUBROUTINE adv_x


   SUBROUTINE adv_y( pdt, pvf, pcrh, pe1e2, psm, ps0, psx, psxx, psy, psyy, psxy )
      !!---------------------------------------------------------------------
      !!                **  routine adv_y  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!                variable on y axis
      !!---------------------------------------------------------------------
      REAL(wp)                , INTENT(in   ) ::   pdt                ! time step
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pvf                ! j-direction ice velocity at U-point [m/s]
      REAL(wp)                , INTENT(in   ) ::   pcrh               ! call adv_x then adv_y (=1) or the opposite (=0)
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pe1e2
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psm                ! area
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ps0                ! field to be advected
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psx , psy          ! 1st moments
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   psxx, psyy, psxy   ! 2nd moments
      !!
      INTEGER  ::   ji, jj                               ! dummy loop indices
      INTEGER  ::   ki0                                  ! dummy loop indices
      REAL(wp) ::   zs1max, zslpmax, ztemp               ! temporary scalars
      REAL(wp) ::   zs1new, zalf , zalfq , zbt           !    -         -
      REAL(wp) ::   zs2new, zalf1, zalf1q, zbt1          !    -         -
      REAL(wp) ::   ze1e2, zpsm, zps0
      REAL(wp) ::   zpsx, zpsy, zpsxx, zpsyy, zpsxy
      REAL(wp) ::   zmask, zV
      REAL(wp), DIMENSION(jpi,jpj) ::   zf0, zfx , zfy , zbet   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zfm, zfxx, zfyy, zfxy   !  -      -
      REAL(wp), DIMENSION(jpi,jpj) ::   zalg, zalg1, zalg1q     !  -      -
      !---------------------------------------------------------------------
      ! in order to avoid lbc_lnk (communications):
      !    ji loop must be 1:jpi   if adv_y is called first
      !                and 2:jpi-1 if adv_y is called second
      ki0 = NINT(pcrh)
      !
      ! Limitation of moments.
      DO_2D( ki0, ki0, 1, 0 )
            !! Only `jj+1` needed
            !! => [0:jpj-1] for `jj` !
            !
            zmask = xmskf(ji,jj)
            zV    = pvf(ji,jj+1)
            !
            ze1e2 = pe1e2(ji,jj)
            zpsm  = psm (ji,jj) ! optimization
            zps0  = ps0 (ji,jj)
            zpsx  = psx (ji,jj)
            zpsxx = psxx(ji,jj)
            zpsy  = psy (ji,jj)
            zpsyy = psyy(ji,jj)
            zpsxy = psxy(ji,jj)

            !  Initialize volumes of boxes (=area if adv_y first called, =psm otherwise)
            zpsm = MAX( pcrh * ze1e2  +  ( 1.0 - pcrh ) * zpsm , epsi20 )
            !
            zslpmax = MAX( 0._wp, zps0 )
            zs1max  = 1.5 * zslpmax
            zs1new  = MIN( zs1max, MAX( -zs1max, zpsy ) )
            zs2new  = MIN( 2.0 * zslpmax - 0.3334 * ABS( zs1new ), MAX( ABS( zs1new ) - zslpmax, zpsyy ) )
            rswitch = ( 1.0 - MAX( 0._wp, SIGN( 1._wp, -zslpmax) ) ) * zmask ! Case of empty boxes & Apply mask

            zps0  = zslpmax
            zpsx  = zpsx   * rswitch
            zpsxx = zpsxx  * rswitch
            zpsy  = zs1new * rswitch
            zpsyy = zs2new * rswitch
            zpsxy = MIN( zslpmax, MAX( -zslpmax, zpsxy ) ) * rswitch

            !  Calculate fluxes and moments between boxes j<-->j+1
            !                                !  Flux from j to j+1 WHEN v GT 0
            zbet(ji,jj)  =  MAX( 0._wp, SIGN( 1._wp, zV ) )
            zalf         =  MAX( 0._wp, zV ) * pdt / zpsm
            zalfq        =  zalf * zalf
            zalf1        =  1.0 - zalf
            zalf1q       =  zalf1 * zalf1
            !
            zfm (ji,jj)  =  zalf  * zpsm
            zf0 (ji,jj)  =  zalf  * ( zps0 + zalf1 * ( zpsy  + (zalf1-zalf) * zpsyy ) )
            zfy (ji,jj)  =  zalfq * ( zpsy + 3.0*zalf1*zpsyy )
            zfyy(ji,jj)  =  zalf  * zalfq * zpsyy
            zfx (ji,jj)  =  zalf  * ( zpsx + zalf1 * zpsxy )
            zfxy(ji,jj)  =  zalfq * zpsxy
            zfxx(ji,jj)  =  zalf  * zpsxx
            !                                !  Readjust moments remaining in the box.
            zpsm  =  zpsm  - zfm(ji,jj)
            zps0  =  zps0  - zf0(ji,jj)
            zpsy  =  zalf1q * ( zpsy -3.0 * zalf * zpsyy )
            zpsyy =  zalf1 * zalf1q * zpsyy
            zpsx  =  zpsx  - zfx(ji,jj)
            zpsxx =  zpsxx - zfxx(ji,jj)
            zpsxy =  zalf1q * zpsxy
            !
            psm (ji,jj) = zpsm ! optimization
            ps0 (ji,jj) = zps0
            psx (ji,jj) = zpsx
            psxx(ji,jj) = zpsxx
            psy (ji,jj) = zpsy
            psyy(ji,jj) = zpsyy
            psxy(ji,jj) = zpsxy
            !
      END_2D
      ! ji+ ji- jj+ jj-
      !LB: because of the F-centric world, need to link here => to fill the j+1
      CALL lbc_lnk( 'adv_x', zfm,'F', 1._wp,  zf0,'F',1._wp, zfxy,'F',1._wp, &
         &                   zfx,'F',-1._wp, zfxx,'F',1._wp, zbet,'F',1._wp, &
         &                   zfy,'F',-1._wp, zfyy,'F',1._wp,                 &
         &                   psm,'F', 1._wp,  ps0,'F',1._wp, psxy,'F',1._wp, &
         &                   psx,'F',-1._wp, psxx,'F',1._wp,                 &
         &                   psy,'F',-1._wp, psyy,'F',1._wp  )
      !! Fix due to the use of lbc_lnk, to avoid division by 0 later on...
      IF( ANY(psm<0._wp) ) CALL ctl_stop('STOP', 'adv_y() [icedyn_adv_pra_bbm_f.F90] : `psm<0` somewhere')
      WHERE( psm < epsi20 ) psm = epsi20
      !
      DO_2D( ki0, ki0, 1, 0 )
            !! Only `jj+1` needed
            !! => [0:jpj-1] for `jj` !
            zV = pvf(ji,jj+1)
            !                                !  Flux from j+1 to j when v LT 0.
            zalf          = MAX( 0._wp, -zV ) * pdt / psm(ji,jj+1)
            zalg  (ji,jj) = zalf
            zalfq         = zalf * zalf
            zalf1         = 1.0 - zalf
            zalg1 (ji,jj) = zalf1
            zalf1q        = zalf1 * zalf1
            zalg1q(ji,jj) = zalf1q
            !
            zfm   (ji,jj) = zfm (ji,jj) + zalf  *    psm (ji,jj+1)
            zf0   (ji,jj) = zf0 (ji,jj) + zalf  * (  ps0 (ji,jj+1) &
               &            - zalf1 * (psy(ji,jj+1) - (zalf1 - zalf ) * psyy(ji,jj+1) ) )
            zfy   (ji,jj) = zfy (ji,jj) + zalfq * (  psy (ji,jj+1) - 3.0 * zalf1 * psyy(ji,jj+1) )
            zfyy  (ji,jj) = zfyy(ji,jj) + zalf  *    psyy(ji,jj+1) * zalfq
            zfx   (ji,jj) = zfx (ji,jj) + zalf  * (  psx (ji,jj+1) - zalf1 * psxy(ji,jj+1) )
            zfxy  (ji,jj) = zfxy(ji,jj) + zalfq *    psxy(ji,jj+1)
            zfxx  (ji,jj) = zfxx(ji,jj) + zalf  *    psxx(ji,jj+1)
      END_2D
      ! ji+ ji- jj+ jj-

      DO_2D( ki0, ki0, 0, 1 )
            !! Only `jj-1` needed
            !! => [2:jpj] for `jj` !
            !                                !  Readjust moments remaining in the box.
            zbt  =       zbet(ji,jj-1)
            zbt1 = 1.0 - zbet(ji,jj-1)
            !
            zpsm  = psm (ji,jj) ! optimization
            zps0  = ps0 (ji,jj)
            zpsx  = psx (ji,jj)
            zpsxx = psxx(ji,jj)
            zpsy  = psy (ji,jj)
            zpsyy = psyy(ji,jj)
            zpsxy = psxy(ji,jj)
            !
            zpsm  = zbt * zpsm + zbt1 * ( zpsm - zfm(ji,jj-1) )
            zps0  = zbt * zps0 + zbt1 * ( zps0 - zf0(ji,jj-1) )
            zpsy  = zalg1q(ji,jj-1) * ( zpsy + 3.0 * zalg(ji,jj-1) * zpsyy )
            zpsyy = zalg1 (ji,jj-1) * zalg1q(ji,jj-1) * zpsyy
            zpsx  = zbt * zpsx  + zbt1 * ( zpsx  - zfx (ji,jj-1) )
            zpsxx = zbt * zpsxx + zbt1 * ( zpsxx - zfxx(ji,jj-1) )
            zpsxy = zalg1q(ji,jj-1) * zpsxy

            !   Put the temporary moments into appropriate neighboring boxes.
            !                                !   Flux from j to j+1 IF v GT 0.
            zbt   =       zbet(ji,jj-1)
            zbt1  = 1.0 - zbet(ji,jj-1)
            zpsm  = zbt * ( zpsm + zfm(ji,jj-1) ) + zbt1 * zpsm
            zalf  = zbt * zfm(ji,jj-1) / zpsm
            zalf1 = 1.0 - zalf
            ztemp = zalf * zps0 - zalf1 * zf0(ji,jj-1)
            !
            zps0  =  zbt  * ( zps0 + zf0(ji,jj-1) ) + zbt1 * zps0
            zpsy  =  zbt  * ( zalf * zfy(ji,jj-1) + zalf1 * zpsy + 3.0 * ztemp ) + zbt1 * zpsy
            zpsyy =  zbt  * ( zalf * zalf * zfyy(ji,jj-1) + zalf1 * zalf1 * zpsyy                           &
               &            + 5.0 * ( zalf * zalf1 * ( zpsy - zfy(ji,jj-1) ) - ( zalf1 - zalf ) * ztemp ) ) &
               &            + zbt1 * zpsyy
            zpsxy =  zbt  * ( zalf * zfxy(ji,jj-1) + zalf1 * zpsxy             &
               &            + 3.0 * (- zalf1 * zfx(ji,jj-1) + zalf * zpsx ) )  &
               &            + zbt1 * zpsxy
            zpsx  =  zbt * ( zpsx  + zfx (ji,jj-1) ) + zbt1 * zpsx
            zpsxx =  zbt * ( zpsxx + zfxx(ji,jj-1) ) + zbt1 * zpsxx

            !                                !  Flux from j+1 to j IF v LT 0.
            zbt   =       zbet(ji,jj)
            zbt1  = 1.0 - zbet(ji,jj)
            zpsm  = zbt * zpsm + zbt1 * ( zpsm + zfm(ji,jj) )
            zalf  = zbt1 * zfm(ji,jj) / zpsm
            zalf1 = 1.0 - zalf
            ztemp = - zalf * zps0 + zalf1 * zf0(ji,jj)
            !
            zps0  = zbt * zps0  + zbt1 * (  zps0 + zf0(ji,jj) )
            zpsy  = zbt * zpsy  + zbt1 * (  zalf * zfy(ji,jj) + zalf1 * zpsy + 3.0 * ztemp )
            zpsyy = zbt * zpsyy + zbt1 * (  zalf * zalf * zfyy(ji,jj) + zalf1 * zalf1 * zpsyy &
               &                         + 5.0 * ( zalf * zalf1 * ( - zpsy + zfy(ji,jj) )     &
               &                         + ( zalf1 - zalf ) * ztemp ) )
            zpsxy = zbt * zpsxy + zbt1 * (  zalf * zfxy(ji,jj) + zalf1 * zpsxy  &
               &                         + 3.0 * ( zalf1 * zfx(ji,jj) - zalf * zpsx ) )
            zpsx  = zbt * zpsx  + zbt1 * ( zpsx  + zfx (ji,jj) )
            zpsxx = zbt * zpsxx + zbt1 * ( zpsxx + zfxx(ji,jj) )
            !
            psm (ji,jj) = zpsm ! optimization
            ps0 (ji,jj) = zps0
            psx (ji,jj) = zpsx
            psxx(ji,jj) = zpsxx
            psy (ji,jj) = zpsy
            psyy(ji,jj) = zpsyy
            psxy(ji,jj) = zpsxy
            !
      END_2D
      ! ji+ ji- jj+ jj-
      !
   END SUBROUTINE adv_y


   SUBROUTINE adv_pra_f_d_init( )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE adv_pra_f_d_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for Prather advection
      !!-------------------------------------------------------------------
      INTEGER ::   ierr
      CHARACTER(len=16), PARAMETER :: crtnm = 'adv_pra_f_d_init'
      !!-------------------------------------------------------------------

      !                             !* allocate prather fields
      ALLOCATE( sx1md(jpi,jpj), sy1md(jpi,jpj), sxx1md(jpi,jpj), syy1md(jpi,jpj), sxy1md(jpi,jpj) , &
         &      sxdd1(jpi,jpj), sydd1(jpi,jpj), sxxdd1(jpi,jpj), syydd1(jpi,jpj), sxydd1(jpi,jpj) , &
         &      sxdd2(jpi,jpj), sydd2(jpi,jpj), sxxdd2(jpi,jpj), syydd2(jpi,jpj), sxydd2(jpi,jpj) , &
         &      sxdd3(jpi,jpj), sydd3(jpi,jpj), sxxdd3(jpi,jpj), syydd3(jpi,jpj), sxydd3(jpi,jpj) , &
         &      STAT = ierr )
      CALL mpp_sum( crtnm, ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', crtnm//' : unable to allocate 1md array for Prather advection scheme')
      !
      CALL adv_pra_f_d_rst( 'READ' )    !* read or initialize all required files
      !
   END SUBROUTINE adv_pra_f_d_init


   SUBROUTINE adv_pra_f_d_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE adv_pra_f_d_rst  ***
      !!
      !! ** Purpose :   Read or write file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER ::   iter     ! local integer
      INTEGER ::   id1      ! local integer
      CHARACTER(len=1), PARAMETER :: cg = 'f'
      !!----------------------------------------------------------------------
      !
      !                                      !==========================!
      IF( TRIM(cdrw) == 'READ' ) THEN        !==  Read or initialize  ==!
         !                                   !==========================!
         !
         IF( ln_rstart ) THEN
            id1 = iom_varid( numrir, 'sx1md'//cg , ldstop = .FALSE. )    ! file exist: id1>0
         ELSE
            id1 = 0                                                  ! no restart: id1=0
         ENDIF
         !
         IF( id1 > 0 ) THEN                     !**  Read the restart file  **!
            !
            !                                                        ! ice damage !#bbm
            CALL iom_get( numrir, jpdom_auto, 'sx1md'//cg , sx1md,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sy1md'//cg , sy1md,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxx1md'//cg, sxx1md )
            CALL iom_get( numrir, jpdom_auto, 'syy1md'//cg, syy1md )
            CALL iom_get( numrir, jpdom_auto, 'sxy1md'//cg, sxy1md )
            !
            CALL iom_get( numrir, jpdom_auto, 'sxdd1'//cg , sxdd1,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sydd1'//cg , sydd1,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxdd1'//cg, sxxdd1 )
            CALL iom_get( numrir, jpdom_auto, 'syydd1'//cg, syydd1 )
            CALL iom_get( numrir, jpdom_auto, 'sxydd1'//cg, sxydd1 )
            !
            CALL iom_get( numrir, jpdom_auto, 'sxdd2'//cg , sxdd2,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sydd2'//cg , sydd2,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxdd2'//cg, sxxdd2 )
            CALL iom_get( numrir, jpdom_auto, 'syydd2'//cg, syydd2 )
            CALL iom_get( numrir, jpdom_auto, 'sxydd2'//cg, sxydd2 )
            !
            CALL iom_get( numrir, jpdom_auto, 'sxdd3'//cg , sxdd3,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sydd3'//cg , sydd3,  psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxdd3'//cg, sxxdd3 )
            CALL iom_get( numrir, jpdom_auto, 'syydd3'//cg, syydd3 )
            CALL iom_get( numrir, jpdom_auto, 'sxydd3'//cg, sxydd3 )
            !
            !
         ELSE                                   !**  start rheology from rest  **!
            !
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest OR previous run without Prather, set moments to 0'
            !
            sx1md = 0._wp   ;   sy1md = 0._wp   ;   sxx1md = 0._wp   ;   syy1md = 0._wp   ;   sxy1md = 0._wp      ! ice damage !#bbm
            sxdd1 = 0._wp   ;   sydd1 = 0._wp   ;   sxxdd1 = 0._wp   ;   syydd1 = 0._wp   ;   sxydd1 = 0._wp      !  !#bbm
            sxdd2 = 0._wp   ;   sydd2 = 0._wp   ;   sxxdd2 = 0._wp   ;   syydd2 = 0._wp   ;   sxydd2 = 0._wp      !  !#bbm
            sxdd3 = 0._wp   ;   sydd3 = 0._wp   ;   sxxdd3 = 0._wp   ;   syydd3 = 0._wp   ;   sxydd3 = 0._wp      !  !#bbm
            !
         ENDIF !IF( id1 > 0 )
         !                                   !=====================================!
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   !==  write in the ice restart file  ==!
         !                                   !=====================================!
         IF(lwp) WRITE(numout,*) '----  ice-adv-rst  ----'
         iter = kt + nn_fsbc - 1             ! ice restarts are written at kt == nitrst - nn_fsbc + 1
         !
         !
         ! In case Prather scheme is used for advection, write second order moments
         ! ------------------------------------------------------------------------
         !
         !                                                           ! ice thickness
         CALL iom_rstput( iter, nitrst, numriw, 'sx1md'//cg , sx1md  )
         CALL iom_rstput( iter, nitrst, numriw, 'sy1md'//cg , sy1md  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxx1md'//cg, sxx1md )
         CALL iom_rstput( iter, nitrst, numriw, 'syy1md'//cg, syy1md )
         CALL iom_rstput( iter, nitrst, numriw, 'sxy1md'//cg, sxy1md )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sxdd1'//cg , sxdd1  )
         CALL iom_rstput( iter, nitrst, numriw, 'sydd1'//cg , sydd1  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxdd1'//cg, sxxdd1 )
         CALL iom_rstput( iter, nitrst, numriw, 'syydd1'//cg, syydd1 )
         CALL iom_rstput( iter, nitrst, numriw, 'sxydd1'//cg, sxydd1 )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sxdd2'//cg , sxdd2  )
         CALL iom_rstput( iter, nitrst, numriw, 'sydd2'//cg , sydd2  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxdd2'//cg, sxxdd2 )
         CALL iom_rstput( iter, nitrst, numriw, 'syydd2'//cg, syydd2 )
         CALL iom_rstput( iter, nitrst, numriw, 'sxydd2'//cg, sxydd2 )
         !
         CALL iom_rstput( iter, nitrst, numriw, 'sxdd3'//cg , sxdd3  )
         CALL iom_rstput( iter, nitrst, numriw, 'sydd3'//cg , sydd3  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxdd3'//cg, sxxdd3 )
         CALL iom_rstput( iter, nitrst, numriw, 'syydd3'//cg, syydd3 )
         CALL iom_rstput( iter, nitrst, numriw, 'sxydd3'//cg, sxydd3 )
         !
      ENDIF
      !
   END SUBROUTINE adv_pra_f_d_rst


#else
   !!----------------------------------------------------------------------
   !!   Default option            Dummy module        NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv_pra_bbm_f
