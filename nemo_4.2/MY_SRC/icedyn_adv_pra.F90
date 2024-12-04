MODULE icedyn_adv_pra
   !!======================================================================
   !!                       ***  MODULE icedyn_adv_pra   ***
   !!   sea-ice : advection => Prather scheme
   !!======================================================================
   !! History :       !  2008-03  (M. Vancoppenolle) original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!--------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv_pra : advection of sea ice using Prather scheme
   !!   adv_x, adv_y    : Prather scheme applied in i- and j-direction, resp.
   !!   adv_pra_init    : initialisation of the Prather scheme
   !!   adv_pra_rst     : read/write Prather field in ice restart file, or initialized to zero
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

   PUBLIC   ice_dyn_adv_pra   ! called by icedyn_adv
   PUBLIC   adv_pra_init      ! called by icedyn_adv

   ! Moments for advection
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxice, syice, sxxice, syyice, sxyice   ! ice thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxsn , sysn , sxxsn , syysn , sxysn    ! snow thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxa  , sya  , sxxa  , syya  , sxya     ! ice concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxsal, sysal, sxxsal, syysal, sxysal   ! ice salinity
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxage, syage, sxxage, syyage, sxyage   ! ice age
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxc0 , syc0 , sxxc0 , syyc0 , sxyc0    ! snow layers heat content
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   sxe  , sye  , sxxe  , syye  , sxye     ! ice layers heat content
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxap , syap , sxxap , syyap , sxyap    ! melt pond fraction
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxvp , syvp , sxxvp , syyvp , sxyvp    ! melt pond volume
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   sxvl , syvl , sxxvl , syyvl , sxyvl    ! melt pond lid volume

   REAL(wp), PARAMETER :: rr_scl_fct = 1.E-6_wp
   REAL(wp), PARAMETER :: r1_scl_fct = 1.E6_wp

   
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icedyn_adv_pra.F90 15049 2021-06-23 16:17:30Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv_pra(         kt, pu_ice, pv_ice, ph_i, ph_s, ph_ip,  &
      &                        pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i )
      !!----------------------------------------------------------------------
      !!                **  routine ice_dyn_adv_pra  **
      !!
      !! ** purpose :   Computes and adds the advection trend to sea-ice
      !!
      !! ** method  :   Uses Prather second order scheme that advects tracers
      !!                but also their quadratic forms. The method preserves
      !!                tracer structures by conserving second order moments.
      !!
      !! Reference:  Prather, 1986, JGR, 91, D6. 6671-6681.
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kt         ! time step
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pu_ice     ! ice i-velocity
      REAL(wp), DIMENSION(:,:)    , INTENT(in   ) ::   pv_ice     ! ice j-velocity
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   ph_i       ! ice thickness
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   ph_s       ! snw thickness
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   ph_ip      ! ice pond thickness
      REAL(wp), DIMENSION(:,:)    , INTENT(inout) ::   pato_i     ! open water area
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_il      ! melt pond lid thickness
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_i       ! ice heat content
      !
      INTEGER  ::   ji, jj, jk, jl, jt      ! dummy loop indices
      INTEGER  ::   icycle                  ! number of sub-timestep for the advection
      REAL(wp) ::   zdt, z1_dt              !   -      -
      REAL(wp), DIMENSION(1)                  ::   zcflprv, zcflnow   ! for global communication
      REAL(wp), DIMENSION(jpi,jpj)            ::   zati1, zati2
      REAL(wp), DIMENSION(jpi,jpj)            ::   zudy, zvdx
      REAL(wp), DIMENSION(jpi,jpj)            ::   z2d
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   zhi_max, zhs_max, zhip_max, zs_i, zsi_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl) ::   ze_i, zei_max
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl) ::   ze_s, zes_max
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   zarea
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   z0ice, z0snw, z0ai, z0smi, z0oi
      REAL(wp), DIMENSION(jpi,jpj,jpl)        ::   z0ap , z0vp, z0vl
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl) ::   z0es
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl) ::   z0ei
      !! diagnostics
      REAL(wp), DIMENSION(jpi,jpj)            ::   zdiag_adv_mass, zdiag_adv_salt, zdiag_adv_heat
      CHARACTER(len=19), PARAMETER :: crtnm = 'ice_dyn_adv_pra'
      CHARACTER(len=1),  PARAMETER :: cgrt    = 'T'
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- '//crtnm//': Prather advection scheme'
      !
      ! --- Record max of the surrounding 9-pts (for call Hbig) --- !
      ! thickness and salinity
      WHERE( pv_i(:,:,:) >= epsi10 )
         zs_i(:,:,:) = psv_i(:,:,:) / pv_i(:,:,:)
      ELSEWHERE
         zs_i(:,:,:) = 0._wp
      END WHERE
      CALL icemax3D( ph_i , zhi_max )
      CALL icemax3D( ph_s , zhs_max )
      CALL icemax3D( ph_ip, zhip_max)
      CALL icemax3D( zs_i , zsi_max )
      CALL lbc_lnk( crtnm, zhi_max, cgrt, 1._wp, zhs_max, cgrt, 1._wp, zhip_max, cgrt, 1._wp, zsi_max, cgrt, 1._wp )
      !
      ! enthalpies
      DO jk = 1, nlay_i
         WHERE( pv_i(:,:,:) >= epsi10 )
            ze_i(:,:,jk,:) = pe_i(:,:,jk,:) / pv_i(:,:,:)
         ELSEWHERE
            ze_i(:,:,jk,:) = 0._wp
         END WHERE
      END DO
      DO jk = 1, nlay_s
         WHERE( pv_s(:,:,:) >= epsi10 )
            ze_s(:,:,jk,:) = pe_s(:,:,jk,:) / pv_s(:,:,:)
         ELSEWHERE
            ze_s(:,:,jk,:) = 0._wp
         END WHERE
      END DO
      CALL icemax4D( ze_i , zei_max )
      CALL icemax4D( ze_s , zes_max )
      CALL lbc_lnk( crtnm, zei_max, cgrt, 1._wp )
      CALL lbc_lnk( crtnm, zes_max, cgrt, 1._wp )
      !
      !
      ! --- If ice drift is too fast, use  subtime steps for advection (CFL test for stability) --- !
      !        Note: the advection split is applied at the next time-step in order to avoid blocking global comm.
      !              this should not affect too much the stability
      zcflnow(1) =                  MAXVAL( ABS( pu_ice(:,:) ) * rDt_ice * r1_e1u(:,:) )
      zcflnow(1) = MAX( zcflnow(1), MAXVAL( ABS( pv_ice(:,:) ) * rDt_ice * r1_e2v(:,:) ) )

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
      zudy(:,:) = pu_ice(:,:) * e2u(:,:) * rr_scl_fct    ! `rr_scl_fct` because we work with km^2 rather than m^2 !!!
      zvdx(:,:) = pv_ice(:,:) * e1v(:,:) * rr_scl_fct    !       "                     "

      DO jt = 1, icycle

         ! diagnostics
         zdiag_adv_mass(:,:) =   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
            &                  + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow
         zdiag_adv_salt(:,:) =   SUM( psv_i(:,:,:) , dim=3 ) * rhoi
         zdiag_adv_heat(:,:) = - SUM(SUM( pe_i(:,:,1:nlay_i,:) , dim=4 ), dim=3 ) &
            &                  - SUM(SUM( pe_s(:,:,1:nlay_s,:) , dim=4 ), dim=3 )

         ! record at_i before advection (for open water)
         zati1(:,:) = SUM( pa_i(:,:,:), dim=3 )

         z2d(:,:) = e1e2t(:,:) * rr_scl_fct * xmskt(:,:)  ! cell area in km^2

         ! --- transported fields --- !
         DO jl = 1, jpl
            zarea(:,:,jl) = z2d(:,:)
            z0snw(:,:,jl) = pv_s (:,:,jl) * z2d(:,:)        ! Snow volume
            z0ice(:,:,jl) = pv_i (:,:,jl) * z2d(:,:)        ! Ice  volume
            z0ai (:,:,jl) = pa_i (:,:,jl) * z2d(:,:)        ! Ice area
            z0smi(:,:,jl) = psv_i(:,:,jl) * z2d(:,:)        ! Salt content
            z0oi (:,:,jl) = poa_i(:,:,jl) * z2d(:,:)        ! Age content
            DO jk = 1, nlay_s
               z0es(:,:,jk,jl) = pe_s(:,:,jk,jl) * z2d(:,:) ! Snow heat content
            END DO
            DO jk = 1, nlay_i
               z0ei(:,:,jk,jl) = pe_i(:,:,jk,jl) * z2d(:,:) ! Ice  heat content
            END DO
            IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               z0ap(:,:,jl) = pa_ip(:,:,jl) * z2d(:,:)      ! Melt pond fraction
               z0vp(:,:,jl) = pv_ip(:,:,jl) * z2d(:,:)      ! Melt pond volume
               IF ( ln_pnd_lids ) THEN
                  z0vl(:,:,jl) = pv_il(:,:,jl) * z2d(:,:)   ! Melt pond lid volume
               ENDIF
            ENDIF
         END DO

         IF( MOD( (kt - 1) / nn_fsbc , 2 ) ==  MOD( (jt - 1) , 2 ) ) THEN   !==  odd ice time step:  adv_x then adv_y  ==!
            !                                                               !--------------------------------------------!
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0ice, sxice, sxxice, syice, syyice, sxyice ) !--- ice volume
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0ice, sxice, sxxice, syice, syyice, sxyice )
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0snw, sxsn , sxxsn , sysn , syysn , sxysn  ) !--- snow volume
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0snw, sxsn , sxxsn , sysn , syysn , sxysn  )
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0smi, sxsal, sxxsal, sysal, syysal, sxysal ) !--- ice salinity
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0smi, sxsal, sxxsal, sysal, syysal, sxysal )
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0ai , sxa  , sxxa  , sya  , syya  , sxya   ) !--- ice concentration
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0ai , sxa  , sxxa  , sya  , syya  , sxya   )
            CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0oi , sxage, sxxage, syage, syyage, sxyage ) !--- ice age
            CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0oi , sxage, sxxage, syage, syyage, sxyage )
            !
            DO jk = 1, nlay_s                                                                           !--- snow heat content
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0es (:,:,jk,:), sxc0(:,:,jk,:),   &
                  &                                 sxxc0(:,:,jk,:), syc0(:,:,jk,:), syyc0(:,:,jk,:), sxyc0(:,:,jk,:) )
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0es (:,:,jk,:), sxc0(:,:,jk,:),   &
                  &                                 sxxc0(:,:,jk,:), syc0(:,:,jk,:), syyc0(:,:,jk,:), sxyc0(:,:,jk,:) )
            END DO
            DO jk = 1, nlay_i                                                                           !--- ice heat content
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0ei(:,:,jk,:), sxe(:,:,jk,:),   &
                  &                                 sxxe(:,:,jk,:), sye(:,:,jk,:), syye(:,:,jk,:), sxye(:,:,jk,:) )
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0ei(:,:,jk,:), sxe(:,:,jk,:),   &
                  &                                 sxxe(:,:,jk,:), sye(:,:,jk,:), syye(:,:,jk,:), sxye(:,:,jk,:) )
            END DO
            !
            IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0ap, sxap, sxxap, syap, syyap, sxyap )    !--- melt pond fraction
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0ap, sxap, sxxap, syap, syyap, sxyap )
               CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )    !--- melt pond volume
               CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )
               IF ( ln_pnd_lids ) THEN
                  CALL adv_x( zdt, zudy, 1._wp, z2d, zarea, z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl ) !--- melt pond lid volume
                  CALL adv_y( zdt, zvdx, 0._wp, z2d, zarea, z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl )
               ENDIF
            ENDIF
            !                                                               !--------------------------------------------!
         ELSE                                                               !== even ice time step:  adv_y then adv_x  ==!
            !                                                               !--------------------------------------------!
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0ice, sxice, sxxice, syice, syyice, sxyice ) !--- ice volume
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0ice, sxice, sxxice, syice, syyice, sxyice )
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0snw, sxsn , sxxsn , sysn , syysn , sxysn  ) !--- snow volume
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0snw, sxsn , sxxsn , sysn , syysn , sxysn  )
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0smi, sxsal, sxxsal, sysal, syysal, sxysal ) !--- ice salinity
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0smi, sxsal, sxxsal, sysal, syysal, sxysal )
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0ai , sxa  , sxxa  , sya  , syya  , sxya   ) !--- ice concentration
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0ai , sxa  , sxxa  , sya  , syya  , sxya   )
            CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0oi , sxage, sxxage, syage, syyage, sxyage ) !--- ice age
            CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0oi , sxage, sxxage, syage, syyage, sxyage )
            DO jk = 1, nlay_s                                                                           !--- snow heat content
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0es (:,:,jk,:), sxc0(:,:,jk,:),   &
                  &                                 sxxc0(:,:,jk,:), syc0(:,:,jk,:), syyc0(:,:,jk,:), sxyc0(:,:,jk,:) )
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0es (:,:,jk,:), sxc0(:,:,jk,:),   &
                  &                                 sxxc0(:,:,jk,:), syc0(:,:,jk,:), syyc0(:,:,jk,:), sxyc0(:,:,jk,:) )
            END DO
            DO jk = 1, nlay_i                                                                           !--- ice heat content
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0ei(:,:,jk,:), sxe(:,:,jk,:),   &
                  &                                 sxxe(:,:,jk,:), sye(:,:,jk,:), syye(:,:,jk,:), sxye(:,:,jk,:) )
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0ei(:,:,jk,:), sxe(:,:,jk,:),   &
                  &                                 sxxe(:,:,jk,:), sye(:,:,jk,:), syye(:,:,jk,:), sxye(:,:,jk,:) )
            END DO
            IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0ap, sxap, sxxap, syap, syyap, sxyap )    !--- melt pond fraction
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0ap, sxap, sxxap, syap, syyap, sxyap )
               CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )    !--- melt pond volume
               CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0vp, sxvp, sxxvp, syvp, syyvp, sxyvp )
               IF ( ln_pnd_lids ) THEN
                  CALL adv_y( zdt, zvdx, 1._wp, z2d, zarea, z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl ) !--- melt pond lid volume
                  CALL adv_x( zdt, zudy, 0._wp, z2d, zarea, z0vl, sxvl, sxxvl, syvl, syyvl, sxyvl )
               ENDIF
            ENDIF
            !
         ENDIF

         ! --- Lateral boundary conditions --- !
         !     caution: for gradients (sx and sy) the sign changes
         CALL lbc_lnk( crtnm, z0ice , cgrt, 1._wp, sxice , cgrt, -1._wp, syice , cgrt, -1._wp  & ! ice volume
            &                          , sxxice, cgrt, 1._wp, syyice, cgrt,  1._wp, sxyice, cgrt,  1._wp  &
            &                          , z0snw , cgrt, 1._wp, sxsn  , cgrt, -1._wp, sysn  , cgrt, -1._wp  & ! snw volume
            &                          , sxxsn , cgrt, 1._wp, syysn , cgrt,  1._wp, sxysn , cgrt,  1._wp  &
            &                          , z0smi , cgrt, 1._wp, sxsal , cgrt, -1._wp, sysal , cgrt, -1._wp  & ! ice salinity
            &                          , sxxsal, cgrt, 1._wp, syysal, cgrt,  1._wp, sxysal, cgrt,  1._wp  &
            &                          , z0ai  , cgrt, 1._wp, sxa   , cgrt, -1._wp, sya   , cgrt, -1._wp  & ! ice concentration
            &                          , sxxa  , cgrt, 1._wp, syya  , cgrt,  1._wp, sxya  , cgrt,  1._wp  &
            &                          , z0oi  , cgrt, 1._wp, sxage , cgrt, -1._wp, syage , cgrt, -1._wp  & ! ice age
            &                          , sxxage, cgrt, 1._wp, syyage, cgrt,  1._wp, sxyage, cgrt,  1._wp  )
         CALL lbc_lnk( crtnm, z0es  , cgrt, 1._wp, sxc0  , cgrt, -1._wp, syc0  , cgrt, -1._wp  & ! snw enthalpy
            &                          , sxxc0 , cgrt, 1._wp, syyc0 , cgrt,  1._wp, sxyc0 , cgrt,  1._wp  )
         CALL lbc_lnk( crtnm, z0ei  , cgrt, 1._wp, sxe   , cgrt, -1._wp, sye   , cgrt, -1._wp  & ! ice enthalpy
            &                          , sxxe  , cgrt, 1._wp, syye  , cgrt,  1._wp, sxye  , cgrt,  1._wp  )
         IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
            IF( ln_pnd_lids ) THEN
               CALL lbc_lnk( crtnm, z0ap , cgrt, 1._wp, sxap , cgrt, -1._wp, syap , cgrt, -1._wp  & ! melt pond fraction
                  &                          , sxxap, cgrt, 1._wp, syyap, cgrt,  1._wp, sxyap, cgrt,  1._wp  &
                  &                          , z0vp , cgrt, 1._wp, sxvp , cgrt, -1._wp, syvp , cgrt, -1._wp  & ! melt pond volume
                  &                          , sxxvp, cgrt, 1._wp, syyvp, cgrt,  1._wp, sxyvp, cgrt,  1._wp  &
                  &                          , z0vl , cgrt, 1._wp, sxvl , cgrt, -1._wp, syvl , cgrt, -1._wp  & ! melt pond lid volume
                  &                          , sxxvl, cgrt, 1._wp, syyvl, cgrt,  1._wp, sxyvl, cgrt,  1._wp  )
            ELSE
               CALL lbc_lnk( crtnm, z0ap , cgrt, 1._wp, sxap , cgrt, -1._wp, syap , cgrt, -1._wp  & ! melt pond fraction
                  &                          , sxxap, cgrt, 1._wp, syyap, cgrt,  1._wp, sxyap, cgrt,  1._wp  &
                  &                          , z0vp , cgrt, 1._wp, sxvp , cgrt, -1._wp, syvp , cgrt, -1._wp  & ! melt pond volume
                  &                          , sxxvp, cgrt, 1._wp, syyvp, cgrt,  1._wp, sxyvp, cgrt,  1._wp  )
            ENDIF
         ENDIF

         ! --- Recover the properties from their contents --- !
         z2d(:,:) = r1_scl_fct * r1_e1e2t(:,:) * tmask(:,:,1)
         DO jl = 1, jpl
            pv_i (:,:,jl) = z0ice(:,:,jl) * z2d(:,:)
            pv_s (:,:,jl) = z0snw(:,:,jl) * z2d(:,:)
            psv_i(:,:,jl) = z0smi(:,:,jl) * z2d(:,:)
            poa_i(:,:,jl) = z0oi (:,:,jl) * z2d(:,:)
            pa_i (:,:,jl) = z0ai (:,:,jl) * z2d(:,:)
            DO jk = 1, nlay_s
               pe_s(:,:,jk,jl) = z0es(:,:,jk,jl) * z2d(:,:)
            END DO
            DO jk = 1, nlay_i
               pe_i(:,:,jk,jl) = z0ei(:,:,jk,jl) * z2d(:,:)
            END DO
            IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               pa_ip(:,:,jl) = z0ap(:,:,jl) * z2d(:,:)
               pv_ip(:,:,jl) = z0vp(:,:,jl) * z2d(:,:)
               IF ( ln_pnd_lids ) THEN
                  pv_il(:,:,jl) = z0vl(:,:,jl) * z2d(:,:)
               ENDIF
            ENDIF
         END DO
         !
         ! derive open water from ice concentration
         zudy(:,:) = zudy(:,:) * r1_scl_fct
         zvdx(:,:) = zvdx(:,:) * r1_scl_fct
         
         zati2(:,:) = SUM( pa_i(:,:,:), dim=3 )
         DO_2D( 0, 0, 0, 0 )
               pato_i(ji,jj) = pato_i(ji,jj) - ( zati2(ji,jj) - zati1(ji,jj) ) &                        !--- open water
                  &                          - ( zudy(ji,jj) - zudy(ji-1,jj) + zvdx(ji,jj) - zvdx(ji,jj-1) ) * r1_e1e2t(ji,jj) * zdt
         END_2D
         CALL lbc_lnk( crtnm, pato_i, cgrt,  1.0_wp )
         !
         ! --- diagnostics --- !
         diag_adv_mass(:,:) = diag_adv_mass(:,:) + (   SUM( pv_i (:,:,:) , dim=3 ) * rhoi + SUM( pv_s (:,:,:) , dim=3 ) * rhos &
            &                                        + SUM( pv_ip(:,:,:) , dim=3 ) * rhow + SUM( pv_il(:,:,:) , dim=3 ) * rhow &
            &                                        - zdiag_adv_mass(:,:) ) * z1_dt
         diag_adv_salt(:,:) = diag_adv_salt(:,:) + (   SUM( psv_i(:,:,:) , dim=3 ) * rhoi &
            &                                        - zdiag_adv_salt(:,:) ) * z1_dt
         diag_adv_heat(:,:) = diag_adv_heat(:,:) + ( - SUM(SUM( pe_i(:,:,1:nlay_i,:) , dim=4 ), dim=3 ) &
            &                                        - SUM(SUM( pe_s(:,:,1:nlay_s,:) , dim=4 ), dim=3 ) &
            &                                        - zdiag_adv_heat(:,:) ) * z1_dt
         !
         ! --- Ensure non-negative fields --- !
         !     Remove negative values (conservation is ensured)
         !     (because advected fields are not perfectly bounded and tiny negative values can occur, e.g. -1.e-20)
         CALL ice_var_zapneg( zdt, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i )
         !
         ! --- Make sure ice thickness is not too big --- !
         !     (because ice thickness can be too large where ice concentration is very small)
         CALL Hbig( zdt, zhi_max, zhs_max, zhip_max, zsi_max, zes_max, zei_max, &
            &            pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i, pe_s, pe_i )
         !
         ! --- Ensure snow load is not too big --- !
         CALL Hsnow( zdt, pv_i, pv_s, pa_i, pa_ip, pe_s )
         !
      END DO !DO jt = 1, icycle
      !
      IF( lrst_ice )   CALL adv_pra_rst( 'WRITE', kt )   !* write Prather fields in the restart file
      !
   END SUBROUTINE ice_dyn_adv_pra


#  include "icedyn_adv_pra_adv_x_3d.h90"

#  include "icedyn_adv_pra_adv_y_3d.h90"



   SUBROUTINE Hbig( pdt, phi_max, phs_max, phip_max, psi_max, pes_max, pei_max, &
      &                  pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i, pe_s, pe_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hbig  ***
      !!
      !! ** Purpose : Thickness correction in case advection scheme creates
      !!              abnormally tick ice or snow
      !!
      !! ** Method  : 1- check whether ice thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by adapting ice concentration
      !!              2- check whether snow thickness is larger than the surrounding 9-points
      !!                 (before advection) and reduce it by sending the excess in the ocean
      !!
      !! ** input   : Max thickness of the surrounding 9-points
      !!-------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdt                                   ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)  , INTENT(in   ) ::   phi_max, phs_max, phip_max, psi_max   ! max ice thick from surrounding 9-pts
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pes_max
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in   ) ::   pei_max
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_i, pv_s, pa_i, pa_ip, pv_ip, psv_i
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_s
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_i
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   z1_dt, zhip, zhi, zhs, zsi, zes, zei, zfra
      !!-------------------------------------------------------------------
      !
      z1_dt = 1._wp / pdt
      !
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         IF ( pv_i(ji,jj,jl) > 0._wp ) THEN
            !
            !                               ! -- check h_ip -- !
            ! if h_ip is larger than the surrounding 9 pts => reduce h_ip and increase a_ip
            IF( ln_pnd_LEV .OR. ln_pnd_TOPO .AND. pv_ip(ji,jj,jl) > 0._wp ) THEN
               zhip = pv_ip(ji,jj,jl) / MAX( epsi20, pa_ip(ji,jj,jl) )
               IF( zhip > phip_max(ji,jj,jl) .AND. pa_ip(ji,jj,jl) < 0.15 ) THEN
                  pa_ip(ji,jj,jl) = pv_ip(ji,jj,jl) / phip_max(ji,jj,jl)
               ENDIF
            ENDIF
            !
            !                               ! -- check h_i -- !
            ! if h_i is larger than the surrounding 9 pts => reduce h_i and increase a_i
            zhi = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)
            IF( zhi > phi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
               pa_i(ji,jj,jl) = pv_i(ji,jj,jl) / MIN( phi_max(ji,jj,jl), hi_max(jpl) )   !-- bound h_i to hi_max (99 m)
            ENDIF
            !
            !                               ! -- check h_s -- !
            ! if h_s is larger than the surrounding 9 pts => put the snow excess in the ocean
            zhs = pv_s(ji,jj,jl) / pa_i(ji,jj,jl)
            IF( pv_s(ji,jj,jl) > 0._wp .AND. zhs > phs_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
               zfra = phs_max(ji,jj,jl) / MAX( zhs, epsi20 )
               !
               wfx_res(ji,jj) = wfx_res(ji,jj) + ( pv_s(ji,jj,jl) - pa_i(ji,jj,jl) * phs_max(ji,jj,jl) ) * rhos * z1_dt
               hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
               !
               pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
               pv_s(ji,jj,jl)          = pa_i(ji,jj,jl) * phs_max(ji,jj,jl)
            ENDIF
            !
            !                               ! -- check s_i -- !
            ! if s_i is larger than the surrounding 9 pts => put salt excess in the ocean
            zsi = psv_i(ji,jj,jl) / pv_i(ji,jj,jl)
            IF( zsi > psi_max(ji,jj,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
               zfra = psi_max(ji,jj,jl) / zsi
               sfx_res(ji,jj) = sfx_res(ji,jj) + psv_i(ji,jj,jl) * ( 1._wp - zfra ) * rhoi * z1_dt
               psv_i(ji,jj,jl) = psv_i(ji,jj,jl) * zfra
            ENDIF
            !
         ENDIF
         END_2D
      END DO
      !
      !                                           ! -- check e_i/v_i -- !
      DO jl = 1, jpl
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_i )
         IF ( pv_i(ji,jj,jl) > 0._wp ) THEN
            ! if e_i/v_i is larger than the surrounding 9 pts => put the heat excess in the ocean
            zei = pe_i(ji,jj,jk,jl) / pv_i(ji,jj,jl)
            IF( zei > pei_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
               zfra = pei_max(ji,jj,jk,jl) / zei
               hfx_res(ji,jj) = hfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
               pe_i(ji,jj,jk,jl) = pe_i(ji,jj,jk,jl) * zfra
            ENDIF
         ENDIF
         END_3D
      END DO
      !                                           ! -- check e_s/v_s -- !
      DO jl = 1, jpl
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_s )
         IF ( pv_s(ji,jj,jl) > 0._wp ) THEN
            ! if e_s/v_s is larger than the surrounding 9 pts => put the heat excess in the ocean
            zes = pe_s(ji,jj,jk,jl) / pv_s(ji,jj,jl)
            IF( zes > pes_max(ji,jj,jk,jl) .AND. pa_i(ji,jj,jl) < 0.15 ) THEN
               zfra = pes_max(ji,jj,jk,jl) / zes
               hfx_res(ji,jj) = hfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
               pe_s(ji,jj,jk,jl) = pe_s(ji,jj,jk,jl) * zfra
            ENDIF
         ENDIF
         END_3D
      END DO
      !
   END SUBROUTINE Hbig


   SUBROUTINE Hsnow( pdt, pv_i, pv_s, pa_i, pa_ip, pe_s )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE Hsnow  ***
      !!
      !! ** Purpose : 1- Check snow load after advection
      !!              2- Correct pond concentration to avoid a_ip > a_i
      !!
      !! ** Method :  If snow load makes snow-ice interface to deplet below the ocean surface
      !!              then put the snow excess in the ocean
      !!
      !! ** Notes :   This correction is crucial because of the call to routine icecor afterwards
      !!              which imposes a mini of ice thick. (rn_himin). This imposed mini can artificially
      !!              make the snow very thick (if concentration decreases drastically)
      !!              This behavior has been seen in Ultimate-Macho and supposedly it can also be true for Prather
      !!-------------------------------------------------------------------
      REAL(wp)                    , INTENT(in   ) ::   pdt   ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_i, pv_s, pa_i, pa_ip
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_s
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      REAL(wp) ::   z1_dt, zvs_excess, zfra
      !!-------------------------------------------------------------------
      !
      z1_dt = 1._wp / pdt
      !
      ! -- check snow load -- !
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         IF ( pv_i(ji,jj,jl) > 0._wp ) THEN
            !
            zvs_excess = MAX( 0._wp, pv_s(ji,jj,jl) - pv_i(ji,jj,jl) * (rho0-rhoi) * r1_rhos )
            !
            IF( zvs_excess > 0._wp ) THEN   ! snow-ice interface deplets below the ocean surface
               ! put snow excess in the ocean
               zfra = ( pv_s(ji,jj,jl) - zvs_excess ) / MAX( pv_s(ji,jj,jl), epsi20 )
               wfx_res(ji,jj) = wfx_res(ji,jj) + zvs_excess * rhos * z1_dt
               hfx_res(ji,jj) = hfx_res(ji,jj) - SUM( pe_s(ji,jj,1:nlay_s,jl) ) * ( 1._wp - zfra ) * z1_dt ! W.m-2 <0
               ! correct snow volume and heat content
               pe_s(ji,jj,1:nlay_s,jl) = pe_s(ji,jj,1:nlay_s,jl) * zfra
               pv_s(ji,jj,jl)          = pv_s(ji,jj,jl) - zvs_excess
            ENDIF
            !
         ENDIF
         END_2D
      END DO
      !
      !-- correct pond concentration to avoid a_ip > a_i -- !
      WHERE( pa_ip(:,:,:) > pa_i(:,:,:) )   pa_ip(:,:,:) = pa_i(:,:,:)
      !
   END SUBROUTINE Hsnow


   SUBROUTINE adv_pra_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE adv_pra_init  ***
      !!
      !! ** Purpose :   allocate and initialize arrays for Prather advection
      !!-------------------------------------------------------------------
      INTEGER ::   ierr
      CHARACTER(len=19), PARAMETER :: crtnm = 'adv_pra_init'
      !!-------------------------------------------------------------------
      !
      !                             !* allocate prather fields
      ALLOCATE( sxice(jpi,jpj,jpl) , syice(jpi,jpj,jpl) , sxxice(jpi,jpj,jpl) , syyice(jpi,jpj,jpl) , sxyice(jpi,jpj,jpl) ,   &
         &      sxsn (jpi,jpj,jpl) , sysn (jpi,jpj,jpl) , sxxsn (jpi,jpj,jpl) , syysn (jpi,jpj,jpl) , sxysn (jpi,jpj,jpl) ,   &
         &      sxa  (jpi,jpj,jpl) , sya  (jpi,jpj,jpl) , sxxa  (jpi,jpj,jpl) , syya  (jpi,jpj,jpl) , sxya  (jpi,jpj,jpl) ,   &
         &      sxsal(jpi,jpj,jpl) , sysal(jpi,jpj,jpl) , sxxsal(jpi,jpj,jpl) , syysal(jpi,jpj,jpl) , sxysal(jpi,jpj,jpl) ,   &
         &      sxage(jpi,jpj,jpl) , syage(jpi,jpj,jpl) , sxxage(jpi,jpj,jpl) , syyage(jpi,jpj,jpl) , sxyage(jpi,jpj,jpl) ,   &
         &      sxap (jpi,jpj,jpl) , syap (jpi,jpj,jpl) , sxxap (jpi,jpj,jpl) , syyap (jpi,jpj,jpl) , sxyap (jpi,jpj,jpl) ,   &
         &      sxvp (jpi,jpj,jpl) , syvp (jpi,jpj,jpl) , sxxvp (jpi,jpj,jpl) , syyvp (jpi,jpj,jpl) , sxyvp (jpi,jpj,jpl) ,   &
         &      sxvl (jpi,jpj,jpl) , syvl (jpi,jpj,jpl) , sxxvl (jpi,jpj,jpl) , syyvl (jpi,jpj,jpl) , sxyvl (jpi,jpj,jpl) ,   &
                                !
         &      sxc0 (jpi,jpj,nlay_s,jpl) , syc0 (jpi,jpj,nlay_s,jpl) , sxxc0(jpi,jpj,nlay_s,jpl) , &
         &      syyc0(jpi,jpj,nlay_s,jpl) , sxyc0(jpi,jpj,nlay_s,jpl)                             , &
                                !
         &      sxe  (jpi,jpj,nlay_i,jpl) , sye  (jpi,jpj,nlay_i,jpl) , sxxe (jpi,jpj,nlay_i,jpl) , &
         &      syye (jpi,jpj,nlay_i,jpl) , sxye (jpi,jpj,nlay_i,jpl)                             , &
         &      STAT = ierr )
      !
      CALL mpp_sum( crtnm, ierr )
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'adv_pra_init : unable to allocate ice arrays for Prather advection scheme')
      !
      CALL adv_pra_rst( 'READ' )    !* read or initialize all required files
      !
   END SUBROUTINE adv_pra_init


   SUBROUTINE adv_pra_rst( cdrw, kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE adv_pra_rst  ***
      !!
      !! ** Purpose :   Read or write file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!----------------------------------------------------------------------
      CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      INTEGER, OPTIONAL, INTENT(in) ::   kt     ! ice time-step
      !
      INTEGER ::   jk, jl   ! dummy loop indices
      INTEGER ::   iter     ! local integer
      INTEGER ::   id1      ! local integer
      CHARACTER(len=25) ::   znam
      CHARACTER(len=2)  ::   zchar, zchar1
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      !                                      !==========================!
      IF( TRIM(cdrw) == 'READ' ) THEN        !==  Read or initialize  ==!
         !                                   !==========================!
         !
         IF( ln_rstart ) THEN
            id1 = iom_varid( numrir, 'sxice' , ldstop = .FALSE. )    ! file exist: id1>0
         ELSE
            id1 = 0                                                  ! no restart: id1=0
         ENDIF
         !
         IF( id1 > 0 ) THEN                     !**  Read the restart file  **!
            !
            !                                                        ! ice thickness
            CALL iom_get( numrir, jpdom_auto, 'sxice' , sxice , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'syice' , syice , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxice', sxxice )
            CALL iom_get( numrir, jpdom_auto, 'syyice', syyice )
            CALL iom_get( numrir, jpdom_auto, 'sxyice', sxyice )
            !                                                        ! snow thickness
            CALL iom_get( numrir, jpdom_auto, 'sxsn'  , sxsn  , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sysn'  , sysn  , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxsn' , sxxsn  )
            CALL iom_get( numrir, jpdom_auto, 'syysn' , syysn  )
            CALL iom_get( numrir, jpdom_auto, 'sxysn' , sxysn  )
            !                                                        ! ice concentration
            CALL iom_get( numrir, jpdom_auto, 'sxa'   , sxa   , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sya'   , sya   , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxa'  , sxxa   )
            CALL iom_get( numrir, jpdom_auto, 'syya'  , syya   )
            CALL iom_get( numrir, jpdom_auto, 'sxya'  , sxya   )
            !                                                        ! ice salinity
            CALL iom_get( numrir, jpdom_auto, 'sxsal' , sxsal , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sysal' , sysal , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxsal', sxxsal )
            CALL iom_get( numrir, jpdom_auto, 'syysal', syysal )
            CALL iom_get( numrir, jpdom_auto, 'sxysal', sxysal )
            !                                                        ! ice age
            CALL iom_get( numrir, jpdom_auto, 'sxage' , sxage , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'syage' , syage , psgn = -1._wp )
            CALL iom_get( numrir, jpdom_auto, 'sxxage', sxxage )
            CALL iom_get( numrir, jpdom_auto, 'syyage', syyage )
            CALL iom_get( numrir, jpdom_auto, 'sxyage', sxyage )
            !                                                        ! snow layers heat content
            DO jk = 1, nlay_s
               WRITE(zchar1,'(I2.2)') jk
               znam = 'sxc0'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sxc0 (:,:,jk,:) = z3d(:,:,:)
               znam = 'syc0'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   syc0 (:,:,jk,:) = z3d(:,:,:)
               znam = 'sxxc0'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   sxxc0(:,:,jk,:) = z3d(:,:,:)
               znam = 'syyc0'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   syyc0(:,:,jk,:) = z3d(:,:,:)
               znam = 'sxyc0'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   sxyc0(:,:,jk,:) = z3d(:,:,:)
            END DO
            !                                                        ! ice layers heat content
            DO jk = 1, nlay_i
               WRITE(zchar1,'(I2.2)') jk
               znam = 'sxe'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sxe (:,:,jk,:) = z3d(:,:,:)
               znam = 'sye'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d, psgn = -1._wp )   ;   sye (:,:,jk,:) = z3d(:,:,:)
               znam = 'sxxe'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   sxxe(:,:,jk,:) = z3d(:,:,:)
               znam = 'syye'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   syye(:,:,jk,:) = z3d(:,:,:)
               znam = 'sxye'//'_l'//zchar1
               CALL iom_get( numrir, jpdom_auto, znam , z3d )   ;   sxye(:,:,jk,:) = z3d(:,:,:)
            END DO
            !
            IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN                                    ! melt pond fraction
               IF( iom_varid( numrir, 'sxap', ldstop = .FALSE. ) > 0 ) THEN
                  CALL iom_get( numrir, jpdom_auto, 'sxap' , sxap , psgn = -1._wp )
                  CALL iom_get( numrir, jpdom_auto, 'syap' , syap , psgn = -1._wp )
                  CALL iom_get( numrir, jpdom_auto, 'sxxap', sxxap )
                  CALL iom_get( numrir, jpdom_auto, 'syyap', syyap )
                  CALL iom_get( numrir, jpdom_auto, 'sxyap', sxyap )
                  !                                                     ! melt pond volume
                  CALL iom_get( numrir, jpdom_auto, 'sxvp' , sxvp , psgn = -1._wp )
                  CALL iom_get( numrir, jpdom_auto, 'syvp' , syvp , psgn = -1._wp )
                  CALL iom_get( numrir, jpdom_auto, 'sxxvp', sxxvp )
                  CALL iom_get( numrir, jpdom_auto, 'syyvp', syyvp )
                  CALL iom_get( numrir, jpdom_auto, 'sxyvp', sxyvp )
               ELSE
                  sxap = 0._wp ;   syap = 0._wp    ;   sxxap = 0._wp    ;   syyap = 0._wp    ;   sxyap = 0._wp   ! melt pond fraction
                  sxvp = 0._wp ;   syvp = 0._wp    ;   sxxvp = 0._wp    ;   syyvp = 0._wp    ;   sxyvp = 0._wp   ! melt pond volume
               ENDIF
               !
               IF ( ln_pnd_lids ) THEN                               ! melt pond lid volume
                  IF( iom_varid( numrir, 'sxvl', ldstop = .FALSE. ) > 0 ) THEN
                     CALL iom_get( numrir, jpdom_auto, 'sxvl' , sxvl , psgn = -1._wp )
                     CALL iom_get( numrir, jpdom_auto, 'syvl' , syvl , psgn = -1._wp )
                     CALL iom_get( numrir, jpdom_auto, 'sxxvl', sxxvl )
                     CALL iom_get( numrir, jpdom_auto, 'syyvl', syyvl )
                     CALL iom_get( numrir, jpdom_auto, 'sxyvl', sxyvl )
                  ELSE
                     sxvl = 0._wp; syvl = 0._wp    ;   sxxvl = 0._wp    ;   syyvl = 0._wp    ;   sxyvl = 0._wp   ! melt pond lid volume
                  ENDIF
               ENDIF
            ENDIF
            !
         ELSE                                   !**  start rheology from rest  **!
            !
            IF(lwp) WRITE(numout,*) '   ==>>   start from rest OR previous run without Prather, set moments to 0'
            !
            sxice = 0._wp   ;   syice = 0._wp   ;   sxxice = 0._wp   ;   syyice = 0._wp   ;   sxyice = 0._wp      ! ice thickness
            sxsn  = 0._wp   ;   sysn  = 0._wp   ;   sxxsn  = 0._wp   ;   syysn  = 0._wp   ;   sxysn  = 0._wp      ! snow thickness
            sxa   = 0._wp   ;   sya   = 0._wp   ;   sxxa   = 0._wp   ;   syya   = 0._wp   ;   sxya   = 0._wp      ! ice concentration
            sxsal = 0._wp   ;   sysal = 0._wp   ;   sxxsal = 0._wp   ;   syysal = 0._wp   ;   sxysal = 0._wp      ! ice salinity
            sxage = 0._wp   ;   syage = 0._wp   ;   sxxage = 0._wp   ;   syyage = 0._wp   ;   sxyage = 0._wp      ! ice age
            sxc0  = 0._wp   ;   syc0  = 0._wp   ;   sxxc0  = 0._wp   ;   syyc0  = 0._wp   ;   sxyc0  = 0._wp      ! snow layers heat content
            sxe   = 0._wp   ;   sye   = 0._wp   ;   sxxe   = 0._wp   ;   syye   = 0._wp   ;   sxye   = 0._wp      ! ice layers heat content
            IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
               sxap = 0._wp ;   syap = 0._wp    ;   sxxap = 0._wp    ;   syyap = 0._wp    ;   sxyap = 0._wp       ! melt pond fraction
               sxvp = 0._wp ;   syvp = 0._wp    ;   sxxvp = 0._wp    ;   syyvp = 0._wp    ;   sxyvp = 0._wp       ! melt pond volume
               IF ( ln_pnd_lids ) THEN
                  sxvl = 0._wp; syvl = 0._wp    ;   sxxvl = 0._wp    ;   syyvl = 0._wp    ;   sxyvl = 0._wp       ! melt pond lid volume
               ENDIF
            ENDIF
         ENDIF
         !
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
         CALL iom_rstput( iter, nitrst, numriw, 'sxice' , sxice  )
         CALL iom_rstput( iter, nitrst, numriw, 'syice' , syice  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxice', sxxice )
         CALL iom_rstput( iter, nitrst, numriw, 'syyice', syyice )
         CALL iom_rstput( iter, nitrst, numriw, 'sxyice', sxyice )
         !                                                           ! snow thickness
         CALL iom_rstput( iter, nitrst, numriw, 'sxsn'  , sxsn   )
         CALL iom_rstput( iter, nitrst, numriw, 'sysn'  , sysn   )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxsn' , sxxsn  )
         CALL iom_rstput( iter, nitrst, numriw, 'syysn' , syysn  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxysn' , sxysn  )
         !                                                           ! ice concentration
         CALL iom_rstput( iter, nitrst, numriw, 'sxa'   , sxa    )
         CALL iom_rstput( iter, nitrst, numriw, 'sya'   , sya    )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxa'  , sxxa   )
         CALL iom_rstput( iter, nitrst, numriw, 'syya'  , syya   )
         CALL iom_rstput( iter, nitrst, numriw, 'sxya'  , sxya   )
         !                                                           ! ice salinity
         CALL iom_rstput( iter, nitrst, numriw, 'sxsal' , sxsal  )
         CALL iom_rstput( iter, nitrst, numriw, 'sysal' , sysal  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxsal', sxxsal )
         CALL iom_rstput( iter, nitrst, numriw, 'syysal', syysal )
         CALL iom_rstput( iter, nitrst, numriw, 'sxysal', sxysal )
         !                                                           ! ice age
         CALL iom_rstput( iter, nitrst, numriw, 'sxage' , sxage  )
         CALL iom_rstput( iter, nitrst, numriw, 'syage' , syage  )
         CALL iom_rstput( iter, nitrst, numriw, 'sxxage', sxxage )
         CALL iom_rstput( iter, nitrst, numriw, 'syyage', syyage )
         CALL iom_rstput( iter, nitrst, numriw, 'sxyage', sxyage )
         !                                                           ! snow layers heat content
         DO jk = 1, nlay_s
            WRITE(zchar1,'(I2.2)') jk
            znam = 'sxc0'//'_l'//zchar1  ;   z3d(:,:,:) = sxc0 (:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'syc0'//'_l'//zchar1  ;   z3d(:,:,:) = syc0 (:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'sxxc0'//'_l'//zchar1 ;   z3d(:,:,:) = sxxc0(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'syyc0'//'_l'//zchar1 ;   z3d(:,:,:) = syyc0(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'sxyc0'//'_l'//zchar1 ;   z3d(:,:,:) = sxyc0(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
         END DO
         !                                                           ! ice layers heat content
         DO jk = 1, nlay_i
            WRITE(zchar1,'(I2.2)') jk
            znam = 'sxe'//'_l'//zchar1   ;   z3d(:,:,:) = sxe (:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'sye'//'_l'//zchar1   ;   z3d(:,:,:) = sye (:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'sxxe'//'_l'//zchar1  ;   z3d(:,:,:) = sxxe(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'syye'//'_l'//zchar1  ;   z3d(:,:,:) = syye(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
            znam = 'sxye'//'_l'//zchar1  ;   z3d(:,:,:) = sxye(:,:,jk,:)
            CALL iom_rstput( iter, nitrst, numriw, znam , z3d )
         END DO
         !
         IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN                                       ! melt pond fraction
            CALL iom_rstput( iter, nitrst, numriw, 'sxap' , sxap  )
            CALL iom_rstput( iter, nitrst, numriw, 'syap' , syap  )
            CALL iom_rstput( iter, nitrst, numriw, 'sxxap', sxxap )
            CALL iom_rstput( iter, nitrst, numriw, 'syyap', syyap )
            CALL iom_rstput( iter, nitrst, numriw, 'sxyap', sxyap )
            !                                                        ! melt pond volume
            CALL iom_rstput( iter, nitrst, numriw, 'sxvp' , sxvp  )
            CALL iom_rstput( iter, nitrst, numriw, 'syvp' , syvp  )
            CALL iom_rstput( iter, nitrst, numriw, 'sxxvp', sxxvp )
            CALL iom_rstput( iter, nitrst, numriw, 'syyvp', syyvp )
            CALL iom_rstput( iter, nitrst, numriw, 'sxyvp', sxyvp )
            !
            IF ( ln_pnd_lids ) THEN                                  ! melt pond lid volume
               CALL iom_rstput( iter, nitrst, numriw, 'sxvl' , sxvl  )
               CALL iom_rstput( iter, nitrst, numriw, 'syvl' , syvl  )
               CALL iom_rstput( iter, nitrst, numriw, 'sxxvl', sxxvl )
               CALL iom_rstput( iter, nitrst, numriw, 'syyvl', syyvl )
               CALL iom_rstput( iter, nitrst, numriw, 'sxyvl', sxyvl )
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE adv_pra_rst

   SUBROUTINE icemax3D( pice , pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax3D ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in ) ::   pice   ! input
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pmax   ! output
      !
      REAL(wp), DIMENSION(Nis0:Nie0) ::   zmax1, zmax2
      REAL(wp)                       ::   zmax3
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      !!----------------------------------------------------------------------
      ! basic version: get the max of epsi20 + 9 neighbours
      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      DO jl = 1, jpl
         DO ji = Nis0, Nie0
            zmax1(ji) = MAX( epsi20, pice(ji,Njs0-1,jl), pice(ji-1,Njs0-1,jl), pice(ji+1,Njs0-1,jl) )
            zmax2(ji) = MAX( epsi20, pice(ji,Njs0  ,jl), pice(ji-1,Njs0  ,jl), pice(ji+1,Njs0  ,jl) )
         END DO
         DO_2D( 0, 0, 0, 0 )
         zmax3 = MAX( epsi20, pice(ji,jj+1,jl), pice(ji-1,jj+1,jl), pice(ji+1,jj+1,jl) )
         pmax(ji,jj,jl) = MAX( epsi20, zmax1(ji), zmax2(ji), zmax3 )
         zmax1(ji) = zmax2(ji)
         zmax2(ji) = zmax3
         END_2D
      END DO
   END SUBROUTINE icemax3D

   SUBROUTINE icemax4D( pice , pmax )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE icemax4D ***
      !! ** Purpose :  compute the max of the 9 points around
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in ) ::   pice   ! input
      REAL(wp), DIMENSION(:,:,:,:), INTENT(out) ::   pmax   ! output
      !
      REAL(wp), DIMENSION(Nis0:Nie0) ::   zmax1, zmax2
      REAL(wp)                       ::   zmax3
      INTEGER  ::   jlay, ji, jj, jk, jl   ! dummy loop indices
      !!----------------------------------------------------------------------
      jlay = SIZE( pice , 3 )   ! size of input arrays
      ! basic version: get the max of epsi20 + 9 neighbours
      ! optimized version : does a little bit more than 2 max of epsi20 + 3 neighbours
      DO jl = 1, jpl
         DO jk = 1, jlay
            DO ji = Nis0, Nie0
               zmax1(ji) = MAX( epsi20, pice(ji,Njs0-1,jk,jl), pice(ji-1,Njs0-1,jk,jl), pice(ji+1,Njs0-1,jk,jl) )
               zmax2(ji) = MAX( epsi20, pice(ji,Njs0  ,jk,jl), pice(ji-1,Njs0  ,jk,jl), pice(ji+1,Njs0  ,jk,jl) )
            END DO
            DO_2D( 0, 0, 0, 0 )
            zmax3 = MAX( epsi20, pice(ji,jj+1,jk,jl), pice(ji-1,jj+1,jk,jl), pice(ji+1,jj+1,jk,jl) )
            pmax(ji,jj,jk,jl) = MAX( epsi20, zmax1(ji), zmax2(ji), zmax3 )
            zmax1(ji) = zmax2(ji)
            zmax2(ji) = zmax3
            END_2D
         END DO
      END DO
   END SUBROUTINE icemax4D

#else
   !!----------------------------------------------------------------------
   !!   Default option            Dummy module        NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv_pra
