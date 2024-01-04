MODULE icedyn_rhg_util
   !!======================================================================
   !!                     ***  MODULE  icedyn_rhg_util  ***
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
   USE lib_mpp
   USE lbclnk         ! lateral boundary conditions (or mpp links)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   strain_rate

   PUBLIC   div_stress_tensor
   PUBLIC   cap_damage        ! called by damage-advection routines
   !PUBLIC   fix_damage
   PUBLIC   rmpT2F
   PUBLIC   rmpF2T
   PUBLIC   rmpU2V
   PUBLIC   rmpV2U

   !PUBLIC   rmpFT2T
   !PUBLIC   rmpTF2F

   PUBLIC   smooth5pT
   !PUBLIC   smooth9p

   PUBLIC clean_small_a_all
   PUBLIC clean_small_a_sgm

   REAL(wp), PARAMETER :: rtol_dmg = 0.1_wp   ! tolerance for damage overshoot (above/below 1/0)

   REAL(wp), PARAMETER :: rclean_below_A = 0.1_wp
   !REAL(wp), PARAMETER :: rclean_below_A = 0.4_wp
   !REAL(wp), PARAMETER :: rclean_below_A = 0.01_wp

   !REAL(wp), PARAMETER :: rrhv_dmg = 0.8_wp

   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE cap_damage( ktb, kts,  cgt, crtn,     pd )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cap_damage  ***
      !!
      !! ** Purpose :   Cap damage and report worying overshoots!
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER,                  INTENT(in)    :: ktb, kts  ! ice time-steps
      CHARACTER(len=1),         INTENT(in)    :: cgt       ! grid point ('T','F')
      CHARACTER(len=*),         INTENT(in)    :: crtn      ! name of routine it is called from !
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pd  ! damage
      !CHARACTER(len=*) , INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      LOGICAL :: l_bad_overshoot, l_bad_undrshoot
      !INTEGER  ::   iter            ! local integer
      !INTEGER  ::   id1, id2, id3   ! local integers
      !!----------------------------------------------------------------------
      IF( cgt == 'T' ) pd(:,:) = pd(:,:) * xmskt(:,:)
      IF( cgt == 'F' ) pd(:,:) = pd(:,:) * xmskf(:,:)

      l_bad_overshoot = ANY( (pd(2:jpi-1,2:jpj-1) >= rn_dmg_max + rtol_dmg) )
      l_bad_undrshoot = ANY( (pd(2:jpi-1,2:jpj-1) <=    0._wp   - rtol_dmg) )

      IF( l_bad_overshoot .OR. l_bad_undrshoot ) THEN
         !! Enter investigation chain:
         IF( l_bad_overshoot ) CALL ctl_warn( 'WARNING', ' "'//TRIM(crtn)//'" => Bad overshoot  for damage @ '//cgt//'-points!' )
         IF( l_bad_undrshoot ) CALL ctl_warn( 'WARNING', ' "'//TRIM(crtn)//'" => Bad undershoot for damage @ '//cgt//'-points!' )
      END IF

      !! The fast way:
      pd(:,:) = MIN( MAX( pd(:,:), 0._wp ) , rn_dmg_max )

   END SUBROUTINE cap_damage


   FUNCTION rmpT2F( pxt,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpT2F
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxt
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv=.FALSE.
      !!
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpT2F(:,:) = 0._wp
      !!
      DO_2D( 1,0, 1,0 )
            !!
            i1 = ji   ; j1 = jj
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj+1
            i4 = ji+1 ; j4 = jj+1
            !!
            zt1 = pxt(i1,j1)*xmskt(i1,j1)
            zt2 = pxt(i2,j2)*xmskt(i2,j2)
            zt3 = pxt(i3,j3)*xmskt(i3,j3)
            zt4 = pxt(i4,j4)*xmskt(i4,j4)
            zfc = xmskf(ji,jj)
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2t(i1,j1)
               zt2 = zt2 * e1e2t(i2,j2)
               zt3 = zt3 * e1e2t(i3,j3)
               zt4 = zt4 * e1e2t(i4,j4)
               zfc = zfc * r1_e1e2f(ji,jj)
            END IF
            !!
            zm = xmskt(i1,j1) + xmskt(i2,j2) + xmskt(i3,j3) + xmskt(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpT2F(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
      END_2D
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpT2F@icedyn_rhg_bbm', rmpT2F, 'F', 1._wp )
      END IF
      !!
   END FUNCTION rmpT2F


   FUNCTION rmpF2T( pxf, lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpF2T
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxf
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zf1, zf2, zf3, zf4, zs, zm, zz, zfc
      LOGICAL  :: lcnsrv=.FALSE.
      !!
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpF2T(:,:) = 0._wp
      !!
      DO_2D( 0,1, 0,1 )
            !!
            i1 = ji   ; j1 = jj
            i2 = ji-1 ; j2 = jj
            i3 = ji-1 ; j3 = jj-1
            i4 = ji   ; j4 = jj-1
            !!
            zf1 = pxf(i1,j1)*xmskf(i1,j1)
            zf2 = pxf(i2,j2)*xmskf(i2,j2)
            zf3 = pxf(i3,j3)*xmskf(i3,j3)
            zf4 = pxf(i4,j4)*xmskf(i4,j4)
            zfc = xmskt(ji,jj)
            IF( lcnsrv ) THEN
               zf1 = zf1 * e1e2f(i1,j1)
               zf2 = zf2 * e1e2f(i2,j2)
               zf3 = zf3 * e1e2f(i3,j3)
               zf4 = zf4 * e1e2f(i4,j4)
               zfc = zfc * r1_e1e2t(ji,jj)
            END IF
            !!
            zm = xmskf(i1,j1) + xmskf(i2,j2) + xmskf(i3,j3) + xmskf(i4,j4)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet F-point, `0` otherwize
            zfc = zfc * zz
            !
            zs = MAX( zm , 1.E-12_wp )
            !!
            rmpF2T(ji,jj) = ( zf1 + zf2 + zf3 + zf4 ) * zfc / zs
      END_2D
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpF2T@icedyn_rhg_bbm', rmpF2T, 'T', 1._wp )
      END IF
      !!
   END FUNCTION rmpF2T


   FUNCTION rmpU2V( pxu,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpU2V
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxu
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv=.FALSE.
      !!
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpU2V(:,:) = 0._wp
      !!
      DO_2D( 0,1, 1,0 )
            !!
            i1 = ji   ; j1 = jj
            i2 = ji   ; j2 = jj+1
            i3 = ji-1 ; j3 = jj+1
            i4 = ji-1 ; j4 = jj
            !!
            zt1 = pxu(i1,j1)*umask(i1,j1,1)
            zt2 = pxu(i2,j2)*umask(i2,j2,1)
            zt3 = pxu(i3,j3)*umask(i3,j3,1)
            zt4 = pxu(i4,j4)*umask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2u(i1,j1)
               zt2 = zt2 * e1e2u(i2,j2)
               zt3 = zt3 * e1e2u(i3,j3)
               zt4 = zt4 * e1e2u(i4,j4)
               zfc = zfc * r1_e1e2v(ji,jj)
            END IF
            !!
            zm = umask(i1,j1,1) + umask(i2,j2,1) + umask(i3,j3,1) + umask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpU2V(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
      END_2D
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpU2V@icedyn_rhg_bbm', rmpU2V, 'V', 1._wp )
      END IF
      !!
   END FUNCTION rmpU2V


   FUNCTION rmpV2U( pxv,  lbcl, lconserv )
      REAL(wp), DIMENSION(jpi,jpj)             :: rmpV2U
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxv
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl, lconserv
      !!
      INTEGER  :: ji, jj, i1, j1, i2, j2, i3, j3, i4, j4
      REAL(wp) :: zt1, zt2, zt3, zt4, zs, zfc, zm, zz
      LOGICAL  :: lcnsrv=.FALSE.
      !!
      IF( PRESENT( lconserv ) ) lcnsrv = lconserv
      !!
      rmpV2U(:,:) = 0._wp
      !!
      DO_2D( 1,0, 0,1 )
            !!
            i1 = ji+1 ; j1 = jj-1
            i2 = ji+1 ; j2 = jj
            i3 = ji   ; j3 = jj
            i4 = ji   ; j4 = jj-1
            !!
            zt1 = pxv(i1,j1)*vmask(i1,j1,1)
            zt2 = pxv(i2,j2)*vmask(i2,j2,1)
            zt3 = pxv(i3,j3)*vmask(i3,j3,1)
            zt4 = pxv(i4,j4)*vmask(i4,j4,1)
            zfc = 1._wp
            IF( lcnsrv ) THEN
               zt1 = zt1 * e1e2v(i1,j1)
               zt2 = zt2 * e1e2v(i2,j2)
               zt3 = zt3 * e1e2v(i3,j3)
               zt4 = zt4 * e1e2v(i4,j4)
               zfc = zfc * r1_e1e2u(ji,jj)
            END IF
            !!
            zm = vmask(i1,j1,1) + vmask(i2,j2,1) + vmask(i3,j3,1) + vmask(i4,j4,1)
            zz = MIN( zm , 1._wp ) ! => `1` if at least a surrounding wet T-point, `0` otherwize
            zfc = zfc * zz
            !!
            zs  = MAX( zm , 1.E-12_wp )
            !!
            rmpV2U(ji,jj) = ( zt1 + zt2 + zt3 + zt4 ) * zfc / zs
            !!
      END_2D
      !!
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'rmpV2U@icedyn_rhg_bbm', rmpV2U, 'U', 1._wp )
      END IF
      !!
   END FUNCTION rmpV2U



   FUNCTION smooth5pT( px, pwij, lbcl )
      !!
      !! Smooth T field by using local value + 4 surounding F points and 4 surounding T points
      REAL(wp), DIMENSION(jpi,jpj)             :: smooth5pT
      !CHARACTER(len=1),         INTENT(in)     :: cgt       ! grid point ('T','F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl !, lconserv
      !!
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4
      REAL(wp) :: zs, zaa, zbb, zfc
      !===================================================================
      !
      zaa = pwij
      zbb = 1._wp - zaa
      smooth5pT(:,:) = 0._wp
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            it1 = ji+1 ; jt1 = jj
            it2 = ji   ; jt2 = jj+1
            it3 = ji-1 ; jt3 = jj
            it4 = ji   ; jt4 = jj-1
            !!
            zm0 = xmskt( ji,jj ) * e1e2t(ji,jj)
            zm1 = xmskt(it1,jt1) * e1e2t(it1,jt1)
            zm2 = xmskt(it2,jt2) * e1e2t(it2,jt2)
            zm3 = xmskt(it3,jt3) * e1e2t(it3,jt3)
            zm4 = xmskt(it4,jt4) * e1e2t(it4,jt4)
            !!
            zt0 = px( ji,jj )*zm0
            zt1 = px(it1,jt1)*zm1
            zt2 = px(it2,jt2)*zm2
            zt3 = px(it3,jt3)*zm3
            zt4 = px(it4,jt4)*zm4
            !!
            zs     =        MAX( zaa*zm0  +  zbb*( zm1 + zm2 + zm3 + zm4 ) , 1.E-12_wp ) ! sum of wheights
            !
            smooth5pT(ji,jj) = (  zaa*zt0  +  zbb*( zt1 + zt2 + zt3 + zt4 ) ) / zs * xmskt(ji,jj)
            !
      END_2D
      !
      IF(PRESENT(lbcl)) THEN
         IF( lbcl ) CALL lbc_lnk( 'smooth5pT@icedyn_rhg_bbm', smooth5pT, 'T', 1._wp )
      END IF
      !
   END FUNCTION smooth5pT


   SUBROUTINE strain_rate( cgt, pU, pV, pUd, pVd, p1_e1e2, pe2X, pe1Y, p1_e2X, p1_e1Y, pe1e1, pe2e2, &
      &                    pmask, pdudx, pdvdy, pshr,                                                &
      &                    lblnk, pdudy, pdvdx, pdiv, pmaxshr )
      !!
      !! Computes the 3 elements of the strain rate tensor, e11, e22 & e12, at either T- or F-points
      !!
      !! Note: when dealing with F-points (cgt='F'), `pmask` must be the actual `fmask` that takes into
      !!       condition the slip/no-slip conditions
      !!       (important for shear strain: `pshr`, `pdudy`, `pdvdx` and `pmaxshr` !)
      !!
      CHARACTER(len=1),         INTENT(in)  :: cgt              ! grid point type: 'T' or 'F'
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pU, pV, pUd, pVd ! u,v of T-point, u,v of F-point                    [m/s]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2          ! T-grid: 1/(e1t*e2t) | F-grid: 1/(e1f*e2f)         [1/m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2X, pe1Y       ! T-grid: e2u,e1v | F-grid: e2v,e1u                 [m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2X, p1_e1Y   ! T-grid: 1/e2u,1/e1v | F-grid: 1/e2v,1/e1u         [1/m]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2     ! T-grid: e1t*e1t,e2t*e2t | F-grid: e1f*e1f,e2f*e2f [m^2]
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pmask            ! 2D land-sea mask for given points...
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdudx, pdvdy, pshr  ! e11, e22 & e12 @ `cgt` points                   [1/s]
      !!
      LOGICAL , OPTIONAL,                 INTENT(in)  :: lblnk
      REAL(wp), OPTIONAL, DIMENSION(:,:), INTENT(out) :: pdudy, pdvdx, pdiv, pmaxshr ! @ `cgt` points                [1/s]
      !!
      LOGICAL  :: l_b_lnk=.FALSE., l_rtrn_dudy, l_rtrn_dvdx, l_rtrn_div, l_rtrn_maxshr
      REAL(wp) :: zE1, zE2, zS1, zS2, z1_e1e2, zzf, ze2e2, ze1e1, zmask
      INTEGER  :: ip, im, jp, jm, ji, jj, k1, k2
      !!
      IF( PRESENT(lblnk) ) l_b_lnk = lblnk

      l_rtrn_dudy = PRESENT( pdudy )
      l_rtrn_dvdx = PRESENT( pdvdx )
      l_rtrn_div  = PRESENT(  pdiv )
      l_rtrn_maxshr = PRESENT( pmaxshr )

      IF ( cgt == 'T' ) THEN
         !! In T-centric cell: dU/dX @ T-point = (U(i,j) - U(i-1,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  0
         im = -1
         jp =  0
         jm = -1
         k1 =  0       ! => DO_2D( 0,1 , 0,1 )
         k2 =  1
      ELSEIF ( cgt == 'F' ) THEN
         !! In F-centric cell: dU/dX @ F-point = (U(i+1,j) - U(i,j))/dx == (U(i+ip,j) - U(i+im,j))/dx
         ip =  1
         im =  0
         jp =  1
         jm =  0
         k1 =  1       ! => DO_2D( 1,0 , 1,0 )
         k2 =  0
      ELSE
         CALL ctl_stop( 'STOP', 'strain_rate(): unknown grid-point type: '//cgt//'!')
      ENDIF

      DO_2D( k1, k2, k1, k2 )

            zmask = pmask(ji,jj)        ! actual mask containing right values for shear boundary conditions

            z1_e1e2 = p1_e1e2(ji,jj) * MIN(zmask, 1._wp)

            ze1e1 = pe1e1(ji,jj)
            ze2e2 = pe2e2(ji,jj)

            !! Divergence at cgt-points, `dU/dx + dV/dy` :
            zE1 = (   pe2X(ji+ip,jj)*pU(ji+ip,jj) - pe2X(ji+im,jj)*pU(ji+im,jj) &
               &    + pe1Y(ji,jj+jp)*pV(ji,jj+jp) - pe1Y(ji,jj+jm)*pV(ji,jj+jm) &
               &  ) * z1_e1e2
            IF( l_rtrn_div ) pdiv(ji,jj)  = zE1

            !! Tension at cgt-points, `dU/dx - dV/dy` :
            zE2 = (  ( pU(ji+ip,jj)*p1_e2X(ji+ip,jj) - pU(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 &
               &    -( pV(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pV(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 &
               &  ) * z1_e1e2

            pdudx(ji,jj) = 0.5_wp * ( zE1 + zE2 )
            pdvdy(ji,jj) = 0.5_wp * ( zE1 - zE2 )

            !! 2 * shear at cgt-points, `dU/dy + dV/dx` :
            zzf = z1_e1e2 * zmask
            zS1 = ( pUd(ji,jj+jp)*p1_e1Y(ji,jj+jp) - pUd(ji,jj+jm)*p1_e1Y(ji,jj+jm) ) * ze1e1 * zzf
            zS2 = ( pVd(ji+ip,jj)*p1_e2X(ji+ip,jj) - pVd(ji+im,jj)*p1_e2X(ji+im,jj) ) * ze2e2 * zzf

            pshr(ji,jj) = 0.5_wp * ( zS1 + zS2 )

            IF( l_rtrn_dudy ) pdudy(ji,jj) = zS1
            IF( l_rtrn_dvdx ) pdvdx(ji,jj) = zS2
            IF( l_rtrn_maxshr ) THEN
               zzf = zS1 + zS2
               pmaxshr(ji,jj)  = SQRT( zE2*zE2 + zzf*zzf )
            ENDIF

      END_2D

      IF( l_b_lnk ) THEN
         CALL lbc_lnk( 'strain_rate@icedyn_adv', pdudx,cgt,1., pdvdy,cgt,1., pshr,cgt,1. )
         !! Could be optimized (gathered) for configuration often used! #fixme!
         IF(l_rtrn_dudy ) CALL lbc_lnk( 'strain_rate@icedyn_adv', pdudy,cgt,1. )
         IF(l_rtrn_dvdx ) CALL lbc_lnk( 'strain_rate@icedyn_adv', pdvdx,cgt,1. )
         IF( l_rtrn_div ) CALL lbc_lnk( 'strain_rate@icedyn_adv',  pdiv,cgt,1. )
         IF( l_rtrn_maxshr ) CALL lbc_lnk( 'strain_rate@icedyn_adv', pmaxshr,cgt,1. )
      ENDIF

   END SUBROUTINE strain_rate



   SUBROUTINE div_stress_tensor( cgt, pe1e1, pe2e2,  pe1e1_d, pe2e2_d,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2X, p1_e1e2Y,  &
      &                               ps11h, ps22h, ps12h,  pdivSx, pdivSy )
      !!----------------------------------------------------------------------------------------------
      !! Computes the vector (pdivSx,pdivSy) = divergence of the h-integrated internal stress tensor
      !!
      !!   depending on the grid: T-centric grid => cgt='T' or F-centric grid => cgt='F'
      !!
      !! INPUT:                                               |     cgt=='T'   |    cgt=='F'    |
      !!   * ps11h, ps22h: sigma11*h, sigma22*h           =>  ! @ point T[i,j] | @ point F[i,j] |
      !!   * ps12h       :       sigma12*h                =>  ! @ point F[i,j] | @ point T[i,j] |
      !!
      !! RETURNS:                                             |     cgt=='T'   |    cgt=='F'    |
      !!   * pdivSx: x-component of the div of the tensor =>  | @ point U[i,j] | @ point V[i,j] |
      !!   * pdivSy: y-component of the div of the tensor =>  | @ point V[i,j] | @ point U[i,j] |
      !!
      !!----------------------------------------------------------------------------------------------
      CHARACTER(len=1),         INTENT(in)  :: cgt
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2, pe1e1_d, pe2e2_d
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2X, p1_e1e2Y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11h, ps22h, ps12h ! components of stress tensors on T- or F-centric grids x h !!!
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdivSx, pdivSy      ! x,y components of the divergence of the tensor
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zt1, zt2, zt3, zt4
      INTEGER  :: ip, im, jp, jm, ji, jj
      !!--------------------------------------------------------------------------------------------
      IF ( cgt == 'T' ) THEN
         ip =  1
         im =  0
         jp =  0
         jm = -1
      ELSEIF ( cgt == 'F' ) THEN
         ip =  0
         im = -1
         jp =  1
         jm =  0
      ELSE
         CALL ctl_stop( 'STOP', 'div_stress_tensor(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      pdivSx(:,:) = 0._wp
      pdivSy(:,:) = 0._wp
      !
      zt1(:,:) = ps11h(:,:) * pe2e2(:,:)
      zt2(:,:) = ps22h(:,:) * pe1e1(:,:)
      !
      zt3(:,:) = ps12h(:,:) * pe1e1_d(:,:)
      zt4(:,:) = ps12h(:,:) * pe2e2_d(:,:)
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
         !                   !--- ds11/dx + ds12/dy
         pdivSx(ji,jj) = ( ( zt1(ji+ip,jj) - zt1(ji+im,jj) ) * p1_e2x(ji,jj) &
            &            + ( zt3(ji,jj+jp) - zt3(ji,jj+jm) ) * p1_e1x(ji,jj) &
            &                 ) * p1_e1e2X(ji,jj)
         !                   !--- ds22/dy + ds12/dx
         pdivSy(ji,jj) = ( ( zt2(ji,jj-jm) - zt2(ji,jj-jp) ) * p1_e1y(ji,jj) &
            &            + ( zt4(ji-im,jj) - zt4(ji-ip,jj) ) * p1_e2y(ji,jj) &
            &                 ) * p1_e1e2Y(ji,jj)
         !
      END_2D
      !
      !! To understand:
      !
      !! T-centric grid:
      !! ---------------
      !DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
      !   !                   !--- Uu points: ds11/dx + ds12/dy
      !   pdivSx(ji,jj) = ( ( zt1(ji+1,jj) - zt1(ji,jj  ) ) * r1_e2u(ji,jj) &
      !      &          + ( zt3(ji  ,jj) - zt3(ji,jj-1) ) * r1_e1u(ji,jj) &
      !      &                 ) * r1_e1e2u(ji,jj)
      !   !                   !--- Vv points: ds22/dy + ds12/dx
      !   pdivSy(ji,jj) = ( ( zt2(ji,jj+1) - zt2(ji,jj  ) ) * r1_e1v(ji,jj) &
      !      &          + ( zt4(ji,jj  ) - zt4(ji-1,jj) ) * r1_e2v(ji,jj) &
      !      &                 ) * r1_e1e2v(ji,jj)
      !   !
      !END_2D
      !
      !! F-centric grid:
      !! ---------------
      !DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
      !   !                   !--- U@v points ( U contribution @ V-points ): ds11/dx + ds12/dy
      !   pdivSx(ji,jj) = ( ( zt1(ji,jj  ) - zt1(ji-1,jj) ) * r1_e2v(ji,jj) &
      !      &            + ( zt3(ji,jj+1) - zt3(ji,jj  ) ) * r1_e1v(ji,jj) &
      !      &                 ) * r1_e1e2v(ji,jj)
      !   !                   !--- V@u points ( V contribution @ U-points ): ds22/dy + ds12/dx
      !   pdivSy(ji,jj) = ( ( zt2(ji  ,jj) - zt2(ji,jj-1) ) * r1_e1u(ji,jj) &
      !      &            + ( zt4(ji+1,jj) - zt4(ji,jj  ) ) * r1_e2u(ji,jj) &
      !      &                 ) * r1_e1e2u(ji,jj)
      !   !
      !END_2D
      !
   END SUBROUTINE div_stress_tensor


   SUBROUTINE clean_small_a_all( pAt, pAf,  pdmgt, pdmgf,  ps11t, ps22t, ps12t,  ps11f, ps22f, ps12f )
      !!
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pAt, pAf             ! ice concentration @T and @F
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pdmgt, pdmgf         ! ice damage @T and @F
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11t, ps22t, ps12t  ! T-centric Sigmas [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11f, ps22f, ps12f  ! F-centric Sigmas [Pa]
      !!
      WHERE( pAt(:,:) < rclean_below_A )
         pdmgt(:,:) = 0._wp
         ps11t(:,:) = 0._wp
         ps22t(:,:) = 0._wp
         ps12t(:,:) = 0._wp
      END WHERE
      WHERE( pAf(:,:) < rclean_below_A )
         pdmgf(:,:) = 0._wp
         ps11f(:,:) = 0._wp
         ps22f(:,:) = 0._wp
         ps12f(:,:) = 0._wp
      END WHERE
      !!
   END SUBROUTINE clean_small_a_all

   SUBROUTINE clean_small_a_sgm( cgt, pAt, pAf,  ps11, ps22, ps12 )
      !!
      !! => clean only the 3 components of the `cgt`-centric stress tensor
      !!   ==> so either `s11t,s22t,s12f` of `s11f,s22f,s12t` !!!
      !!   ==> `ps12` is not at the same point as `ps11` & `ps22` by definition !!!
      !!
      CHARACTER(len=1),         INTENT(in)    :: cgt                  ! For the `cgt`-centric stress tensor (cgt:'T','F')
      REAL(wp), DIMENSION(:,:), INTENT(in)    :: pAt, pAf             ! ice concentration @T and @F
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: ps11, ps22, ps12
      !!
      IF(cgt == 'T') THEN
         WHERE( pAt(:,:) < rclean_below_A )
            ps11(:,:) = 0._wp
            ps22(:,:) = 0._wp
         ENDWHERE
         WHERE( pAf(:,:) < rclean_below_A ) ps12(:,:) = 0._wp
         !!
      ELSEIF(cgt == 'F') THEN
         WHERE( pAf(:,:) < rclean_below_A )
            ps11(:,:) = 0._wp
            ps22(:,:) = 0._wp
         ENDWHERE
         WHERE( pAt(:,:) < rclean_below_A ) ps12(:,:) = 0._wp
         !!
      ENDIF
      !!
   END SUBROUTINE clean_small_a_sgm

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_rhg_util

