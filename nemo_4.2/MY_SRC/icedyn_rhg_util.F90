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
   !!    cap_damage:
   !!    rmpT2F:
   !!    ...
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
   PUBLIC   div_stress_tensor_v2
   PUBLIC   cap_damage        ! called by damage-advection routines
   PUBLIC   rmpT2F
   PUBLIC   rmpF2T
   PUBLIC   rmpU2V
   PUBLIC   rmpV2U

   !PUBLIC   rmpFT2T
   !PUBLIC   rmpTF2F

   PUBLIC   smooth5p
   PUBLIC   smooth9p

   PUBLIC   smoothCrossTF
   
   PUBLIC clean_small_a_all
   PUBLIC clean_small_a_sgm

   PUBLIC sigmaII

   REAL(wp), PARAMETER, PUBLIC :: rclean_below_A = 0.01_wp
   
   REAL(wp), PARAMETER :: rtol_dmg = 0.1_wp   ! tolerance for damage overshoot (above/below 1/0)   

   !REAL(wp), PARAMETER :: rrhv_dmg = 0.8_wp

   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE cap_damage( cgt, crtn,     pd )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE cap_damage  ***
      !!
      !! ** Purpose :   Cap damage and report worying overshoots!
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      CHARACTER(len=1),         INTENT(in)    :: cgt       ! grid point ('T','F')
      CHARACTER(len=*),         INTENT(in)    :: crtn      ! name of routine it is called from !
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: pd  ! damage
      !!
      LOGICAL :: l_bad_overshoot, l_bad_undrshoot
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



   FUNCTION smooth5p( cgt, px, pwij,  lbcl )
      REAL(wp), DIMENSION(jpi,jpj) :: smooth5p
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt       ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4
      REAL(wp) :: zs, zaa, zbb, zfc
      LOGICAL  :: l_bcl
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      smooth5p(:,:) = 0._wp
      !
      IF    ( cgt=='T') THEN
         zma(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSEIF( cgt=='F') THEN
         zma(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smooth5p(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            it1 = ji+1 ; jt1 = jj
            it2 = ji   ; jt2 = jj+1
            it3 = ji-1 ; jt3 = jj
            it4 = ji   ; jt4 = jj-1
            !!
            zm0 = zma( ji,jj )
            zm1 = zma(it1,jt1)
            zm2 = zma(it2,jt2)
            zm3 = zma(it3,jt3)
            zm4 = zma(it4,jt4)
            !!
            zt0 = px( ji,jj )*zm0
            zt1 = px(it1,jt1)*zm1
            zt2 = px(it2,jt2)*zm2
            zt3 = px(it3,jt3)*zm3
            zt4 = px(it4,jt4)*zm4
            !!
            zs     =        MAX( zaa*zm0  +  zbb*( zm1 + zm2 + zm3 + zm4 ) , 1.E-12_wp ) ! sum of wheights
            !
            smooth5p(ji,jj) = (  zaa*zt0  +  zbb*( zt1 + zt2 + zt3 + zt4 ) ) / zs
            !
      END_2D
      !
      IF( cgt=='T') THEN
         smooth5p(:,:) = smooth5p(:,:)*xmskt(:,:)
      ELSE
         smooth5p(:,:) = smooth5p(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smooth5p@icedyn_rhg_bbm', smooth5p, cgt, 1._wp )
      !
   END FUNCTION smooth5p


   FUNCTION smooth9p( cgt, px, pwij, lbcl )
      REAL(wp), DIMENSION(jpi,jpj)             :: smooth9p
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt       ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4
      INTEGER  :: it5, jt5, it6, jt6, it7, jt7, it8, jt8
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4, zm5, zm6, zm7, zm8
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8
      REAL(wp) :: zs, zaa, zbb, zd
      LOGICAL  :: l_bcl
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      zd  = 0.7071067811865475_wp ! 1/sqrt(2)
      !
      IF    ( cgt=='T') THEN
         zma(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSEIF( cgt=='F') THEN
         zma(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smooth9p(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      smooth9p(:,:) = 0._wp
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            it1 = ji+1 ; jt1 = jj
            it2 = ji   ; jt2 = jj+1
            it3 = ji-1 ; jt3 = jj
            it4 = ji   ; jt4 = jj-1
            !!
            it5 = ji+1 ; jt5 = jj+1
            it6 = ji-1 ; jt6 = jj+1
            it7 = ji-1 ; jt7 = jj-1
            it8 = ji+1 ; jt8 = jj-1
            !!
            zm0 = zma( ji,jj )
            zm1 = zma(it1,jt1) 
            zm2 = zma(it2,jt2) 
            zm3 = zma(it3,jt3) 
            zm4 = zma(it4,jt4)
            !!
            zm5 = zma(it5,jt5) * zd
            zm6 = zma(it6,jt6) * zd
            zm7 = zma(it7,jt7) * zd
            zm8 = zma(it8,jt8) * zd
            !!
            zt0 = px( ji,jj )*zm0
            zt1 = px(it1,jt1)*zm1 ; zt5 = px(it5,jt5)*zm5
            zt2 = px(it2,jt2)*zm2 ; zt6 = px(it6,jt6)*zm6
            zt3 = px(it3,jt3)*zm3 ; zt7 = px(it7,jt7)*zm7
            zt4 = px(it4,jt4)*zm4 ; zt8 = px(it8,jt8)*zm8
            !
            zs     =        MAX( zaa * zm0  +  zbb * ( zm1 + zm2 + zm3 + zm4 + zm5 + zm6 + zm7 + zm8 ) , 1.E-12_wp ) ! sum of wheights
            !
            smooth9p(ji,jj) = (  zaa * zt0  +  zbb * ( zt1 + zt2 + zt3 + zt4 + zt5 + zt6 + zt7 + zt8 ) ) / zs
            !
      END_2D
      !
      IF( cgt=='T') THEN
         smooth9p(:,:) = smooth9p(:,:)*xmskt(:,:)
      ELSE
         smooth9p(:,:) = smooth9p(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smooth9p@icedyn_rhg_bbm', smooth9p, cgt, 1._wp )
      !
   END FUNCTION smooth9p



   FUNCTION smoothCrossTF( cgt, px, pxE, pwij,  lbcl )
      REAL(wp), DIMENSION(jpi,jpj) :: smoothCrossTF
      !!
      CHARACTER(len=1),             INTENT(in) :: cgt   ! grid ('T' or 'F')
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: px    ! field on `cgt` grid
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pxE   ! counterpart field !
      REAL(wp)                    , INTENT(in) :: pwij
      LOGICAL,            OPTIONAL, INTENT(in) :: lbcl
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zma0, zmaE
      INTEGER  :: ji, jj
      INTEGER  :: it1, jt1, it2, jt2, it3, jt3, it4, jt4, idxp
      REAL(wp) :: zm0, zm1, zm2, zm3, zm4
      REAL(wp) :: zt0, zt1, zt2, zt3, zt4
      REAL(wp) :: zs, zaa, zbb, zfc
      LOGICAL  :: l_bcl
      !===================================================================
      l_bcl = .FALSE.
      IF(PRESENT(lbcl)) l_bcl = lbcl
      !
      zaa = pwij
      zbb = 1._wp - zaa
      smoothCrossTF(:,:) = 0._wp
      !
      IF    ( cgt=='T') THEN
         idxp = 0
         zma0(:,:) = xmskt(:,:)*e1e2t(:,:)
         zmaE(:,:) = xmskf(:,:)*e1e2f(:,:)
      ELSEIF( cgt=='F') THEN
         idxp = 1
         zma0(:,:) = xmskf(:,:)*e1e2f(:,:)
         zmaE(:,:) = xmskt(:,:)*e1e2t(:,:)
      ELSE
         CALL ctl_stop( 'STOP', 'smoothCrossTF(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
            it1 = ji+idxp   ; jt1 = jj+idxp
            it2 = ji-1+idxp ; jt2 = jj+idxp
            it3 = ji-1+idxp ; jt3 = jj-1+idxp
            it4 = ji+idxp   ; jt4 = jj-1+idxp
            !IF on F-grid:
            !it1 = ji+1 ; jt1 = jj+1
            !it2 = ji   ; jt2 = jj+1
            !it3 = ji   ; jt3 = jj
            !it4 = ji+1 ; jt4 = jj
            !!
            zm0 = zma0( ji,jj )
            zm1 = zmaE(it1,jt1)
            zm2 = zmaE(it2,jt2)
            zm3 = zmaE(it3,jt3)
            zm4 = zmaE(it4,jt4)
            !!
            zt0 =  px( ji,jj )*zm0
            zt1 = pxE(it1,jt1)*zm1
            zt2 = pxE(it2,jt2)*zm2
            zt3 = pxE(it3,jt3)*zm3
            zt4 = pxE(it4,jt4)*zm4
            !!
            zs     =        MAX( zaa*zm0  +  zbb*( zm1 + zm2 + zm3 + zm4 ) , 1.E-12_wp ) ! sum of wheights
            !
            smoothCrossTF(ji,jj) = (  zaa*zt0  +  zbb*( zt1 + zt2 + zt3 + zt4 ) ) / zs
            !
      END_2D
      !
      IF( cgt=='T') THEN
         smoothCrossTF(:,:) = smoothCrossTF(:,:)*xmskt(:,:)
      ELSE
         smoothCrossTF(:,:) = smoothCrossTF(:,:)*xmskf(:,:)
      ENDIF
      !
      IF(l_bcl) CALL lbc_lnk( 'smoothCrossTF@icedyn_rhg_bbm', smoothCrossTF, cgt, 1._wp )
      !
   END FUNCTION smoothCrossTF





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



   SUBROUTINE div_stress_tensor( cgt, pe1e1, pe2e2,  pe1e1_e, pe2e2_e,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y,  &
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
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2, pe1e1_e, pe2e2_e
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2x, p1_e1e2y
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
      zt3(:,:) = ps12h(:,:) * pe1e1_e(:,:)
      zt4(:,:) = ps12h(:,:) * pe2e2_e(:,:)
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
      !                   !--- ds11/dx + ds12/dy
      pdivSx(ji,jj) = ( ( zt1(ji+ip,jj) - zt1(ji+im,jj) ) * p1_e2x(ji,jj) &
         &            + ( zt3(ji,jj+jp) - zt3(ji,jj+jm) ) * p1_e1x(ji,jj) &
         &                 ) * p1_e1e2x(ji,jj)
      !                   !--- ds22/dy + ds12/dx
      pdivSy(ji,jj) = ( ( zt2(ji,jj-jm) - zt2(ji,jj-jp) ) * p1_e1y(ji,jj) &
         &            + ( zt4(ji-im,jj) - zt4(ji-ip,jj) ) * p1_e2y(ji,jj) &
         &                 ) * p1_e1e2y(ji,jj)
      !
      END_2D
      !
   END SUBROUTINE div_stress_tensor

   SUBROUTINE div_stress_tensor_v2( cgt,  pe2x, pe1y, pe1e1, pe2e2,  pe1e1_e, pe2e2_e,  p1_e2x, p1_e1x, p1_e1y, p1_e2y, p1_e1e2x, p1_e1e2y,  &
      &                               ps11h, ps22h, ps12h,  pdivSx, pdivSy )
      !! => use the same discretization approach as what's being done in EVP, using `sigma_1` and `sigma_2` rather than `sigma_11` and `sigma_22`...
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
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe2x, pe1y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pe1e1, pe2e2, pe1e1_e, pe2e2_e
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e2x, p1_e1x, p1_e1y, p1_e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: p1_e1e2x, p1_e1e2y
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: ps11h, ps22h, ps12h ! components of stress tensors on T- or F-centric grids x h !!!
      REAL(wp), DIMENSION(:,:), INTENT(out) :: pdivSx, pdivSy      ! x,y components of the divergence of the tensor
      !!
      REAL(wp), DIMENSION(jpi,jpj) :: zs1, zs2, zs3
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
         CALL ctl_stop( 'STOP', 'div_stress_tensor_v2(): unknown grid-point type: '//cgt//'!')
      ENDIF
      !
      pdivSx(:,:) = 0._wp
      pdivSy(:,:) = 0._wp
      !
      zs1(:,:) = ps11h(:,:) + ps22h(:,:) ! h * sigma_1
      zs2(:,:) = ps11h(:,:) - ps22h(:,:) ! h * sigma_2
      zs3(:,:) = 2._wp * ps12h(:,:)      ! 2 * h * sigma_12
      !
      DO_2D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1 )
      !                   !--- U points
      pdivSx(ji,jj) = 0.5_wp * ( (( zs1(ji+ip,jj)                   - zs1(ji+im,jj)                  ) *   pe2x(ji,jj)  &
         &                      + ( zs2(ji+ip,jj)*pe2e2(ji+ip,jj)   - zs2(ji+im,jj)*pe2e2(ji+im,jj)  ) * p1_e2x(ji,jj)) &
         &                      + ( zs3(ji,jj+jp)*pe1e1_e(ji,jj+jp) - zs3(ji,jj+jm)*pe1e1_e(ji,jj+jm)) * p1_e1x(ji,jj)  &
         &                      ) * p1_e1e2x(ji,jj)
      !
      !                !--- V points
      pdivSy(ji,jj) = 0.5_wp * ( (( zs1(ji,jj-jm)                   - zs1(ji,jj-jp)                  ) *   pe1y(ji,jj)  &
         &                      - ( zs2(ji,jj-jm)*pe1e1(ji,jj-jm)   - zs2(ji,jj-jp)*pe1e1(ji,jj-jp)  ) * p1_e1y(ji,jj)) &
         &                      + ( zs3(ji-im,jj)*pe2e2_e(ji-im,jj) - zs3(ji-ip,jj)*pe2e2_e(ji-ip,jj)) * p1_e2y(ji,jj)  &
         &                      ) * p1_e1e2y(ji,jj)
      !
      END_2D
      !
   END SUBROUTINE div_stress_tensor_v2



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



   ELEMENTAL FUNCTION sigmaII( ps11, ps22, ps12 )
      !!------------------------------------------------------------------------------------
      !! Compute `sigma_II`: maximum shearing stress aka 2nd invariant of stress tensor [Pa]
      !!------------------------------------------------------------------------------------
      REAL(wp)             :: sigmaII
      REAL(wp), INTENT(in) :: ps11, ps22, ps12  ! sigma_11, sigma_22, sigma_12 [Pa]
      !!
      REAL(wp) :: ztmp, ztA, ztB
      !!------------------------------------------------------------------------------------
      ztmp  = 0.5_wp * (ps11 - ps22)
      ztA   = ztmp*ztmp
      ztB   = ps12*ps12
      !!
      sigmaII = SQRT( ztA + ztB )
      !!
   END FUNCTION sigmaII





#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_rhg_util
