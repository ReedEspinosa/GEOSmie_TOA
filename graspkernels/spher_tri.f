C**********************************************************************
C  SPHER_TRI: GRASP-type sphere kernels from Mishchenko's SPHER.
C
C  For every refractive index and every radius node r_k of a grid, it
C  computes the triangle-basis (hat function in ln r) average
C
C     K_k = INT T_k(u) * (C/V)(u) du ,  u = ln r ,
C
C  where T_k rises linearly from 0 at u_lo to 1 at u_k and falls back to
C  0 at u_hi, and C is the extinction or scattering cross section (or
C  C_sca times a normalized scattering-matrix element). K_k is the GRASP
C  text-kernel value, in 1/um when radii and wavelength are in um.
C
C  Each half of the hat is one call of SPHER, with the triangle size
C  distribution NDISTR=7 (the GEOSmie addition in spher.f). Because the
C  volume-weighted triangle in ln r is used as SPHER's size
C  distribution, <C>/<V> from SPHER is the T-weighted mean of C/V over
C  that half, and
C     K_k = SUM_halves (h/2) * <C>/<V> ,  h = width of the half in ln r.
C  Splitting at the apex keeps the kink out of the Gaussian quadrature.
C
C  Check mode (MODE=6): each half is instead computed with the
C  unmodified NDISTR=6 (bimodal volume log normal) with GAM=-1 and very
C  broad, offset modes, whose difference is linear in ln r to ~1e-9.
C  This exercises only upstream code and verifies the NDISTR=7 branch.
C
C  Usage:  spher_tri <job file> <output file>
C
C  Job file (list directed):
C     LAM  NK  DDELT  MODE         wavelength [um], Gauss points per
C                                  subinterval, SPHER accuracy, 7 or 6
C     NANG                         number of scattering angles
C     ANG(1:NANG)                  scattering angles [deg]
C     NM                           number of refractive indices
C     MRR MRI   (NM lines)         MRI >= 0
C     NR                           number of radius nodes
C     RK RLO RHI NLO NHI (NR lines)  apex and edges of the hat [um] and
C                                  number of quadrature subintervals in
C                                  each half (0 = omit that half)
C
C  Output (unformatted stream): integers NM, NR, NANG, then for each
C  refractive index (outer) and node (inner) the real*8 values
C     K_ext, K_sca, K_11(NANG), K_12, K_22, K_33, K_34, K_44
C  where K_ij = SUM_halves (h/2) * <C_sca>/<V> * F_ij and F_ij are
C  Mishchenko's normalized scattering-matrix elements (F11 = a1,
C  F22 = a2, F33 = a3, F44 = a4, F12 = b1, F34 = b2).
C
C  SPHER's own printout (unit 6) and coefficient file (unit 10) are sent
C  to /dev/null; progress goes to standard error.
C**********************************************************************
      PROGRAM SPHER_TRI
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMIE=10000, NPL=2*NMIE, NANGMX=2000, NRMX=2000)
      REAL*8 AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BET1(NPL),BET2(NPL)
      REAL*8 LAM,MRR,MRI
      REAL*8 ANG(NANGMX),RK(NRMX),RLO(NRMX),RHI(NRMX)
      INTEGER NLO(NRMX),NHI(NRMX)
      REAL*8 F(NANGMX,6),OUT(2+6*NANGMX)
      REAL*8, ALLOCATABLE :: RMR(:),RMI(:)
      CHARACTER*1000 JOBFN, OUTFN

      IF (IARGC().NE.2) THEN
         WRITE (0,*) 'Usage: spher_tri <job file> <output file>'
         STOP 1
      ENDIF
      CALL GETARG(1,JOBFN)
      CALL GETARG(2,OUTFN)

      OPEN (11,FILE=TRIM(JOBFN),STATUS='OLD')
      READ (11,*) LAM,NK,DDELT,MODE
      READ (11,*) NANG
      IF (NANG.GT.NANGMX) STOP 'NANG > NANGMX'
      READ (11,*) (ANG(I),I=1,NANG)
      READ (11,*) NM
      ALLOCATE (RMR(NM),RMI(NM))
      DO I=1,NM
         READ (11,*) RMR(I),RMI(I)
      ENDDO
      READ (11,*) NR
      IF (NR.GT.NRMX) STOP 'NR > NRMX'
      DO I=1,NR
         READ (11,*) RK(I),RLO(I),RHI(I),NLO(I),NHI(I)
      ENDDO
      CLOSE (11)
      IF (MODE.NE.7 .AND. MODE.NE.6) STOP 'MODE must be 7 or 6'

      OPEN (6,FILE='/dev/null')
      OPEN (10,FILE='/dev/null')
      OPEN (20,FILE=TRIM(OUTFN),ACCESS='STREAM',FORM='UNFORMATTED',
     &      STATUS='REPLACE')
      WRITE (20) NM,NR,NANG
      NOUT=2+6*NANG

      DO IM=1,NM
         MRR=RMR(IM)
         MRI=RMI(IM)
         DO IR=1,NR
            DO J=1,NOUT
               OUT(J)=0D0
            ENDDO
            DO IH=1,2
               IF (IH.EQ.1) THEN
                  R1=RLO(IR)
                  R2=RK(IR)
                  N=NLO(IR)
               ELSE
                  R1=RK(IR)
                  R2=RHI(IR)
                  N=NHI(IR)
               ENDIF
               IF (N.LE.0 .OR. R2.LE.R1) CYCLE
               H=DLOG(R2/R1)
               CALL HALF (MODE,IH,RK(IR),R1,R2,H,LAM,MRR,MRI,N,NK,
     &                    DDELT,AL1,AL2,AL3,AL4,BET1,BET2,L1,
     &                    CEXT,CSCA,VOLUME)
               CALL MATR_ANG (AL1,AL2,AL3,AL4,BET1,BET2,L1-1,
     &                        ANG,NANG,NANGMX,F)
               W=0.5D0*H/VOLUME
               OUT(1)=OUT(1)+W*CEXT
               OUT(2)=OUT(2)+W*CSCA
               DO K=1,6
                  DO J=1,NANG
                     OUT(2+(K-1)*NANG+J)=OUT(2+(K-1)*NANG+J)
     &                                   +W*CSCA*F(J,K)
                  ENDDO
               ENDDO
            ENDDO
            WRITE (20) (OUT(J),J=1,NOUT)
         ENDDO
         WRITE (0,'(A,I5,A,I5,A,F9.5,A,ES11.4)') ' spher_tri: m ',IM,
     &         ' of ',NM,'  done (mr=',MRR,', mi=',MRI
      ENDDO
      CLOSE (20)
      END

C**********************************************************************
C  One half of the hat on [R1,R2] (apex at RKK), via SPHER.
C  MODE=7: NDISTR=7 (triangle, GEOSmie addition to spher.f).
C  MODE=6: NDISTR=6 with GAM=-1 and BB1=BB2=BIG, centres placed so that
C          the difference of the two near-flat modes is proportional to
C          (u-u_lo) on the rising half (IH=1) or (u_hi-u) on the falling
C          half (IH=2).
C**********************************************************************
      SUBROUTINE HALF (MODE,IH,RKK,R1,R2,H,LAM,MRR,MRI,N,NK,DDELT,
     &                 AL1,AL2,AL3,AL4,BET1,BET2,L1,CEXT,CSCA,VOLUME)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMIE=10000, NPL=2*NMIE)
      REAL*8 AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BET1(NPL),BET2(NPL)
      REAL*8 LAM,MRR,MRI
      PARAMETER (BIG=1D6, D=1D0)
      NP=1
      IF (MODE.EQ.7) THEN
         NDISTR=7
         AA=RKK
         BB=H
         GAM=0D0
         AA1=0D0
         BB1=0D0
         AA2=0D0
         BB2=0D0
      ELSE
         NDISTR=6
         AA=0D0
         BB=0D0
         GAM=-1D0
         BB1=BIG
         BB2=BIG
         IF (IH.EQ.1) THEN
            AA1=R1*DEXP(D)
            AA2=R1*DEXP(-D)
         ELSE
            AA1=R2*DEXP(-D)
            AA2=R2*DEXP(D)
         ENDIF
      ENDIF
      CALL SPHER (AA,BB,GAM,LAM,MRR,MRI,R1,R2,N,NP,NDISTR,
     &            NK,L1,AL1,AL2,AL3,AL4,BET1,BET2,CEXT,CSCA,AREA,
     &            VOLUME,RVW,RMEAN,AA1,BB1,AA2,BB2,DDELT)
      RETURN
      END

C**********************************************************************
C  Scattering matrix at arbitrary angles from SPHER's expansion
C  coefficients. Copied from subroutine MATR in GEOSmie's
C  gsf/spher_expan.f (itself adapted from Mishchenko's MATR), renamed,
C  with array sizes and outputs passed as arguments instead of
C  params.h and COMMON /MATR_OUT/.
C     F(:,1..6) = F11, F12, F22, F33, F34, F44 (this order)
C**********************************************************************
      SUBROUTINE MATR_ANG (A1,A2,A3,A4,B1,B2,LMAX,ANGD,NANG,NANGMX,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NMIE=10000, NPL=2*NMIE)
      REAL*8 A1(NPL),A2(NPL),A3(NPL),A4(NPL),B1(NPL),B2(NPL)
      REAL*8 ANGD(NANGMX),F(NANGMX,6)
      D2R=DACOS(-1D0)/180D0
      L1MAX=LMAX+1
      D6=DSQRT(6D0)*0.25D0
      DO 500 I1=1,NANG
         TAA=ANGD(I1)*D2R
         U=DCOS(TAA)
         F11=0D0
         F2=0D0
         F3=0D0
         F44=0D0
         F12=0D0
         F34=0D0
         P1=0D0
         P2=0D0
         P3=0D0
         P4=0D0
         PP1=1D0
         PP2=0.25D0*(1D0+U)*(1D0+U)
         PP3=0.25D0*(1D0-U)*(1D0-U)
         PP4=D6*(U*U-1D0)
         DO 400 L1=1,L1MAX
            L=L1-1
            DL=DFLOAT(L)
            DL1=DFLOAT(L1)
            F11=F11+A1(L1)*PP1
            F44=F44+A4(L1)*PP1
            IF(L.EQ.LMAX) GO TO 350
            PL1=DFLOAT(2*L+1)
            P=(PL1*U*PP1-DL*P1)/DL1
            P1=PP1
            PP1=P
  350       IF(L.LT.2) GO TO 400
            F2=F2+(A2(L1)+A3(L1))*PP2
            F3=F3+(A2(L1)-A3(L1))*PP3
            F12=F12+B1(L1)*PP4
            F34=F34+B2(L1)*PP4
            IF(L.EQ.LMAX) GO TO 400
            PL2=DL*DL1*U
            PL3=DL1*(DL*DL-4D0)
            PL4=1D0/(DL*(DL1*DL1-4D0))
            P=(PL1*(PL2-4D0)*PP2-PL3*P2)*PL4
            P2=PP2
            PP2=P
            P=(PL1*(PL2+4D0)*PP3-PL3*P3)*PL4
            P3=PP3
            PP3=P
            P=(PL1*U*PP4-DSQRT(DL*DL-4D0)*P4)/DSQRT(DL1*DL1-4D0)
            P4=PP4
            PP4=P
  400    CONTINUE
         F(I1,1)=F11
         F(I1,2)=F12
         F(I1,3)=(F2+F3)*0.5D0
         F(I1,4)=(F2-F3)*0.5D0
         F(I1,5)=F34
         F(I1,6)=F44
  500 CONTINUE
      RETURN
      END
