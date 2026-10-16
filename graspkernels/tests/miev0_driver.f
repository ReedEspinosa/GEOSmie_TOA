C  Test driver for W. Wiscombe's MIEV0 (not part of GEOSmie; MIEV0.f and
C  ErrPack.f are taken from the user's copy, see test_independent_mie.py).
C  Build in double precision:
C    gfortran -O2 -fdefault-real-8 -o miev0_driver miev0_driver.f
C             MIEV0.f ErrPack.f
C  Input (stdin): NANG / angles [deg] / NX / NX lines "x mr mi"
C  Output (stdout), per size: x qext qsca, then for each angle
C    S11 S12 S33 S34 (Bohren & Huffman 1983 definitions)
      PROGRAM MIEV0DRV
      IMPLICIT NONE
      INTEGER MAXANG
      PARAMETER (MAXANG=501)
      INTEGER NANG, NX, I, J, IPOLZN, MOMDIM, NMOM
      REAL XMU(MAXANG), ANG(MAXANG), XX, MR, MI, MIMCUT, QEXT, QSCA
      REAL GQSC, PMOM(0:1,1), SPIKE, D2R, A1, A2
      COMPLEX CREFIN, SFORW, SBACK, S1(MAXANG), S2(MAXANG)
      COMPLEX TFORW(2), TBACK(2)
      LOGICAL PERFCT, ANYANG, PRNT(2)
      D2R = ACOS(-1.0)/180.0
      READ (*,*) NANG
      READ (*,*) (ANG(I),I=1,NANG)
      DO I=1,NANG
         XMU(I) = COS(ANG(I)*D2R)
      ENDDO
      PERFCT = .FALSE.
      MIMCUT = 1.0E-12
      ANYANG = .TRUE.
      NMOM = 0
      IPOLZN = 0
      MOMDIM = 1
      PRNT(1) = .FALSE.
      PRNT(2) = .FALSE.
      READ (*,*) NX
      DO J=1,NX
         READ (*,*) XX, MR, MI
         CREFIN = CMPLX(MR, -MI)
         CALL MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NANG, XMU,
     &        NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC, PMOM,
     &        SFORW, SBACK, S1, S2, TFORW, TBACK, SPIKE )
         WRITE (*,'(3ES25.16)') XX, QEXT, QSCA
         DO I=1,NANG
            A1 = ABS(S1(I))**2
            A2 = ABS(S2(I))**2
            WRITE (*,'(4ES25.16)') 0.5*(A1+A2), 0.5*(A2-A1),
     &           REAL(S2(I)*CONJG(S1(I))), AIMAG(S2(I)*CONJG(S1(I)))
         ENDDO
      ENDDO
      END
