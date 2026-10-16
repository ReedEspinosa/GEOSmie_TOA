C  Demo main program of M. I. Mishchenko's spher.f (version of 08/06/2005),
C  moved verbatim out of spher.f so that spher.f can be linked as a library.
C  Build:  gfortran -O2 -o spher_main spher_main.f spher.f
C  Running it must reproduce tests/spher_upstream_test_output.txt
C  (the reference output shipped at the top of the upstream file).

      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      PARAMETER (NMIE=10000, NPL=2*NMIE)
      REAL*8 AL1(NPL),AL2(NPL),AL3(NPL),AL4(NPL),BET1(NPL),           
     &       BET2(NPL),MRR,MRI,LAM                                            
      open (6, file='spher.print')
      open (10, file='spher.write')

      AA=0.6D0
      BB=0.2 D0
      AA1=0.17D0
      AA2=3.44D0
      BB1=DLOG(1.96D0)*DLOG(1.96D0)
      BB2=DLOG(2.37D0)*DLOG(2.37D0)
      GAM=1D0                                                               
      LAM=0.63D0                                            
      MRR=1.53 D0                                               
      MRI=0.008 D0                                            
      NDISTR=3
      NK=100
      N=100
      NP=4
      R1=0 D0
      R2=35D0       
      DDELT=1D-7
      WRITE (6,1000) DDELT
 1000 FORMAT ('NUMERICAL ACCURACY =',D8.2)
      CALL SPHER (AA,BB,GAM,LAM,MRR,MRI,R1,R2,N,NP,NDISTR,            
     *       NK,L1,AL1,AL2,AL3,AL4,BET1,BET2,CEXT,CSCA,AREA,VOL,RVW,
     *       RMEAN,AA1,BB1,AA2,BB2,DDELT)
      QE=CEXT/AREA
      LMAX=L1-1
      WRITE (6,1001) AREA,VOL,RVW,RMEAN
      CALL MATR (AL1,AL2,AL3,AL4,BET1,BET2,LMAX)
 1001 FORMAT ('<G> = ',D12.6,'  <V> = ',D12.6,'  Rvw = ',D12.6,
     *        '  <R> = ',D12.6)
      STOP                                                         
      END                                     
