      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
      IMPLICIT NONE

C Declare parameters:
C Note: important that MXNANG be consistent with dimension of S1 and S2
C       in calling routine!

      INTEGER MXNANG,NMXX
C      PARAMETER(MXNANG=1000,NMXX=15000)
      PARAMETER(MXNANG=1000,NMXX=150000)

C Arguments:

      INTEGER NANG
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)

C Local variables:

      LOGICAL SINGLE
      INTEGER J,JJ,N,NSTOP,NMX,NN
      DOUBLE PRECISION CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,
     &                 THETA,XSTOP,YMOD
      DOUBLE PRECISION
     &   AMU(MXNANG),
     &   PI(MXNANG),
     &   PI0(MXNANG),
     &   PI1(MXNANG),
     &   TAU(MXNANG)
      DOUBLE COMPLEX
     &   DCXS1(2*MXNANG-1),
     &   DCXS2(2*MXNANG-1)

C***********************************************************************
C
C Subroutine BHMIE is derived from the Bohren-Huffman Mie scattering
C     subroutine to calculate scattering and absorption by a homogenous
C     isotropic sphere.
C Given:
C    X = 2*pi*a/lambda
C    REFREL = (complex refr. index of sphere)/(real index of medium)
C    NANG = number of angles between 0 and 90 degrees
C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
C           if called with NANG<2, will set NANG=2 and will compute
C           scattering for theta=0,90,180.
C Returns:
C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
C                                scatt. E perp. to scatt. plane)
C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
C                                scatt. E parr. to scatt. plane)
C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
C    QBACK = 4.*pi*(dC_sca/domega)/pi*a**2
C          = backscattering efficiency
C    GSCA = <cos(theta)> for scattering
C
C S1 and S2 are the diagonal elements of the "amplitude scattering matrix"
C (see eq. 3.12 of Bohren & Huffman 1983) -- the off-diagonal elements
C vanish for a spherical target.
C For unpolarized incident light, the intensity of scattered light a
C distance r from the sphere is just
C          1
C  I_s = ------ * I_in * S_11
C        (kr)^2
C
C where k=2*pi/lambda 
C and the "Muller matrix element" S_11 = 0.5*( |S_1|^2 + |S_2|^2 )
C
C for incident light polarized perp to the scattering plane,
C the scattered light is polarized perp to the scattering plane
C with intensity I_s = I_in * |S_1|^2 / (kr)^2
C
C for incident light polarized parallel to the scattering plane,
C the scattered light is polarized parallel to the scattering plane
C with intensity I_s = I_in * |S_2|^2 / (kr)^2
C
C History:
C Original program taken from Bohren and Huffman (1983), Appendix A
C Modified by B.T.Draine, Princeton Univ. Obs., 90.10.26
C in order to compute <cos(theta)>
C 91.05.07 (BTD): Modified to allow NANG=1
C 91.08.15 (BTD): Corrected error (failure to initialize P)
C 91.08.15 (BTD): Modified to enhance vectorizability.
C 91.08.15 (BTD): Modified to make NANG=2 if called with NANG=1
C 91.08.15 (BTD): Changed definition of QBACK.
C 92.01.08 (BTD): Converted to full double precision and double complex
C                 eliminated 2 unneed lines of code
C                 eliminated redundant variables (e.g. APSI,APSI0)
C                 renamed RN -> EN = double precision N
C                 Note that DOUBLE COMPLEX and DCMPLX are not part
C                 of f77 standard, so this version may not be fully
C                 portable.  In event that portable version is
C                 needed, use src/bhmie_f77.f
C 93.06.01 (BTD): Changed AMAX1 to generic function MAX
C 98.09.17 (BTD): Added variable "SINGLE" and warning in event that
C                 code is used with single-precision arithmetic (i.e.,
C                 compiler does not support DOUBLE COMPLEX)
C 99.02.17 (BTD): Replaced calls to REAL() and IMAG() by
C                 REALPART() and IMAGPART() for compatibility with g77
C                 Note that when code is used with standard f77 
C                 compilers, it is now necessary to enable two lines
C                 defining functions REALPART(X) and IMAGPART(X)
C 99.02.19 (BTD): added lines to be enabled to properly define
C                 REALPART() and IMAGPART() if NOT using g77
C                 ***see below!!***
C 01.02.16 (BTD): added IMPLICIT NONE
C 01.02.27 (BTD): changed definition of QBACK back to convention of
C                 Bohren & Huffman and others:
C                 Q_back = 4.*pi*(dC_sca/dOmega)/(pi*a^2) in backward
C                          direction
c 02.03.09 (BTD): defined statement function REALPART_SP to
c                 avoid warning regarding type conversion when taking
c                 real part of S1(1) to evaluate QEXT
c                 some cleanup regarding type conversion
c 02.05.30 (BTD): introduced internal double complex arrays DCXS1,DCXS2
c                 to possibly increase accuracy during summations.
c                 After summations, output scattering amplitudes
c                 via single complex arrays S1,S2 as before.
c                 Usage of this routine is unaffected by change.
c                 Note: no longer need statement function REALPART_SP
c 02.09.18 (BTD): Error in evaluation of QBACK reported by Okada Yasuhiko
c                 Was calculating QBACK using S1 rather than DCXS1
c                 Corrected.
c 02.10.16 (BTD): Added comments explaining definition of S_1 and S_2 .
C end history
C
C***********************************************************************
C 
C This module is dependent on whether compiler supports double precision
C complex variables:
C
C If your compiler does NOT support double complex, comment out following
C three lines, and uncomment corresponding 3 lines further below
C
      DOUBLE COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
      DOUBLE COMPLEX D(NMXX)
      PARAMETER(SINGLE=.FALSE.)

C      COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
C      COMPLEX D(NMXX)
C      PARAMETER(SINGLE=.TRUE.)

C**********************************************************************

C Following five statements should be enabled if NOT using g77.
C They assume that the compiler supports double complex, since the
C statements DBLE and DIMAG are used.  If double complex is not available
C (see above) you will need to change DIMAG to AIMAG
C
C If using g77, following statements could be commented out, as 
C REALPART and IMAGPART are g77 intrinsic functions
C However, they do not need to be commented out.

      DOUBLE COMPLEX DPCX
      DOUBLE PRECISION REALPART
      DOUBLE PRECISION IMAGPART
      REALPART(DPCX)=(DBLE(DPCX))
      IMAGPART(DPCX)=(DIMAG(DPCX))
      
C***********************************************************************
C*** Safety checks

      IF(SINGLE)WRITE(0,*)'Warning: this version of bhmie uses only ',
     &          'single precision complex numbers!'
      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2

C*** Obtain pi:

      PII=4.D0*ATAN(1.D0)
      DX=X
      DREFRL=REFREL
      Y=X*DREFRL
      YMOD=ABS(Y)

C*** Series expansion terminated after NSTOP terms
C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.*X**0.3333+2.
      NMX=NINT(MAX(XSTOP,YMOD))+15
C BTD experiment 91.1.15: add one more term to series and compare results
C      NMX=MAX(XSTOP,YMOD)+16
C test: compute 7001 wavelengths between .0001 and 1000 micron
C for a=1.0micron SiC grain.  When NMX increased by 1, only a single
C computed number changed (out of 4*7001) and it only changed by 1/8387
C conclusion: we are indeed retaining enough terms in series!

      NSTOP=NINT(XSTOP)

      IF(NMX.GT.NMXX)THEN
         WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
         STOP
      ENDIF

C*** Require NANG.GE.1 in order to calculate scattering intensities

      DANG=0.
      IF(NANG.GT.1)DANG=.5*PII/DBLE(NANG-1)
      DO J=1,NANG
         THETA=DBLE(J-1)*DANG
         AMU(J)=COS(THETA)
      ENDDO
      DO J=1,NANG
         PI0(J)=0.
         PI1(J)=1.
      ENDDO
      NN=2*NANG-1
      DO J=1,NN
         DCXS1(J)=(0.D0,0.D0)
         DCXS2(J)=(0.D0,0.D0)
      ENDDO

C*** Logarithmic derivative D(J) calculated by downward recurrence
C    beginning with initial value (0.,0.) at J=NMX

      D(NMX)=(0.,0.)
      NN=NMX-1
      DO N=1,NN
         EN=NMX-N+1
         D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
      ENDDO

C*** Riccati-Bessel functions with real argument X
C    calculated by upward recurrence

      PSI0=COS(DX)
      PSI1=SIN(DX)
      CHI0=-SIN(DX)
      CHI1=COS(DX)
      XI1=DCMPLX(PSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO N=1,NSTOP
         EN=N
         FN=(2.E0*EN+1.)/(EN*(EN+1.))

C for given N, PSI  = psi_n        CHI  = chi_n
C              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
C              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
C Calculate psi_n and chi_n

         PSI=(2.E0*EN-1.)*PSI1/DX-PSI0
         CHI=(2.E0*EN-1.)*CHI1/DX-CHI0
         XI=DCMPLX(PSI,-CHI)

C*** Store previous values of AN and BN for use
C    in computation of g=<cos(theta)>

         IF(N.GT.1)THEN
            AN1=AN
            BN1=BN
         ENDIF

C*** Compute AN and BN:

         AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
         AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
         BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
         BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)

C*** Augment sums for Qsca and g=<cos(theta)>

         QSCA=QSCA+REAL((2.*EN+1.)*(ABS(AN)**2+ABS(BN)**2))
         GSCA=GSCA+REAL(((2.*EN+1.)/(EN*(EN+1.)))*
     &        (REALPART(AN)*REALPART(BN)+IMAGPART(AN)*IMAGPART(BN)))
         IF(N.GT.1)THEN
            GSCA=GSCA+REAL(((EN-1.)*(EN+1.)/EN)*
     &      (REALPART(AN1)*REALPART(AN)+IMAGPART(AN1)*IMAGPART(AN)+
     &      REALPART(BN1)*REALPART(BN)+IMAGPART(BN1)*IMAGPART(BN)))
         ENDIF

C*** Now calculate scattering intensity pattern
C    First do angles from 0 to 90

         DO J=1,NANG
            JJ=2*NANG-J
            PI(J)=PI1(J)
            TAU(J)=EN*AMU(J)*PI(J)-(EN+1.)*PI0(J)
            DCXS1(J)=DCXS1(J)+FN*(AN*PI(J)+BN*TAU(J))
            DCXS2(J)=DCXS2(J)+FN*(AN*TAU(J)+BN*PI(J))
         ENDDO

C*** Now do angles greater than 90 using PI and TAU from
C    angles less than 90.
C    P=1 for N=1,3,...; P=-1 for N=2,4,...

         P=-P
         DO J=1,NANG-1
            JJ=2*NANG-J
            DCXS1(JJ)=DCXS1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
            DCXS2(JJ)=DCXS2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
         ENDDO
         PSI0=PSI1
         PSI1=PSI
         CHI0=CHI1
         CHI1=CHI
         XI1=DCMPLX(PSI1,-CHI1)

C*** Compute pi_n for next value of n
C    For each angle J, compute pi_n+1
C    from PI = pi_n , PI0 = pi_n-1

         DO J=1,NANG
            PI1(J)=((2.*EN+1.)*AMU(J)*PI(J)-(EN+1.)*PI0(J))/EN
            PI0(J)=PI(J)
         ENDDO
      ENDDO

C*** Have summed sufficient terms.
C    Now compute QSCA,QEXT,QBACK,and GSCA

      GSCA=REAL(2.D0*GSCA/QSCA)
      QSCA=REAL((2.D0/(DX*DX))*QSCA)
      QEXT=REAL((4.D0/(DX*DX))*REALPART(DCXS1(1)))
      QBACK=REAL(4.D0*(ABS(DCXS1(2*NANG-1))/DX)**2)

C prepare single precision complex scattering amplitude for output

      DO J=1,2*NANG-1
         S1(J)=CMPLX(DCXS1(J))
         S2(J)=CMPLX(DCXS2(J))
      ENDDO

      RETURN
      END
