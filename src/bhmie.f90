module bhmie_routine

  use types
  implicit none
  save

  real(dp),parameter :: pii = 3.14159265358979323846264338327950288419716939937510582097494_dp

contains

  subroutine bhmie(x,refrel,nang,s1,s2,qext,qsca,qback,gsca,angles)

    implicit none

    ! Declare parameters:

    integer,parameter :: nmxx=1000000

    ! Arguments:

    integer, intent(in) :: nang
    real(dp), intent(in) :: x
    complex(dp), intent(in) :: refrel
    complex(dp), dimension(:),intent(out) :: s1,s2
    real(dp), intent(out) :: gsca,qback,qext,qsca
    real(dp),intent(in),optional :: angles(:)

    ! Local variables:

    integer :: j,jj,n,nstop,nmx,nn
    real(dp) :: chi,chi0,chi1,dang,dx,en,fn,p,psi,psi0,psi1,theta,xstop,ymod
    real(dp), dimension(nang) :: amu, pi, pi0, pi1, tau
    complex(dp) :: an,an1,bn,bn1,drefrl,xi,xi1,y

    complex(dp),allocatable,dimension(:) :: d
    allocate(d(nmxx))

    !***********************************************************************
    !
    ! Subroutine BHMIE is derived from the Bohren-Huffman Mie scattering
    !     subroutine to calculate scattering and absorption by a homogenous
    !     isotropic sphere.
    ! Given:
    !    X = 2*pi*a/lambda
    !    REFREL = (complex refr. index of sphere)/(real index of medium)
    !    NANG = number of angles between 0 and 90 degrees
    !           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
    !           if called with NANG<2, will set NANG=2 and will compute
    !           scattering for theta=0,90,180.
    ! Returns:
    !    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
    !                                scatt. E perp. to scatt. plane)
    !    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
    !                                scatt. E parr. to scatt. plane)
    !    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
    !    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
    !    QBACK = 4.*pi*(dC_sca/domega)/pi*a**2
    !          = backscattering efficiency
    !    GSCA = <cos(theta)> for scattering
    !
    ! S1 and S2 are the diagonal elements of the "amplitude scattering matrix"
    ! (see eq. 3.12 of Bohren & Huffman 1983) -- the off-diagonal elements
    ! vanish for a spherical target.
    ! For unpolarized incident light, the intensity of scattered light a
    ! distance r from the sphere is just
    !          1
    !  I_s = ------ * I_in * S_11
    !        (kr)^2
    !
    ! where k=2*pi/lambda 
    ! and the "Muller matrix element" S_11 = 0.5*( |S_1|^2 + |S_2|^2 )
    !
    ! for incident light polarized perp to the scattering plane,
    ! the scattered light is polarized perp to the scattering plane
    ! with intensity I_s = I_in * |S_1|^2 / (kr)^2
    !
    ! for incident light polarized parallel to the scattering plane,
    ! the scattered light is polarized parallel to the scattering plane
    ! with intensity I_s = I_in * |S_2|^2 / (kr)^2
    !
    ! History:
    ! Original program taken from Bohren and Huffman (1983), Appendix A
    ! Modified by B.T.Draine, Princeton Univ. Obs., 90.10.26
    ! in order to compute <cos(theta)>
    ! 91.05.07 (BTD): Modified to allow NANG=1
    ! 91.08.15 (BTD): Corrected error (failure to initialize P)
    ! 91.08.15 (BTD): Modified to enhance vectorizability.
    ! 91.08.15 (BTD): Modified to make NANG=2 if called with NANG=1
    ! 91.08.15 (BTD): Changed definition of QBACK.
    ! 92.01.08 (BTD): Converted to full real(dp) and complex(dp)
    !                 eliminated 2 unneed lines of code
    !                 eliminated redundant variables (e.g. APSI,APSI0)
    !                 renamed RN -> EN = real(dp) N
    !                 Note that complex(dp) and DCMPLX are not part
    !                 of f77 standard, so this version may not be fully
    !                 portable.  In event that portable version is
    !                 needed, use src/bhmie_f77.f
    ! 93.06.01 (BTD): Changed AMAX1 to generic function MAX
    ! 98.09.17 (BTD): Added variable "SINGLE" and warning in event that
    !                 code is used with single-precision arithmetic (i.e.,
    !                 compiler does not support complex(dp))
    ! 99.02.17 (BTD): Replaced calls to REAL() and IMAG() by
    !                 REAL() and AIMAG() for compatibility with g77
    !                 Note that when code is used with standard f77 
    !                 compilers, it is now necessary to enable two lines
    !                 defining functions REAL(X) and AIMAG(X)
    ! 99.02.19 (BTD): added lines to be enabled to properly define
    !                 REAL() and AIMAG() if NOT using g77
    !                 ***see below!!***
    ! 01.02.16 (BTD): added IMPLICIT NONE
    ! 01.02.27 (BTD): changed definition of QBACK back to convention of
    !                 Bohren & Huffman and others:
    !                 Q_back = 4.*pi*(dC_sca/dOmega)/(pi*a^2) in backward
    !                          direction
    ! 02.03.09 (BTD): defined statement function REAL_SP to
    !                 avoid warning regarding type conversion when taking
    !                 real part of S1(1) to evaluate QEXT
    !                 some cleanup regarding type conversion
    ! 02.05.30 (BTD): introduced internal complex(dp) arrays S1,S2
    !                 to possibly increase accuracy during summations.
    !                 After summations, output scattering amplitudes
    !                 via single complex arrays S1,S2 as before.
    !                 Usage of this routine is unaffected by change.
    !                 Note: no longer need statement function REAL_SP
    ! 02.09.18 (BTD): Error in evaluation of QBACK reported by Okada Yasuhiko
    !                 Was calculating QBACK using S1 rather than S1
    !                 Corrected.
    ! 02.10.16 (BTD): Added comments explaining definition of S_1 and S_2 .
    ! 10.07.23 (TPR): Converted to Fortran 95, and inputs and outputs have
    !                 64-bit precision. Added option to specify directly
    !                 which angles to compute the scattering matrix for.
    ! end history
    !
    !***********************************************************************

    !***********************************************************************
    !*** Safety checks

    if(nang.lt.2) stop "nang should be > 1"

    !*** Obtain pi:

    dx=x
    drefrl=refrel
    y=x*drefrl
    ymod=abs(y)

    !*** Series expansion terminated after NSTOP terms
    !    Logarithmic derivatives calculated from NMX on down

    xstop=x+4._dp*x**0.3333_dp+2._dp
    nmx=nint(max(xstop,ymod))+15
    nstop=nint(xstop)

    ! BTD experiment 91.1.15: add one more term to series and compare results
    !      NMX=MAX(XSTOP,YMOD)+16
    ! test: compute 7001 wavelengths between .0001 and 1000 micron
    ! for a=1.0micron SiC grain.  When NMX increased by 1, only a single
    ! computed number changed (out of 4*7001) and it only changed by 1/8387
    ! conclusion: we are indeed retaining enough terms in series!

    if(nmx.gt.nmxx)then
       write(0,*)'error: nmx > nmxx=',nmxx,' for |m|x=',ymod
       stop
    endif

    !*** Require NANG.GE.1 in order to calculate scattering intensities

    if(present(angles)) then

       if(size(angles).ne.nang) stop "The number of angles specified should be equal to nang"
       if(angles(1).ne.0.) stop "angles(1) should be 0."
       if(angles(nang).ne.pii/2.) stop "angles(nang) should be pi/2."

       amu = cos(angles)

    else

       dang=0.

       if(nang.gt.1) dang=.5_dp*pii/real(nang-1,dp)

       do j=1,nang
          theta=real(j-1,dp)*dang
          amu(j)=cos(theta)
       end do

    end if

    pi0=0._dp
    pi1=1._dp

    s1=(0._dp,0._dp)
    s2=(0._dp,0._dp)

    !*** Logarithmic derivative D(J) calculated by downward recurrence
    !    beginning with initial value (0.,0.) at J=NMX

    d(nmx)=(0._dp,0._dp)
    nn=nmx-1

    do n=1,nn
       en=nmx-n+1
       d(nmx-n)=(en/y)-(1._dp/(d(nmx-n+1)+en/y))
    end do

    !*** Riccati-Bessel functions with real argument X
    !    calculated by upward recurrence

    psi0=cos(dx)
    psi1=sin(dx)
    chi0=-sin(dx)
    chi1=cos(dx)
    xi1=cmplx(psi1,-chi1)
    qsca=0._dp
    gsca=0._dp
    p=-1._dp
    do n=1,nstop
       en=n
       fn=(2._dp*en+1._dp)/(en*(en+1._dp))

       ! for given N, PSI  = psi_n        CHI  = chi_n
       !              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
       !              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
       ! Calculate psi_n and chi_n

       psi=(2._dp*en-1._dp)*psi1/dx-psi0
       chi=(2._dp*en-1._dp)*chi1/dx-chi0
       xi=cmplx(psi,-chi)

       !*** Store previous values of AN and BN for use
       !    in computation of g=<cos(theta)>

       if(n.gt.1)then
          an1=an
          bn1=bn
       endif

       !*** Compute AN and BN:

       an=(d(n)/drefrl+en/dx)*psi-psi1
       an=an/((d(n)/drefrl+en/dx)*xi-xi1)
       bn=(drefrl*d(n)+en/dx)*psi-psi1
       bn=bn/((drefrl*d(n)+en/dx)*xi-xi1)

       !*** Augment sums for Qsca and g=<cos(theta)>

       qsca=qsca+(2._dp*en+1._dp)*(abs(an)**2+abs(bn)**2)
       gsca=gsca+(2._dp*en+1._dp)/(en*(en+1._dp))*(real(an)*real(bn)+aimag(an)*aimag(bn))
       if(n.gt.1)then
          gsca=gsca+(en-1._dp)*(en+1._dp)/en* &
               &      (real(an1)*real(an)+aimag(an1)*aimag(an)+ &
               &      real(bn1)*real(bn)+aimag(bn1)*aimag(bn))
       endif

       !*** Now calculate scattering intensity pattern
       !    First do angles from 0 to 90

       do j=1,nang
          jj=2*nang-j
          pi(j)=pi1(j)
          tau(j)=en*amu(j)*pi(j)-(en+1.)*pi0(j)
          s1(j)=s1(j)+fn*(an*pi(j)+bn*tau(j))
          s2(j)=s2(j)+fn*(an*tau(j)+bn*pi(j))
       end do

       !*** Now do angles greater than 90 using PI and TAU from
       !    angles less than 90.
       !    P=1 for N=1,3,...; P=-1 for N=2,4,...

       p=-p
       do j=1,nang-1
          jj=2*nang-j
          s1(jj)=s1(jj)+fn*p*(an*pi(j)-bn*tau(j))
          s2(jj)=s2(jj)+fn*p*(bn*pi(j)-an*tau(j))
       end do
       psi0=psi1
       psi1=psi
       chi0=chi1
       chi1=chi
       xi1=cmplx(psi1,-chi1)

       !*** Compute pi_n for next value of n
       !    For each angle J, compute pi_n+1
       !    from PI = pi_n , PI0 = pi_n-1

       do j=1,nang
          pi1(j)=((2._dp*en+1._dp)*amu(j)*pi(j)-(en+1._dp)*pi0(j))/en
          pi0(j)=pi(j)
       end do
    end do

    !*** Have summed sufficient terms.
    !    Now compute QSCA,QEXT,QBACK,and GSCA

    gsca=2._dp*gsca/qsca
    qsca=(2._dp/(dx*dx))*qsca
    qext=(4._dp/(dx*dx))*real(s1(1))
    qback=4._dp*(abs(s1(2*nang-1))/dx)**2

  end subroutine bhmie

end module bhmie_routine
