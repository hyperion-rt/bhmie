module bhmie_wrapper

  use types
  use materials
  use distributions
  use bhmie_routine, pi => pii

  implicit none
  save

  real(dp),parameter :: refractive_index = 1.0_dp
  ! refractive index of medium, in this case vacuum

contains

  subroutine compute_dust_properties(prefix, abundance, m, d, density, amin, amax, na, wavelengths, n_angles, n_small_angles)

    implicit none

    character(len=*),intent(in) :: prefix
    ! dust model name

    type(material),intent(in) :: m(:)
    ! list of materials

    type(size_distribution),intent(in) :: d(:)
    ! list of size distributions

    real(dp),intent(in) :: amin, amax
    ! min and max of the dust size distribution

    integer,intent(in) :: na
    ! number of size bins to use

    real(dp),intent(in) :: abundance(:), density(:), wavelengths(:)
    ! the density of the grains

    integer,intent(inout) :: n_angles, n_small_angles
    real(dp) :: angles(n_angles + n_small_angles)
    ! the number of angles to compute, and optionally the actual angles to use

    real(dp) :: a1, a2, a
    ! variables used to find the limits and center of size bins

    integer :: ic
    ! loop variable for components

    ! running totals
    real(dp) :: cext(size(wavelengths))
    real(dp) :: csca(size(wavelengths))
    real(dp) :: cback(size(wavelengths))
    real(dp) :: kappa_ext(size(wavelengths))
    real(dp) :: gsca(size(wavelengths))
    real(dp) :: s11(size(wavelengths), 2*(n_angles+n_small_angles)-1)
    real(dp) :: s12(size(wavelengths), 2*(n_angles+n_small_angles)-1)
    real(dp) :: s33(size(wavelengths), 2*(n_angles+n_small_angles)-1)
    real(dp) :: s34(size(wavelengths), 2*(n_angles+n_small_angles)-1)

    ! individual
    real(dp) :: qext_i
    real(dp) :: qsca_i
    real(dp) :: qback_i
    real(dp) :: cext_i(size(wavelengths))
    real(dp) :: csca_i(size(wavelengths))
    real(dp) :: cback_i(size(wavelengths))
    real(dp) :: gsca_i(size(wavelengths))
    complex(dp) :: s1_i(size(wavelengths), 2*(n_angles+n_small_angles)-1)
    complex(dp) :: s2_i(size(wavelengths), 2*(n_angles+n_small_angles)-1)

    real(dp) :: x
    ! size parameter

    real(dp) :: cross_section
    real(dp) :: volume
    real(dp) :: weight_number
    real(dp) :: weight_mass
    ! variables used in the summation

    integer :: iw,ia
    ! loop variables

    real(dp) :: logamin, logastep
    ! convenience variables for size distribution

    integer :: current

    ! Compute angles
    angles(1) = 0.
    do ia=2,n_angles
       angles(n_small_angles + ia) = real(ia - 1) / real(n_angles - 1) * pi / 2._dp
    end do

    ! Add small angles at lower end (will also be mirrored at high end)
    do ia=1,n_small_angles
       current = n_small_angles + 2 - ia
       angles(current) = angles(current+1) / sqrt(10.) 
    end do

    ! Initialize running totals
    s11 = 0._dp
    s12 = 0._dp
    s33 = 0._dp
    s34 = 0._dp
    cext = 0._dp
    csca = 0._dp
    cback = 0._dp
    gsca = 0._dp
    kappa_ext = 0._dp

    ! Initialize convenience variables for size distribution
    logamin = log10(amin)
    logastep = (log10(amax) - log10(amin)) / real(na, dp)

    do ia=1,na

       ! Find lower, central, and upper size for current bin
       a1 = 10._dp**(logamin + logastep * (real(ia, dp)-1._dp))
       a  = 10._dp**(logamin + logastep * (real(ia, dp)-0.5_dp))
       a2 = 10._dp**(logamin + logastep * (real(ia, dp)))

       write(*,'("Computing ",ES11.4," to ",ES11.4)') a1, a2

       ! Find the cross section and volume for the current grains
       cross_section = pi*a*a*1.e-8
       volume = 4._dp/3._dp*pi*a*a*a*1.e-12

       do ic=1,size(m)

          ! Find the weights to be applied to this size, both in number and mass
          weight_number = distribution_weight_number(d(ic),a1,a2) * abundance(ic)
          weight_mass = distribution_weight_mass(d(ic),a1,a2) * abundance(ic)

          if(weight_number > 0._dp .or. weight_mass > 0._dp) then

             ! Loop over wavelengths, and find properties for each
             do iw=1,size(m(ic)%wavelengths)

                ! Compute dimensionless size parameter
                x = 2._dp*pi*a/m(ic)%wavelengths(iw)

                ! Compute properties for size/wavelength
                call bhmie(x, m(ic)%refractive_indices(iw), size(angles), &
                     & s1_i(iw, :), s2_i(iw, :), qext_i, qsca_i, qback_i, gsca_i(iw), angles)

                ! Compute extinction, scattering, and backscattering cross-sections
                cext_i(iw) = qext_i * cross_section
                csca_i(iw) = qsca_i * cross_section
                cback_i(iw) = qback_i * cross_section

             end do

             ! Compute scattering matrix elements
             s11 = s11 + 0.5_dp*(+ abs(s1_i)*abs(s1_i) + abs(s2_i)*abs(s2_i)) * weight_number
             s12 = s12 + 0.5_dp*(- abs(s1_i)*abs(s1_i) + abs(s2_i)*abs(s2_i)) * weight_number
             s33 = s33 + real(s2_i * conjg(s1_i), dp) * weight_number
             s34 = s34 + aimag(s2_i * conjg(s1_i)) * weight_number

             ! Add values for single size to the totals
             cext = cext + cext_i * weight_number
             csca = csca + csca_i * weight_number
             cback = cback + cback_i * weight_number
             gsca = gsca + gsca_i * weight_number

             ! Compute and add opacities to the totals
             kappa_ext = kappa_ext + cext_i / (volume * density(ic)) * weight_mass

          end if

       end do

    end do

    ! Output

    ! Wavelengths
    open(unit=20,file=trim(prefix)//'.wav')
    do iw=1,size(wavelengths)
       write(20,'(ES23.16)') wavelengths(iw)
    end do
    close(unit=20)

    ! Angles
    open(unit=20,file=trim(prefix)//'.mu')
    do ia=1,size(angles)
       write(20,'(ES23.16)') cos(angles(ia))
    end do
    do ia=1,size(angles)-1
       write(20,'(ES23.16)') -cos(angles(size(angles)-ia))
    end do
    close(unit=20)

    ! Albedo
    open(unit=20,file=trim(prefix)//'.alb')
    do iw=1,size(wavelengths)
       write(20,'(ES11.4)') csca(iw)/cext(iw)
    end do
    close(unit=20)

    ! Chi (Kappa to Extinction)
    open(unit=20,file=trim(prefix)//'.chi')
    do iw=1,size(wavelengths)
       write(20,'(ES11.4)') kappa_ext(iw)
    end do
    close(unit=20)

    ! g
    open(unit=20,file=trim(prefix)//'.g')
    do iw=1,size(wavelengths)
       write(20,'(ES11.4)') gsca(iw)
    end do
    close(unit=20)

    ! P11
    open(unit=20,file=trim(prefix)//'.f11')
    do iw=1,size(wavelengths)
       do ia=1,size(angles)*2-1
          write(20,'(ES11.4," ")',advance='no') s11(iw,ia)
       end do
       write(20,*)
    end do
    close(unit=20)

    ! P12
    open(unit=20,file=trim(prefix)//'.f12')
    do iw=1,size(wavelengths)
       do ia=1,size(angles)*2-1
          write(20,'(ES11.4," ")',advance='no') s12(iw,ia)
       end do
       write(20,*)
    end do
    close(unit=20)

    ! P33
    open(unit=20,file=trim(prefix)//'.f33')
    do iw=1,size(wavelengths)
       do ia=1,size(angles)*2-1
          write(20,'(ES11.4," ")',advance='no') s33(iw,ia)
       end do
       write(20,*)
    end do
    close(unit=20)

    ! P34
    open(unit=20,file=trim(prefix)//'.f34')
    do iw=1,size(wavelengths)
       do ia=1,size(angles)*2-1
          write(20,'(ES11.4," ")',advance='no') s34(iw,ia)
       end do
       write(20,*)
    end do
    close(unit=20)

  end subroutine compute_dust_properties

end module bhmie_wrapper

