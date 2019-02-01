module bhmie_wrapper

  use types
  use materials
  use distributions
  use bhmie_routine, pi => pii
  !$use omp_lib

  implicit none
  save

  real(dp),parameter :: refractive_index = 1.0_dp
  ! refractive index of medium, in this case vacuum

contains

  subroutine compute_dust_properties(prefix, output_format, abundance_mass, m, d, &
       & density, gas_to_dust, amin, amax, na, wavelengths, n_angles, n_small_angles,n_threads)

    implicit none

    character(len=*),intent(in) :: prefix
    ! dust model name

    integer,intent(in) :: output_format
    ! the file format for outputting

    integer,intent(in) :: n_threads

    type(material),intent(in) :: m(:)
    ! list of materials

    type(size_distribution),intent(in) :: d(:)
    ! list of size distributions

    real(dp),intent(in) :: amin, amax
    ! min and max of the dust size distribution

    integer,intent(in) :: na
    ! number of size bins to use

    real(dp),intent(in) :: abundance_mass(:), density(:), gas_to_dust, wavelengths(:)
    ! the density of the grains

    integer,intent(inout) :: n_angles, n_small_angles
    real(dp) :: angles(n_angles + n_small_angles)
    ! the number of angles to compute, and optionally the actual angles to use

    real(dp) :: a1, a2, a
    ! variables used to find the limits and center of size bins

    integer :: ic
    ! loop variable for components

    character(len=100) :: suffix

    ! running totals
    real(dp) :: cext(size(wavelengths))
    real(dp) :: csca(size(wavelengths))
    real(dp) :: cback(size(wavelengths))
    real(dp) :: kappa_ext(size(wavelengths))
    real(dp) :: gsca(size(wavelengths))
    real(dp) :: angles_full(2*(n_angles+n_small_angles)-1)
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
    complex(dp) :: s1_i(2*(n_angles+n_small_angles)-1)
    complex(dp) :: s2_i(2*(n_angles+n_small_angles)-1)
    real(dp) :: s11_i(2*(n_angles+n_small_angles)-1)
    real(dp) :: s12_i(2*(n_angles+n_small_angles)-1)
    real(dp) :: s33_i(2*(n_angles+n_small_angles)-1)
    real(dp) :: s34_i(2*(n_angles+n_small_angles)-1)
    real(dp) :: x
    ! size parameter

    real(dp) :: cross_section
    real(dp) :: volume
    real(dp) :: weight_number
    ! variables used in the summation

    real(dp),allocatable :: abundance_number(:)
    real(dp),allocatable :: average_particle_mass(:)
    ! number abundance and average particle mass

    integer :: iw,ia
    ! loop variables

    real(dp) :: logamin, logastep
    ! convenience variables for size distribution

    integer :: current

    ! Calculate average particle mass, and number abundance

    allocate(abundance_number(size(m)))
    allocate(average_particle_mass(size(m)))

    do ic=1,size(m)

       ! Average particle mass is average volume times density
       average_particle_mass(ic) = average_volume(d(ic)) / 1.e12 * density(ic)

       ! Number abundance
       abundance_number(ic) = abundance_mass(ic) / average_particle_mass(ic)

    end do

    ! Renormalize abundances
    abundance_number = abundance_number / sum(abundance_number)

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

    call print_progress_bar(1,na)

    !$call omp_set_num_threads(n_threads)

    !$omp parallel do default(firstprivate) shared(cext,csca,cback,gsca,kappa_ext,s11,s12,s33,s34) num_threads(n_threads)
    do ia=1,na
       !$if ia .eq. 1 print*,"Number of threads in use = ", omp_get_num_threads()

       call delete_progress_bar(ia,na)
       call print_progress_bar(ia,na)

       ! Find lower, central, and upper size for current bin
       a1 = 10._dp**(logamin + logastep * (real(ia, dp)-1._dp))
       a  = 10._dp**(logamin + logastep * (real(ia, dp)-0.5_dp))
       a2 = 10._dp**(logamin + logastep * (real(ia, dp)))

       ! Find the cross section and volume for the current grains
       cross_section = pi*a*a*1.e-8
       volume = 4._dp/3._dp*pi*a*a*a*1.e-12

       
       do ic=1,size(m)

          ! Find the weights to be applied to this size (in number)
          weight_number = distribution_weight_number(d(ic),a1,a2) * abundance_number(ic)

          if(weight_number > 0._dp) then

             ! Loop over wavelengths, and find properties for each
             do iw=1,size(m(ic)%wavelengths)

                ! Compute dimensionless size parameter
                x = 2._dp*pi*a/m(ic)%wavelengths(iw)

                ! Compute properties for size/wavelength
                call bhmie(x, m(ic)%refractive_indices(iw), size(angles), &
                     & s1_i, s2_i, qext_i, qsca_i, qback_i, gsca_i(iw), angles)

                ! Compute extinction, scattering, and backscattering cross-sections
                cext_i(iw) = qext_i * cross_section
                csca_i(iw) = qsca_i * cross_section
                cback_i(iw) = qback_i * cross_section

                ! Compute scattering matrix elements
                s11_i = 0.5_dp*(+ abs(s1_i)*abs(s1_i) + abs(s2_i)*abs(s2_i))
                s12_i = 0.5_dp*(- abs(s1_i)*abs(s1_i) + abs(s2_i)*abs(s2_i))
                s33_i = real(s2_i * conjg(s1_i), dp)
                s34_i = aimag(s2_i * conjg(s1_i))

                ! Add to running total
                !$omp critical (inner_loop)
                s11(iw, :) = s11(iw, :) + s11_i * weight_number
                s12(iw, :) = s12(iw, :) + s12_i * weight_number
                s33(iw, :) = s33(iw, :) + s33_i * weight_number
                s34(iw, :) = s34(iw, :) + s34_i * weight_number
                !$omp end critical (inner_loop)

             end do

             ! Add values for single size to the totals
             !$omp critical (middle_loop)
             cext = cext + cext_i * weight_number
             csca = csca + csca_i * weight_number
             cback = cback + cback_i * weight_number
             gsca = gsca + gsca_i * csca_i * weight_number
             !$omp end critical (middle_loop)

          end if

       end do

    end do
    !$omp end parallel do

    kappa_ext = cext * sum(abundance_mass / average_particle_mass)

    ! Normalize g
    gsca = gsca / csca

    ! Create array with all the angles
    do ia=1,size(angles)
       angles_full(ia) = angles(ia)
    end do
    do ia=1,size(angles)-1
       angles_full(size(angles) + ia) = pi - angles(size(angles)-ia)
    end do

    ! Correct opacity for gas-to-dust
    kappa_ext = kappa_ext / (1. + gas_to_dust)

    ! Output

    select case(output_format)
    case(1)

       open(unit=20,file=trim(prefix)//'.summary')
       do iw=1,size(wavelengths)
          write(20,'(6(ES11.4,2X))') wavelengths(iw), cext(iw), csca(iw), &
               & kappa_ext(iw), gsca(iw), -s12(iw,size(angles))/s11(iw,size(angles))
       end do
       close(unit=20)

    case(2)

       call write_1d_array(trim(prefix)//'.wav','ES23.16',wavelengths) ! Wavelengths
       call write_1d_array(trim(prefix)//'.mu','ES23.16',cos(angles_full)) ! Angles
       call write_1d_array(trim(prefix)//'.alb','ES11.4',csca/cext)! Albedo
       call write_1d_array(trim(prefix)//'.chi','ES11.4',kappa_ext) ! Chi (Kappa to Extinction)
       call write_1d_array(trim(prefix)//'.g','ES11.4',gsca) ! Average cos(theta)

       call write_2d_array(trim(prefix)//'.f11','ES11.4',s11)
       call write_2d_array(trim(prefix)//'.f12','ES11.4',s12)
       call write_2d_array(trim(prefix)//'.f33','ES11.4',s33)
       call write_2d_array(trim(prefix)//'.f34','ES11.4',s34)

    case(3)

       do iw=1,size(wavelengths)
          write(suffix,'(ES11.4)') wavelengths(iw)
          open(unit=20,file=trim(adjustl(prefix))//'.'//trim(adjustl(suffix)))
          write(20,'("Dust properties calculated using Bohren and Huffman subroutine")')
          write(20,'("kappa calculated using gas-to-dust ratio of ",F6.2)') gas_to_dust
          write(20,'(F8.4," = wavelength (microns)")') wavelengths(iw)
          write(20,'(ES11.4," = <cext>")') cext(iw)
          write(20,'(ES11.4," = <csca>")') csca(iw)
          write(20,'(ES11.4," = kappa (cm^2/g)")') kappa_ext(iw)*(1._dp - csca(iw)/cext(iw))
          write(20,'(F7.4," = <cos(theta)>")') gsca(iw)
          write(20,*)
          write(20,'(3X,"angle",5X,"S11",9X,"S22",9X,"S33",9X,"S44",9X,"S12",9X,"S34",9X)')
          do ia=1,size(angles)*2-1
             write(20,'(F7.2,6(1X,ES11.4))') angles_full(ia)*180._dp/pi, &
                  & s11(iw,ia),s11(iw,ia),s33(iw,ia),s33(iw,ia),s12(iw,ia),s34(iw,ia)
          end do
          write(20,*)
          close(unit=20)
       end do

    end select

  end subroutine compute_dust_properties

  subroutine print_progress_bar(i, imax)
    implicit none
    integer,intent(in) :: i, imax
    character(len=1), parameter :: bar = '='
    integer :: k
    write(6,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') 100*i/imax,'%','|', (bar, k =1,50*i/imax)
    close(6)
    open(6)
    if(i==imax) write(6,'(a)') '| done.'
  end subroutine print_progress_bar

  subroutine delete_progress_bar(i, imax)
    implicit none
    integer,intent(in) :: i, imax
    character(len=1), parameter :: back = char(8)
    integer :: k
    write(6,'(256a1)', advance='no') (back, k =1,(50*i/imax)+9)
  end subroutine delete_progress_bar

  subroutine write_1d_array(filename,fmt,array)
    implicit none
    character(len=*),intent(in) :: filename,fmt
    real(dp),intent(in) :: array(:)
    integer :: i
    open(unit=20,file=filename)
    do i=1,size(array)
       write(20,'('//trim(fmt)//')') array(i)
    end do
    close(unit=20)
  end subroutine write_1d_array

  subroutine write_2d_array(filename,fmt,array)
    implicit none
    character(len=*),intent(in) :: filename,fmt
    real(dp),intent(in) :: array(:,:)
    integer i1,i2
    open(unit=20,file=filename)
    do i1=1,size(array,1)
       do i2=1,size(array,2)
          write(20,'('//trim(fmt)//'," ")',advance='no') array(i1,i2)
       end do
       write(20,*)
    end do
    close(unit=20)
  end subroutine write_2d_array

end module bhmie_wrapper

