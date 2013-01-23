module distributions

  use lib_array
  use types
  implicit none
  save

  type size_distribution
     integer :: type
     ! type = 1: power-law
     ! type = 2: numerical
     real(dp) :: amin, amax, apower, aturn
     real(dp),allocatable :: a(:), n(:), m(:)
  end type size_distribution

contains

  subroutine read_distribution(unit, d)

    implicit none

    integer,intent(in) :: unit
    type(size_distribution),intent(out) :: d
    character(len=100) :: distribution_type, filename
    real(dp) :: norm
    integer :: i, ioerr, na

    ! Read in the distribution type
    read(unit,*) distribution_type

    select case(trim(distribution_type))
    case('power')

       d%type = 1

       ! Read in min, max, and power of power-law distribution
       read(unit,*) d%amin, d%amax, d%apower

    case('ped')

       d%type = 2

       ! Read in min, turnoff, and power of PED distribution
       read(unit,*) d%amin, d%aturn, d%apower

       ! Define number of a values to compute the PED for
       na = 1000

       ! Allocate arrays
       allocate(d%a(na), d%n(na), d%m(na))

       do i=1,na
          d%a(i) = 10._dp ** (-5._dp + real(i, dp) / real(na, dp) * 10._dp)
          d%n(i) = d%a(i) ** d%apower * exp(-d%a(i)/d%aturn)
       end do

       ! Compute mass distribution
       d%m = d%n * d%a**3

       ! Set max
       d%amax = d%a(na)

       ! Normalize number distribution
       norm = integral_loglog(d%a, d%n)
       d%n = d%n / norm

       ! Normalize mass distribution
       norm = integral_loglog(d%a, d%m)
       d%m = d%m / norm

    case('table')

       d%type = 2

       ! Read in name of file containing the table
       read(unit,*) filename

       ! Find number of lines in file
       open(unit=35,file=filename,status='old')
       na = 0
       do
          read(35,*,iostat=ioerr)
          if(ioerr.ne.0) exit
          na = na + 1
       end do
       close(unit=35)

       ! Allocate arrays
       allocate(d%a(na), d%n(na), d%m(na))

       ! Read in distribution
       open(unit=35,file=filename,status='old')
       do i=1,na
          read(35,*,iostat=ioerr) d%a(i), d%n(i)
       end do
       close(unit=35)

       ! Compute mass distribution
       d%m = d%n * d%a**3

       ! Set min and max
       d%amin = d%a(1)
       d%amax = d%a(na)

       ! Normalize number distribution
       norm = integral_loglog(d%a, d%n)
       d%n = d%n / norm

       ! Normalize mass distribution
       norm = integral_loglog(d%a, d%m)
       d%m = d%m / norm

    case default
       stop "Not implemented"
    end select

  end subroutine read_distribution

  real(dp) function distribution_weight_number(d, amin, amax)

    implicit none

    type(size_distribution),intent(in) :: d
    real(dp),intent(in) :: amin, amax
    real(dp) :: a1, a2

    if(amax < d%amin .or. amin > d%amax) then
       distribution_weight_number = 0._dp
    else if(amin < d%amin .and. amax > d%amax) then
       distribution_weight_number = 1._dp
    else
       a1 = max(amin, d%amin)
       a2 = min(amax, d%amax)
       select case(d%type)
       case(1)
          distribution_weight_number = (a2**(d%apower+1._dp) - a1**(d%apower+1._dp)) &
               &                     / (d%amax**(d%apower+1._dp) - d%amin**(d%apower+1._dp))
       case(2)
          distribution_weight_number = integral_loglog(d%a, d%n, a1, a2)
       case default
          stop "Not implemented"
       end select
    end if

  end function distribution_weight_number

  real(dp) function average_volume(d)

    implicit none

    type(size_distribution),intent(in) :: d

       select case(d%type)
       case(1)
           average_volume = 4. / 3. * 3.1415926 * (d%amax**(d%apower+4._dp) - d%amin**(d%apower+4._dp)) &
                                                / (d%amax**(d%apower+1._dp) - d%amin**(d%apower+1._dp)) &
                                                * (d%apower + 1._dp) / (d%apower + 4._dp)
       case(2)
           average_volume = integral_loglog(d%a, d%n * 4. / 3. * 3.1415926 * d%a ** 3) / integral_loglog(d%a, d%n)
       case default
          stop "Not implemented"
       end select

  end function average_volume

  real(dp) function distribution_weight_mass(d, amin, amax)

    implicit none

    type(size_distribution),intent(in) :: d
    real(dp),intent(in) :: amin, amax
    real(dp) :: a1, a2

    if(amax < d%amin .or. amin > d%amax) then
       distribution_weight_mass = 0._dp
    else if(amin < d%amin .and. amax > d%amax) then
       distribution_weight_mass = 1._dp
    else
       a1 = max(amin, d%amin)
       a2 = min(amax, d%amax)
       select case(d%type)
       case(1)
          distribution_weight_mass = (a2**(d%apower+4._dp) - a1**(d%apower+4._dp)) &
               &                   / (d%amax**(d%apower+4._dp) - d%amin**(d%apower+4._dp))
       case(2)
          distribution_weight_mass = integral_loglog(d%a, d%m, a1, a2)
       case default
          stop "Not implemented"
       end select

    end if

  end function distribution_weight_mass

end module distributions
