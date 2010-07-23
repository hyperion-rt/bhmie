module distributions

  use types
  implicit none
  save

  type size_distribution
     integer :: type
     ! type = 1: power-law
     ! type = 2: PED (not implemented)
     ! type = 3: numerical (not implemented)
     real(dp) :: amin, amax, apower
  end type size_distribution

contains

  subroutine read_distribution(unit, d)
    implicit none
    integer,intent(in) :: unit
    type(size_distribution),intent(out) :: d
    character(len=100) :: distribution_type
    read(unit,*) distribution_type
    backspace(unit)
    select case(trim(distribution_type))
    case('power')
       read(unit,*) distribution_type, d%amin, d%amax, d%apower
       d%type = 1
    case default
       stop "Not implemented"
    end select
  end subroutine read_distribution

  real(dp) function distribution_weight_number(d, a1, a2)
    implicit none
    type(size_distribution),intent(in) :: d
    real(dp),intent(in) :: a1, a2
    select case(d%type)
    case(1)
       if(a2 < d%amin .or. a1 > d%amax) then
          distribution_weight_number = 0._dp
       else if(a1 < d%amin .and. a2 > d%amax) then
          distribution_weight_number = 1._dp
       else if(a1 < d%amin) then
          distribution_weight_number = (a2**(d%apower+1._dp) - d%amin**(d%apower+1._dp)) &
               &                     / (d%amax**(d%apower+1._dp) - d%amin**(d%apower+1._dp))    
       else if(a2 > d%amax) then
          distribution_weight_number = (d%amax**(d%apower+1._dp) - a1**(d%apower+1._dp)) &
               &                     / (d%amax**(d%apower+1._dp) - d%amin**(d%apower+1._dp))
       else
          distribution_weight_number = (a2**(d%apower+1._dp) - a1**(d%apower+1._dp)) &
               &                     / (d%amax**(d%apower+1._dp) - d%amin**(d%apower+1._dp))
       end if
    case default
       stop "Not implemented"
    end select
  end function distribution_weight_number

  real(dp) function distribution_weight_mass(d, a1, a2)
    implicit none
    type(size_distribution),intent(in) :: d
    real(dp),intent(in) :: a1, a2
    select case(d%type)
    case(1)
       if(a2 < d%amin .or. a1 > d%amax) then
          distribution_weight_mass = 0._dp
       else if(a1 < d%amin .and. a2 > d%amax) then
          distribution_weight_mass = 4._dp
       else if(a1 < d%amin) then
          distribution_weight_mass = (a2**(d%apower+4._dp) - d%amin**(d%apower+4._dp)) &
               &                   / (d%amax**(d%apower+4._dp) - d%amin**(d%apower+4._dp))    
       else if(a2 > d%amax) then
          distribution_weight_mass = (d%amax**(d%apower+4._dp) - a1**(d%apower+4._dp)) &
               &                   / (d%amax**(d%apower+4._dp) - d%amin**(d%apower+4._dp))
       else
          distribution_weight_mass = (a2**(d%apower+4._dp) - a1**(d%apower+4._dp)) &
               &                   / (d%amax**(d%apower+4._dp) - d%amin**(d%apower+4._dp))
       end if
    case default
       stop "Not implemented"
    end select
  end function distribution_weight_mass

end module distributions
