program test

  use distributions
  use materials
  use types
  use bhmie_wrapper

  implicit none

  character(len=1000) :: input_file
  ! the parameter file

  type(material),allocatable :: m(:)
  ! materials used to compute the physical properties

  type(size_distribution),allocatable :: d(:)
  ! dust size distributions used to compute the physical properties

  integer :: n_angles, n_small_angles
  ! number of angles to compute scattering matrix for, and number of additional fine angles

  integer :: n_components
  ! the number of components

  character(len=100) :: prefix
  ! output filename

  integer :: na
  real(dp) :: amin, amax
  ! overall size range and number of size bins to use 

  real(dp),allocatable :: density(:), abundance(:), wavelengths(:)
  ! density of the various components

  integer :: ic
  ! loop variable for components

  integer :: output_format

  ! retrieve parameter file from command-line
  call get_command_argument(1, input_file)
  if(trim(input_file)=='') stop "Usage: bhmie input_file"

  ! open parameter file, and read
  open(unit=32, file=input_file)

  ! read header line
  read(32,*) prefix
  read(32,*) output_format
  read(32,*) amin
  read(32,*) amax
  read(32,*) na
  read(32,*) n_angles
  read(32,*) n_small_angles
  read(32,*) n_components

  ! allocate arrays
  allocate(abundance(n_components))
  allocate(density(n_components))
  allocate(m(n_components))
  allocate(d(n_components))

  do ic=1,n_components
     read(32,*)
     read(32,*) abundance(ic)
     read(32,*) density(ic)
     call read_material(32, m(ic))
     call read_distribution(32, d(ic))
  end do

  close(unit=32)

  do ic=2,n_components
     if(size(m(1)%wavelengths).ne.size(m(ic)%wavelengths)) then
        print *,'Size of wavelength arrays do not match for different materials'
        stop
     end if
     if(any(m(1)%wavelengths.ne.m(ic)%wavelengths)) then
        print *,'Wavelength arrays do not match for different materials'
        stop
     end if
  end do

  allocate(wavelengths(size(m(1)%wavelengths)))
  wavelengths = m(1)%wavelengths

  call compute_dust_properties(prefix,output_format,abundance,m,d,density,amin,amax,na,wavelengths,n_angles,n_small_angles)

end program test
