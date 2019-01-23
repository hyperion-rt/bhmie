program main

  use lib_array
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

  real(dp),allocatable :: density(:), abundance_mass(:), wavelengths(:)
  ! density of the various components

  integer :: ic
  ! loop variable for components

  integer :: output_format

  real(dp) :: gas_to_dust

  integer :: n_wav
  real(dp) :: wav_min, wav_max


  integer :: n_threads=1 !needed for omp parallel version

  ! retrieve parameter file from command-line
  call get_command_argument(1, input_file)
  if(trim(input_file)=='') stop "Usage: bhmie input_file"

  ! open parameter file, and read
  open(unit=32, file=input_file, status='old')

  ! read header line
  read(32,*) prefix
  read(32,*) output_format
  read(32,*) amin
  read(32,*) amax
  read(32,*) na
  read(32,*) n_angles
  read(32,*) n_small_angles
  read(32,*) n_components
  read(32,*) gas_to_dust
  read(32,*) wav_min, wav_max, n_wav
  !$ read(32,*) n_threads

  !$ print*,"Running in parallel with ",n_threads," threads"

  allocate(wavelengths(n_wav))
  call logspace(wav_min, wav_max, wavelengths)

  ! allocate arrays
  allocate(abundance_mass(n_components))
  allocate(density(n_components))
  allocate(m(n_components))
  allocate(d(n_components))

  do ic=1,n_components
     read(32,*)
     read(32,*) abundance_mass(ic)
     read(32,*) density(ic)
     call read_material(32, m(ic))
     call interpolate_material(m(ic), wavelengths)
     call read_distribution(32, d(ic))
  end do

  close(unit=32)

  ! Re-normalize mass abundance
  abundance_mass = abundance_mass / sum(abundance_mass)

  call compute_dust_properties(prefix,output_format,abundance_mass,m,d,density, &
       & gas_to_dust,amin,amax,na,wavelengths,n_angles,n_small_angles,n_threads)

end program main
