module materials

  use lib_array
  use types
  implicit none
  save

  type material
     character(len=100) :: filename
     real(dp),allocatable :: wavelengths(:)
     complex(dp),allocatable :: refractive_indices(:)
  end type material

contains

  subroutine read_material(unit, m)

    implicit none

    integer,intent(in) :: unit
    type(material),intent(out) :: m

    character(len=1000) :: filename
    integer :: n_wav,ioerr,iw
    real(dp) :: ref_real, ref_imag

    read(unit, *) filename

    ! Find number of wavelengths
    n_wav = 0
    open(unit=33,file=filename, status='old')
    do
       read(33,*,iostat=ioerr)
       if(ioerr.ne.0) exit
       n_wav = n_wav + 1
    end do
    close(unit=33)

    m%filename = filename

    ! Allocate arrays
    allocate(m%wavelengths(n_wav))
    allocate(m%refractive_indices(n_wav))

    ! Read in values
    open(unit=33,file=filename, status='old')
    do iw=1,n_wav
       read(33,*) m%wavelengths(iw), ref_real, ref_imag
       m%refractive_indices(iw) = cmplx(ref_real, ref_imag)
    end do
    close(unit=33)

  end subroutine read_material

  subroutine interpolate_material(m, wav)

    implicit none

    type(material),intent(inout) :: m
    real(dp),intent(in) :: wav(:)
    real(dp) :: ref_real(size(wav)), ref_imag(size(wav))
    integer :: iw

    if(maxval(wav) > maxval(m%wavelengths)) then
       write(*,'("ERROR: Refractive index file ",A," only goes up to ",F10.4," microns")') trim(m%filename), maxval(m%wavelengths)
       stop
    end if

    if(minval(wav) < minval(m%wavelengths)) then
       write(*,'("ERROR: Refractive index file ",A," only goes down to ",F10.4," microns")') trim(m%filename), maxval(m%wavelengths)
       stop
    end if

    ! Interpolate
    ref_real = interp1d_loglog(m%wavelengths, real(m%refractive_indices), wav)
    ref_imag = interp1d_loglog(m%wavelengths, aimag(m%refractive_indices), wav)

    ! Deallocate current arrays
    deallocate(m%wavelengths)
    deallocate(m%refractive_indices)

    ! Reallocate with new size
    allocate(m%wavelengths(size(wav)))
    allocate(m%refractive_indices(size(wav)))

    ! Repopulate arrays
    do iw=1,size(wav)
       m%wavelengths(iw) = wav(iw)
       m%refractive_indices(iw) = cmplx(ref_real(iw), ref_imag(iw))
    end do

  end subroutine interpolate_material

end module materials
