module materials

  use types
  implicit none
  save

  type material
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
    open(unit=33,file=filename)
    do
       read(33,*,iostat=ioerr)
       if(ioerr.ne.0) exit
       n_wav = n_wav + 1
    end do
    close(unit=33)

    ! Allocate arrays
    allocate(m%wavelengths(n_wav))
    allocate(m%refractive_indices(n_wav))

    ! Read in values
    open(unit=33,file=filename)
    do iw=1,n_wav
       read(33,*) m%wavelengths(iw), ref_real, ref_imag
       m%refractive_indices(iw) = cmplx(ref_real, ref_imag)
    end do
    close(unit=33)

  end subroutine read_material

end module materials
