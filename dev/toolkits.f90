module toolkits
use ifport

contains
      
! Export 2d double array in text file
  subroutine Export2DArray(DArray, DataFileName)
  real*8, dimension(:,:), intent(in) :: DArray
  character*(*), intent(in) :: DataFileName
  integer, dimension(2) :: DArrayShape
  integer :: i, j
! The DArray must have no more than 999 columns, else an error mark writes in data file
  character*3 :: ColNum_str
  character*10 :: fmt_str
  open(11, file = DataFileName, status='replace')  
  DArrayShape = shape(DArray)
  if (DArrayShape(2) .LE. 999) then
    write(ColNum_str, '(i3)') DArrayShape(2)
    fmt_str = '('//trim(ColNum_str)//'E12.5)'
    do i = 1, DArrayShape(1)
      write(11, fmt_str) (DArray(i,j), j = 1, DArrayShape(2))
    end do
  else
    write(1, *) 'The given array has more than 999 columns.'
  end if
  write(11, *)
  write(11, *) 'Successfully generate data file at '//clock()//' '//date()
  close(11)
  end subroutine
  
! Profile exporter
  subroutine ProfileExport(XLoc, YLoc, ZVal, DataFileName)
  real*8, dimension(:), intent(in) :: XLoc, YLoc
  real*8, dimension(:,:), intent(in) :: ZVal
  character*(*), intent(in) :: DataFileName
  integer, dimension(2) :: ZVal_shape
  integer :: NumPtX, NumPtY, i, j
  open(11, file = DataFileName, status = 'replace')
  NumPtX = size(XLoc)
  NumPtY = size(YLoc)
  ZVal_shape = shape(ZVal)
! Check the size of location arrays
  if ((NumPtX .eq. ZVal_shape(1)) .and. (NumPtY .eq. ZVal_shape(2))) then
    do i = 1, NumPtX
      do j = 1, NumPtY
        write(11, '(3e12.5)') XLoc(i), YLoc(j), ZVal(i,j)
      end do
    end do
    write(11, *)
    write(11, *) 'Successfully generate data file at '//clock()//' '//date()
  else
    write(11, *) 'The specified location arrays are not consistent with the given profile data!'
  end if
  close(11) 
  end subroutine

  Logical function CheckRange(eval, lb, ub) ! if the value in the range of (lb, ub)
    real*8, intent(in) :: eval, lb, ub
    if ((eval .ge. lb) .and. (eval .le. ub)) then
        CheckRange = .true.
    else
        CheckRange = .false.
    end if
  end function
  
  ! Check the unit of input temperature as declaration of unit 
  subroutine CheckTempUnit(unit, Tin, Tout)
    character*1, intent(in) :: unit
    real*8, intent(in) :: Tin
    real*8, intent(out) :: Tout
    real*8 :: T0 = 273.15
    select case(unit)
      case('K')
        if (Tin .LE. T0) then
          Tout = Tin+T0
        else
          Tout = Tin
        end if
      case('C')
        if (Tin .LE. T0) then
          Tout = Tin
        else
          Tout = Tin-T0
        end if
    end select
  end subroutine

end module