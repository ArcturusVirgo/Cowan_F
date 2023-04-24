program main
  implicit none
  integer, parameter :: N = 1024
  real(kind = 8), dimension(N) :: x1, y1, x2, y2
  integer :: i, j
  real(kind = 8) :: SS_reg, SS_tot, R2
  
  open(21, file = 'EXP.dat', status = 'unknown')
  open(22, file = 'result-0.3998-1.dat', status = 'unknown')
  do i = 1, N
    read(21, *) x1(i), y1(i)
    read(22, *) x2(i), y2(i)
  end do
  
  SS_reg = 0
  SS_tot = 0
  y1 = y1 / maxval(y1)
  y2 = y2 / maxval(y2)
!  do j = 1, N
!    write(*, *) x2(j), y2(j), y1(j)
!  end do
  do j = 1, N
    SS_reg = SS_reg + (y2(j) - sum(y1) / N)**2
    SS_tot = SS_tot + (y1(j) - sum(y1) / N)**2
  end do
  R2 = SS_reg / SS_tot
  if (R2 > 1) then
    R2 = 1 / R2
  end if
  
  write(*, *) R2

end program main