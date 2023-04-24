module ion_abundance
  use DFile_Mod
  implicit none
  type :: t_init_data
    integer :: N
    character(len = 2) :: atomic_name
    real, dimension(:), allocatable :: ion_num, ion_energy, electron_num
  end type t_init_data
  
  type :: t_abu
    integer :: N
    character(len = 2) :: atomic_name
    real, dimension(:), allocatable :: percentage
  end type t_abu

contains
  function get_init_data(filename) result(init_data)
    implicit none
    character(len = *), intent(in) :: filename
    type(t_init_data) :: init_data
    integer :: i
    
    open(21, file = filename, status = 'old')
    init_data%N = GetFileN(21)
    allocate(init_data%ion_num(init_data%N), init_data%ion_energy(init_data%N), init_data%electron_num(init_data%N))
    do i = 1, init_data%N
      read(21, *) init_data%ion_num(i), init_data%ion_energy(i), init_data%electron_num(i)
    end do
  end function get_init_data
  
  function cal_abu(init_data) result(abu)
    implicit none
    type(t_init_data), intent(in) :: init_data
    type(t_abu) :: abu
    
    real, dimension(:), allocatable :: ratio, S, Ar, A3r, temp
    integer :: i, j, index, Te
    real(kind=8) :: Ne
    
    abu%N = init_data%N
    abu%atomic_name = init_data%atomic_name
    !    allocate(abu%percentage(abu%N))
    allocate(ratio(abu%N - 1), S(abu%N), Ar(abu%N), A3r(abu%N), temp(abu%N))
    
    do Te = 15, 45, 1
      do index = 17, 22, 1
        Ne = 10.0**index
        do i = 1, init_data%N
          S(i) = (9 * (1e-6) * init_data%electron_num(i) * sqrt(Te / init_data%ion_energy(i)) * &
              exp(-init_data%ion_energy(i) / Te)) / &
              (init_data%ion_energy(i)**(1.5) * (4.88 + Te / init_data%ion_energy(i)))
          Ar(i) = 5.2 * (1e-14) * sqrt(init_data%ion_energy(i) / Te) * init_data%ion_num(i) * (0.429 + 0.5 * &
              LOG(init_data%ion_energy(i) / Te) + 0.469 * sqrt(Te / init_data%ion_energy(i)))
          A3r(i) = 2.97 * (1e-27) * init_data%electron_num(i) / (Te * init_data%ion_energy(i)**2 * &
              (4.88 + Te / init_data%ion_energy(i)))
        end do
        do i = 1, init_data%N - 1
          ratio(i) = S(i) / (Ar(i + 1) + Ne * A3r(i + 1))
        end do
        abu%percentage = get_percentage(ratio, abu%N)
        open(21, file = 'abu.dat', status = 'unknown')
        write(21, '(F6.3,E14.5,2x,48f17.5)') real(Te), Ne, abu%percentage(5:14)
      end do
    end do
  
  end function cal_abu
  
  function get_percentage(ratio, N) result(percentage)
    implicit none
    integer, intent(in) :: N
    real, dimension(N), intent(in) :: ratio
    real, dimension(N) :: percentage, temp
    integer :: i, j
    
    temp(1) = 1
    do i = 2, N
      do j = 1, i - 1
        temp(i) = temp(j) * ratio(j)
      end do
    end do
    do i = 1, N
      percentage(i) = 1 / sum(temp) * temp(i)
    end do
  end function get_percentage
end module ion_abundance