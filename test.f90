program main
  use ion_abundance
  implicit none
  type(t_init_data) :: init_data
  type(t_abu) :: abu
  init_data = get_init_data('ipinput.dat')
  abu = cal_abu(init_data)
end program main

