program Simulations
  use DFile_Mod
  use cowan_data_processing
  implicit none
  
  type(t_exp_data) :: exp_data
  type(t_cal_data) :: cal_data
  type(t_widen_data) :: result_data
  type(t_config) :: config
  type(t_abundance_data) :: abundance_data
  type(t_add_spectra) :: res
  integer :: i, j  ! 循环变量
  
  config = get_config('config')
  !  write(*, *) config%abundance_filename
  exp_data = get_exp_data('EXP.dat')
  abundance_data = get_abundance_data(config)
  !  do i = 1, abundance_data%num
  !    if (sum(abundance_data%abundance(i, :)) < 0.98) then
  !
  !      write(*, *) i, sum(abundance_data%abundance(i, :))
  !    end if
  !  end do
  !  cal_data = get_cal_data('spectra-Ge4.dat')
  !  result_data = widen(config, exp_data, cal_data, 15d0, 1d22)
  !  do i=1,result_data%num
  !    write(*,*) result_data%wavelength_nm(i), result_data%gaussian(i)
  !  end do
  
    res = add_spectra(exp_data, config)

end program Simulations

