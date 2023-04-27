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
	character(len = 100) :: filename, write_filename


	!  Ge
	!	  config = get_config('config')
	!	  exp_data = get_exp_data('EXP.dat')
	!		abundance_data = get_abundance_data(config)
	!		res = add_spectra(exp_data, config)

	!		Al
		config = get_config('config2')
		exp_data = get_exp_data('Al_EXP.dat')

		abundance_data = get_abundance_data(config)
		res = add_spectra(exp_data, config)

	!	 Al 单个测试
!	config = get_config('config2')
!	exp_data = get_exp_data('Al_EXP.dat')
!	cal_data = get_cal_data('spectra.dat')
!	result_data = widen(config, exp_data, cal_data, 34.6d0, 2.43d21)
!	open(22, file = 'widen_Al3.dat', status = 'unknown')
!	do j = 1, result_data%num
!		write(22, *) result_data%wavelength_nm(j), result_data%gaussian(j)
!	end do
!	close(22)
!		do i = 3, 6
!			write(filename, '(A10, I1, A4)') 'spectra_Al', i, '.dat'
!			write(write_filename, '(A8, I1, A4)') 'widen_Al', i, '.dat'
!			cal_data = get_cal_data(filename)
!			result_data = widen(config, exp_data, cal_data, 34.6d0, 2.43d21)
!			open(22, file = write_filename, status = 'unknown')
!			do j = 1, result_data%num
!				write(22, *) result_data%wavelength_nm(j), result_data%cross_P(j)
!			end do
!			close(22)
!		end do

	!	result_data = widen(config, exp_data, cal_data, 21d0, 1d21)
	!	open(22, file = 'widen_Al_3.dat', status = 'unknown')
	!	do i = 1, result_data%num
	!		write(22, *) result_data%wavelength_nm(i), result_data%cross_P(i)
	!	end do

end program Simulations

