module cowan_data_processing
  use DFile_Mod
  implicit none
  
  ! =========== 结构体定义 ===========
  ! 实验数据
  type :: t_exp_data
    integer :: num
    real(kind = 8), dimension(:), allocatable :: wavelength_nm, wavelength_ev, intensity
  end type t_exp_data
  ! 计算数据
  type :: t_cal_data
    integer :: num
    integer, dimension(:), allocatable :: index_l, index_h
    real(kind = 8), dimension(:), allocatable :: energy_l, energy_h, wavelength_ev, intensity, J_l, J_h
  end type t_cal_data
  ! 结果数据
  type :: t_widen_data
    integer :: num
    real(kind = 8), dimension(:), allocatable :: wavelength_ev, wavelength_nm, gaussian, cross_P, cross_NP
  end type t_widen_data
  ! 粒子丰度数据
  type :: t_abundance_data
    integer :: num
    real(kind = 8), dimension(:), allocatable :: temperature, Ne
    real(kind = 8), dimension(:, :), allocatable :: abundance
  end type t_abundance_data
  ! 基础配置
  type :: t_config
    integer :: ion_num
    real(kind = 8) :: fwhmgauss, delta, min_wavelength_nm, max_wavelength_nm
    character(len = 20), dimension(:), allocatable :: ion_filename
    character(len = 20) :: abundance_filename
  end type t_config
  ! 叠加光谱
  type t_add_spectra
    integer :: num, col_num
    real(kind = 8), dimension(:, :), allocatable :: wavelength_ev, wavelength_nm, intensity
  end type t_add_spectra
  
  real(kind = 8) :: pi = 3.1415926535897


contains
  
  function get_exp_data(filename) result(exp_data)
    implicit none
    character(len = *), intent(in) :: filename
    type(t_exp_data) :: exp_data
    integer :: i
    
    open(21, file = filename, status = 'old')
    exp_data%num = GetFileN(21)
    allocate(exp_data%wavelength_nm(exp_data%num))
    allocate(exp_data%wavelength_ev(exp_data%num))
    allocate(exp_data%intensity(exp_data%num))
    do i = 1, exp_data%num
      read(21, *) exp_data%wavelength_nm(i), exp_data%intensity(i)
      exp_data%wavelength_ev(i) = 1239.85 / exp_data%wavelength_nm(i)
    end do
    close(21)
  
  end function get_exp_data
  
  function get_cal_data(filename) result(cal_data)
    implicit none
    character(len = *), intent(in) :: filename
    type(t_cal_data) :: cal_data
    integer :: i
    
    open(21, file = filename, status = 'old')
    cal_data%num = GetFileN(21)
    allocate(cal_data%index_l(cal_data%num))
    allocate(cal_data%index_h(cal_data%num))
    allocate(cal_data%energy_l(cal_data%num))
    allocate(cal_data%energy_h(cal_data%num))
    allocate(cal_data%wavelength_ev(cal_data%num))
    allocate(cal_data%intensity(cal_data%num))
    allocate(cal_data%J_l(cal_data%num))
    allocate(cal_data%J_h(cal_data%num))
    do i = 1, cal_data%num
      read(21, *) cal_data%energy_l(i), cal_data%energy_h(i), cal_data%wavelength_ev(i), &
          cal_data%intensity(i), cal_data%index_l(i), cal_data%index_h(i), cal_data%J_l(i), cal_data%J_h(i)
    end do
    close(21)
  
  end function get_cal_data
  
  function get_abundance_data(config) result(abu_data)
    implicit none
    type(t_config), intent(in) :: config
    integer :: i
    
    type(t_abundance_data) :: abu_data
    open(21, file = config%abundance_filename, status = 'old')
    abu_data%num = GetFileN(21)
    allocate(abu_data%temperature(abu_data%num))
    allocate(abu_data%Ne(abu_data%num))
    allocate(abu_data%abundance(abu_data%num, config%ion_num))
    do i = 1, abu_data%num
      read(21, *) abu_data%temperature(i), abu_data%Ne(i), abu_data%abundance(i, :)
    end do
  end function get_abundance_data
  
  function get_config(filename) result(config)
    implicit none
    character(len = *), intent(in) :: filename
    type(t_config) :: config
    integer :: i
    
    open(21, file = filename, status = 'old')
    read(21, *)
    read(21, *) config%min_wavelength_nm, config%max_wavelength_nm
    read(21, *)
    read(21, *) config%fwhmgauss
    read(21, *)
    read(21, *) config%delta
    read(21, *)
    read(21, *) config%ion_num
    allocate(config%ion_filename(config%ion_num))
    read(21, *)
    do i = 1, config%ion_num
      read(21, *) config%ion_filename(i)
    end do
    read(21, *)
    read(21, *) config%abundance_filename
    
    close(21)
  
  end function get_config
  
  function widen(config, exp_data, cal_data, temperature, Ne) result(result_data)
    ! 实验数据的结构体 计算数据结构体 温度 半高宽度 偏移(nm)
    implicit none
    type(t_exp_data), intent(in) :: exp_data
    type(t_cal_data), intent(in) :: cal_data
    type(t_config), intent(in) :: config
    real(kind = 8), intent(in) :: temperature, Ne
    
    type(t_widen_data) :: result_data
    real(kind = 8) :: max_wavelength_ev, min_wavelength_ev
    real(kind = 8), dimension(:), allocatable :: new_wavelength, new_itenstiy, new_energy, new_J, population
    real(kind = 8) :: uu, ss, tt, fwhmgauss, delta
    integer :: num
    integer :: i, j
    real(kind = 8) :: min_energy, min_J
    
    ! 获取特定值
    fwhmgauss = config%fwhmgauss
    delta = config%delta
    
    ! 确定展宽范围
    max_wavelength_ev = 1239.85 / config%min_wavelength_nm
    min_wavelength_ev = 1239.85 / config%max_wavelength_nm
    
    ! 找到下态最小能量和最小能量对应的J
    min_energy = cal_data%energy_l(1)
    min_J = cal_data%J_l(1)
    do i = 2, cal_data%num
      if (min_energy > cal_data%energy_l(i)) then
        min_energy = cal_data%energy_l(i)
        min_J = cal_data%J_l(i)
      end if
    end do

    ! 筛选波长范围在实验数据范围内的跃迁正例个数
    num = 0
    do i = 1, cal_data%num
      if (abs(cal_data%wavelength_ev(i)) < max_wavelength_ev .and. abs(cal_data%wavelength_ev(i)) > min_wavelength_ev) then
        num = num + 1
      end if
    end do

    ! 分配内存
    allocate(new_wavelength(num))
    allocate(new_itenstiy(num))
    allocate(new_energy(num))
    allocate(new_J(num))
    allocate(population(num))
    
    ! 重新赋值
    num = 0
    do i = 1, cal_data%num
      ! 筛选波长范围在实验数据范围内的跃迁正例
      if (abs(cal_data%wavelength_ev(i)) < max_wavelength_ev .and. abs(cal_data%wavelength_ev(i)) > min_wavelength_ev) then
        num = num + 1
        new_wavelength(num) = abs(1239.85 / (1239.85 / cal_data%wavelength_ev(i) - delta))
        new_itenstiy(num) = abs(cal_data%intensity(i))
        ! 取上态和下态中能量较高的那个和对应的J
        if (cal_data%energy_l(i) > cal_data%energy_h(i)) then
          new_energy(num) = cal_data%energy_l(i)
          new_J(num) = cal_data%J_l(i)
        elseif (cal_data%energy_h(i) > cal_data%energy_l(i)) then
          new_energy(num) = cal_data%energy_h(i)
          new_J(num) = cal_data%J_h(i)
        end if
      end if
    end do


    
    ! 计算布居
    do i = 1, num
      population(i) = (2 * new_J(i) + 1) * exp(-abs(new_energy(i) - min_energy) * 0.124 / temperature) / (2 * min_J + 1)
    enddo
    
    ! 初始化，准备展宽
    result_data%num = exp_data%num
    allocate(result_data%wavelength_nm(result_data%num))
    allocate(result_data%wavelength_ev(result_data%num))
    allocate(result_data%gaussian(result_data%num))
    allocate(result_data%cross_NP(result_data%num))
    allocate(result_data%cross_P(result_data%num))
    result_data%wavelength_nm = exp_data%wavelength_nm
    result_data%wavelength_ev = exp_data%wavelength_ev
    ! 开始展宽
    do i = 1, result_data%num
      tt = 0
      ss = 0
      uu = 0
      do j = 1, num
        tt = tt + new_itenstiy(j) / sqrt(2 * pi) / fwhmgauss * 2.355 * &
            exp(-2.355**2 * (new_wavelength(j) - result_data%wavelength_ev(i))**2 / fwhmgauss**2 / 2)
        ss = ss + (new_itenstiy(j) / (2 * new_J(j) + 1)) * (2 * fwhmgauss) / &
            (2 * pi * (((new_wavelength(j) - result_data%wavelength_ev(i))**2 + (2 * fwhmgauss)**2 / 4)))
        uu = uu + (new_itenstiy(j) * population(j) / (2 * new_J(j) + 1)) * (2 * fwhmgauss) / &
            (2 * pi * (((new_wavelength(j) - result_data%wavelength_ev(i))**2 + (2 * fwhmgauss)**2 / 4)))
      end do
      result_data%gaussian(i) = tt
      result_data%cross_NP(i) = ss
      result_data%cross_P(i) = uu

    end do
  end function widen
  
  function add_spectra(exp_data, config)
    implicit none
    type(t_exp_data), intent(in) :: exp_data
    type(t_config), intent(in) :: config
    
    type(t_abundance_data) :: abu_data
    type(t_widen_data) :: temp_widen_data
    type(t_cal_data) :: temp_cal_data
    
    type(t_add_spectra) :: add_spectra
    
    integer :: i, j
    real(kind = 8) :: temperature, Ne, sigma
    character(len = 50) :: filename, save_filename
    real(kind = 8), dimension(:), allocatable :: y1, y2
    abu_data = get_abundance_data(config)
    
    add_spectra%num = exp_data%num
    add_spectra%col_num = abu_data%num
    allocate(add_spectra%wavelength_nm(add_spectra%num, add_spectra%col_num))
    allocate(add_spectra%wavelength_ev(add_spectra%num, add_spectra%col_num))
    allocate(add_spectra%intensity(add_spectra%num, add_spectra%col_num))
    
    do i = 1, abu_data%num
      add_spectra%wavelength_nm(:, i) = exp_data%wavelength_nm
      add_spectra%wavelength_ev(:, i) = exp_data%wavelength_ev
      add_spectra%intensity(:, i) = 0
      temperature = abu_data%temperature(i)
      Ne = abu_data%Ne(i)
      do j = 1, config%ion_num
        filename = config%ion_filename(j)
        temp_cal_data = get_cal_data(filename)
        temp_widen_data = widen(config, exp_data, temp_cal_data, temperature, Ne)
        add_spectra%intensity(:, i) = add_spectra%intensity(:, i) + temp_widen_data%cross_P * abu_data%abundance(i, j)
      end do
      ! 拟合的R2计算
      allocate(y1(add_spectra%num), y2(add_spectra%num))
      y1 = exp_data%intensity/ maxval(exp_data%intensity)
      y2 = add_spectra%intensity(:, i) / maxval(add_spectra%intensity(:, i))
      sigma = sqrt(sum((y2-y1)**2)/add_spectra%num)
      deallocate(y1)
      deallocate(y2)
      
      if (i<10)then
        write(save_filename, '(A, F6.4, A, I1, A)') 'result-', sigma, '-', i, '.dat'
      elseif(i>=10 .and. i<100)then
        write(save_filename, '(A, F6.4, A, I2, A)') 'result-', sigma, '-', i, '.dat'
      elseif(i>=100 .and. i<1000)then
        write(save_filename, '(A, F6.4, A, I3, A)') 'result-', sigma, '-', i, '.dat'
      elseif(i>=1000 .and. i<10000)then
        write(save_filename, '(A, F6.4, A, I3, A)') 'result-', sigma, '-', i, '.dat'
      end if
      open(22, file = save_filename, status = 'unknown')
      do j = 1, add_spectra%num
        write(22, *) add_spectra%wavelength_nm(j, i), add_spectra%intensity(j, i)
      end do
      close(22)
      write(*, *) i, save_filename, '已保存'
    
    end do
  
  end function add_spectra

end module cowan_data_processing