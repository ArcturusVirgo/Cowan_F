program ion_fraction
    implicit none
    
    ! input the Atomic index and initial temperature
    integer, parameter :: N = 32
    integer :: i,Te
    real :: pNe,Ratsum
    
    integer Z(N), pz(N), k, mm
    real pez(N)
    real S(N), Ar(N), A3r(N), Rat(N + 1), CoM(N), CoMM(N), Cosum, &
            Co(N), Ratio(N + 1), Rsum(N), RRR(10000, N + 1), &
            TTT(N + 1, 10000), aveion, kk
    
    ! 定义输入文件（电荷态，电离能，最外层电子数），输出文件（比率）
    open(1, file = 'ipinput.dat', status = 'old')
    open(2, file = 'rationoutput.dat', status = 'unknown')
    
    ! 读入初始数据
    read(1, *) (Z(i), pez(i), pz(i), i = 1, N)
    ! write(2, 40) (Z(i), pez(i), pz(i), i = 1, N)
    40 format(1x, I2, 2x, f11.4, 2x, I2)
    
    ! write(2, *) 'output the result:'
    
    do Te = 15, 45, 1
        k = 0
        do kk = 17, 22, 1
            pNe = 10.0**kk
            ! 定义函数S(Z, Te), Ar and A3r
            do i = 1, N
                S(i) = (9 * (1e-6) * pz(i) * sqrt(Te / pez(i)) * exp(-pez(i) / Te)) / &
                        (pez(i)**(1.5) * (4.88 + Te / pez(i)))
                Ar(i) = 5.2 * (1e-14) * sqrt(Pez(i) / Te) * Z(i) * (0.429 + 0.5 * &
                        LOG(Pez(i) / Te) + 0.469 * sqrt(Te / Pez(i)))
                A3r(i) = 2.97 * (1e-27) * pz(i) / (Te * pez(i)**2 * (4.88 + Te / pez(i)))
            end do
            
            ! 求解过程-----------
            ! 比例系数
            do i = 1, N
                Co(i) = S(i) / (Ar(i) + PNe * A3r(i))
            end do
            CoM(1) = Co(1)
            do i = 2, N
                CoM(i) = Co(i) * CoM(i - 1)
            end do
            Cosum = 0
            do i = 1, N
                CoMM(i) = i * CoM(i)
                Cosum = Cosum + CoMM(i)
            end do
            Cosum = Cosum
            ! 求解离子数密度
            Rat(1) = PNe / Cosum
            
            do i = 2, N + 1
                Rat(i) = Co(i - 1) * Rat(i - 1)
            enddo
            ! write(2, 50)Te,(Rat(l), l = 2, N+1)
            ! 对离子数密度求和，得到总的离子密度
            Ratsum = 0
            do i = 1, N + 1
                Ratsum = Ratsum + Rat(i)
            enddo
            aveion = PNe / Ratsum
            !	write(2,*) PNe,aveion
            ! 得到离子丰度
            
            k = k + 1
            do i = 1, N + 1
                Ratio(i) = Rat(i) / Ratsum
                if (Ratio(i) < 1e-6) then
                    Ratio(i) = 0
                endif
            enddo
            write(2, 50)Te, PNe, (Ratio(i), i = 5, 14)
        
        end do
    end do
    
    48 format(I4, 2x, 47f17.5)
    49 format(47f17.5)
    50 format(f6.3, e14.5, 2x, 48f17.5)

end program ion_fraction