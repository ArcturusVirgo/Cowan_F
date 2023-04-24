      program Ion Ratio


!     input the Atomic index and initial temperature
      parameter (N=32)
!

!      #########################################################
!	 使用说明：                                             #
!	    1）需创建初始输入文件ipinput.dat                    #
!        （电荷态Z，电离能Pez和最外层电子数Pz）；        #
c    2）修改元素的电子数N和电子密度PNe；                 #
c    3）输出结果文件为Rationoutput.dat;                  #
c                     2005年10月23日                     #
c#########################################################
      
      integer Z(N),pz(N),k,mm
	real pez(N)
	real S(N),Ar(N),A3r(N),Rat(N+1),CoM(N),CoMM(N),Cosum,
     X      Co(N),Ratio(N+1),Rsum(N),RRR(10000,N+1),
     X     TTT(N+1,10000),aveion,kk

c     定义输入文件（电荷态，电离能，最外层电子数），输出文件（比率）
      open(1,file='ipinput.dat',status='old')
	open(2,file='rationoutput.dat',status='unknown')
      
c     读入初始数据      
	read(1,*) (Z(i),pez(i),pz(i),i=1,N)
!      write(2,40) (Z(i),pez(i),pz(i),i=1,N)
   40 format(1x,I2,2x,f11.4,2x,I2)
      
c      write(2,*) 'output the result:'

      do 2000 mm=15,35,1
      k=0
      do 1000 kk=17,22,1
	    pNe=10.0**kk
c     定义函数S(Z,Te),Ar and A3r
	do 100 i=1,N
	S(i)=(9*(1e-6)*pz(i)*sqrt(Te/pez(i))*exp(-pez(i)/Te))/
     X   (pez(i)**(1.5)*(4.88+Te/pez(i)))
      Ar(i)=5.2*(1e-14)*sqrt(Pez(i)/Te)*Z(i)*(0.429+0.5*
     X   LOG(Pez(i)/Te)+0.469*sqrt(Te/Pez(i)))
      A3r(i)=2.97*(1e-27)*pz(i)/(Te*pez(i)**2*(4.88+Te/pez(i)))
100   continue

c求解过程-----------
c比例系数      
	do i=1,N
      Co(i)=S(i)/(Ar(i)+PNe*A3r(i))
      enddo
      CoM(1)=Co(1)
	do i=2,N
      CoM(i)=Co(i)*CoM(i-1)
      enddo
      Cosum=0
      do i=1,N
      CoMM(i)=i*CoM(i)
	Cosum=Cosum+CoMM(i)
      enddo
       Cosum=Cosum
c求解离子数密度
      Rat(1)=PNe/Cosum




	do i=2,N+1
      Rat(i)=Co(i-1)*Rat(i-1)
      enddo
c      write(2,50)Te,(Rat(l),l=2,N+1)
c对离子数密度求和，得到总的离子密度
      Ratsum=0
      do i=1,N+1
      Ratsum=Ratsum+Rat(i)
      enddo
      aveion=PNe/Ratsum
!	write(2,*) PNe,aveion
c得到离子丰度

      k=k+1
      do i=1,N+1
      Ratio(i)=Rat(i)/Ratsum
	if (Ratio(i) .LT.1e-6 ) then
	Ratio(i)=0
	endif
      RRR(k,i)=Ratio(i)
	enddo
      write(2,50)PNe,(Ratio(l),l=1,N+1)
     
1000  continue
2000  continue



      do i=1,N+1
	do j=1,k
	  TTT(i,j)=RRR(j,i)
      enddo
	enddo
	do i=2,N+1
!	write(2,*) (TTT(i,j),j=1,k)
      enddo


48    format(I4,2x,47f17.5)	

49    format(47f17.5)

50    format(e14.5,2x,48f17.5)
	end









