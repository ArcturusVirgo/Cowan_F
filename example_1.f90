
    program Resonance


	parameter(m=180000,n=1000,delta=0.03)
!	parameter(m=180000,n=1000,delta=0.04)
    real wavelength(m),gf(m),specen(m),fwhm(m),wfwhm(m),lj,&
         aversum,averintensity,average,wspecen(m),fwspecen(m),lowk(m),highk(m),&
         wave(n),intensity(n),fwhmnew(m),wavenew(m),newwspecen(m),&
         intensitynew(m),tt,emin,emax,cross(n),ss,flowk(m),fwhmgauss,&
         fhighk(m),fjt(m),fjtp(m),ffjtp(m),ffjt(m),lowenergy(m),&
         flowenergy(m),minimum,crossb(n),population(m),uu,ffwhm(m),newfjt(m),newfjtp(m)


	emin=70
	emax=300
    Te=25
    fwhmgauss=0.2
!   open(10,file="Jenergy-fwhm.dat",status="old")
	open(11,file="spectra.dat",status="old")

!	open(11,file="test 1- 5.txt",status="old")
!	open(18,file="test 1- 5.dat",sta tus="unknown")


!   open(15,file="tt",status="unknown")
!   open(16,file="out.dat",status="unknown")

    open(17,file="gauss-25eV.dat",status="unknown")
    open(18,file="cross-Np-25eV.dat",status="unknown")
    open(19,file="cross-P-25eV.dat",status="unknown")

!   Read initial data
!	i=0
!10	i=i+1
!     read(10,"(f9.3,2x,f9.3)",end=20) specen(i),fwhm(i)
!     goto 10
!20	  continue
!       energynum=i-1
!       write(15,"(f9.3,2x,f9.3)") (specen(i),fwhm(i),i=1, energynum)
	i=0
30	i=i+1
     read(11,*,end=40) lowenergy(i),&
     ! read(11,"(3f9.3,f9.4,2I3,2f5.1)",end=40) lowenergy(i),
     wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
     goto 30
40	continue
      trans=i-1


!      find the minimum energy in the lower levels for boltzman law

      minimum=lowenergy(1)
      do j=2,trans
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!      fwhmgauss=0.27

!       write(15,*) minimum

!      write(16,"(3f9.3)") (wspecen(i),wavelength(i),gf(i),i=1,trans)
      do k=1,trans
!       do j=1,energynum
!          if(wspecen(k).eq.specen(j)) then
!	   if(abs(specen(k)-wspecen(j)).lt.0.001) then
!	         wfwhm(k)=fwhm(j)
!                ffwhm(k)=fwhm(j)
!                  if (wfwhm(k).lt.fwhmgauss) then
                  wfwhm(k)=fwhmgauss
!                   endif
!	    endif
!	 enddo
!	write(16,"(f9.3,f9.4,f9.3,2x,2I3,2f5.1)")wavelength(k),gf(k),
!    X   ffwhm(k)*2,lowk(k),highk(k),fjt(k),fjtp(k)
      enddo

!     Line by line broaden code
!     # define both two-dimension variables: wave and intensity
!     # define one dimension of sum of intensity: intensiysum
      mm=0
!	emin=70
!	emax=200
      do 100 i=1,trans
      if(abs(wavelength(i)).gt. emin .and. abs(wavelength(i)).lt.emax)then
		      mm=mm+1
			wavenew(mm)=abs(1239.85/(1239.85/wavelength(i)-delta))
			intensitynew(mm)=abs(gf(i))
            !fwhmnew(mm)=wfwhm(i)*2
            fwhmnew(mm)=wfwhm(i)
	        flowk(mm)=lowk(i)
	        fhighk(mm)=highk(i)
	        newfjtp(mm)=fjtp(i)
            newfjt(mm)=fjt(i)
            fwspecen(mm)=wspecen(i)
            flowenergy(mm)=lowenergy(i)
		end if
100   continue

       !     consider the Boltzman contribution

!      do Te=15,15
         do i=1,mm
       if(flowenergy(i)>fwspecen(i))then
           newwspecen(i)=flowenergy(i)
           ffjtp(mm)=newfjt(i)
	else
           newwspecen(i)=fwspecen(i)
           ffjtp(mm)=newfjtp(i)
	end if

      population(i)=(2*ffjtp(i)+1)*exp(-abs(newwspecen(i)-minimum)*0.124/Te)/(2*lj+1)

	    write(15,"(f10.4)") population(i)
         enddo
!	  enddo


      do i=1,n
       wave(i)=emin+(emax-emin)*i/n
            intensity(i)=0
            cross(i)=0
            crossb(i)=0
	enddo



      do i=1,n
       tt=0
	   ss=0
	   uu=0
      do j=1,mm
!       tt=tt+(intensitynew(j)/(2*ffjtp(j)+1))/sqrt(2*3.14158)
!     	 /fwhmgauss*2.355*exp(-2.355**2*(wavenew(j)-wave(i))**2
!         /fwhmgauss**2/2)

       tt=tt+intensitynew(j)/sqrt(2*3.14158)&
     	 &/fwhmgauss*2.355*exp(-2.355**2*(wavenew(j)-wave(i))**2&
          &/fwhmgauss**2/2)




	 ss=ss+(intensitynew(j)/(2*ffjtp(j)+1))*fwhmnew(j)/&
           &(2*3.1415928*(((wavenew(j)-wave(i))**2+fwhmnew(j)**2/4)))

       uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	 &fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	 &fwhmnew(j)**2/4)))

      enddo
            intensity(i)=tt
            cross(i)=ss
            crossb(i)=uu

	write(6,*)  i,crossb(i),cross(i),intensity(i)
	enddo



      do i=1,n
	     if(intensity(i).lt.0.00001 ) then
           intensity(i)=0
	     endif
	     write(17,*) 1239.85/wave(i), intensity(i)
	enddo

      do i=1,n
           if(cross(i).lt.0.00001 ) then
           cross(i)=0
	     endif
	     write(18,*) 1239.85/wave(i), cross(i)
	enddo


      do i=1,n
           if(crossb(i).lt.0.00001 ) then
           crossb(i)=0
	     endif
	     write(19,*) 1239.85/wave(i), crossb(i)
	enddo

      end
