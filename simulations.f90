
    program Simulations


	parameter(m=100000,n=2000,delta=-0)
    real wavelength(m),gf(m),specen(m),fwhm(m),wfwhm(m),lj,&
         aversum,averintensity,average,wspecen(m),fwspecen(m),lowk(m),highk(m),&
         wave(n),intensity(n),fwhmnew(m),wavenew(m),newwspecen(m),&
         intensitynew(m),tt,emin,emax,cross(n),ss,flowk(m),fwhmgauss,&
         fhighk(m),fjt(m),fjtp(m),ffjtp(m),ffjt(m),lowenergy(m),&
         flowenergy(m),minimum,crossb(n),population(m),uu,ffwhm(m),newfjt(m),newfjtp(m),&
		 Exp_inten(3000),Exp_wave(3000),Theo_wave(10,3000),Theo_inten(10,3000),&
		 Ion4_fraction(i),Ion5_fraction(i),Ion6_fraction(i),&
         Ion7_fraction(i),Ion8_fraction(i),Ion9_fraction(i),Ion10_fraction(i),Ion11_fraction(i),&
	     Ion12_fraction(i), Ion13_fraction(i)
   integer Exp_num,fraction_num,Ion4_num


	emin=49
!   49 eV=25.3 nm
	emax=250
!   250 eV=4.96 nm
!    Te=25
    fwhmgauss=0.27

open(9,file="Ion-fraction.dat",status="unknown") 
open(10,file="EXP.dat",status="unknown")  

open(11,file="spectra-4.dat",status="unknown")
open(12,file="spectra-5.dat",status="unknown")
open(13,file="spectra-6.dat",status="unknown")
open(14,file="spectra-7.dat",status="unknown")
open(15,file="spectra-8.dat",status="unknown")
open(16,file="spectra-9.dat",status="unknown")
open(17,file="spectra-10.dat",status="unknown")
open(18,file="spectra-11.dat",status="unknown")
open(19,file="spectra-12.dat",status="unknown")
open(20,file="spectra-13.dat",status="unknown")

open(11,file="spectra-temp.dat",status="unknown")
open(19,file="cross-P.dat",status="unknown") 

! read EXP data
	i=0
91	i=i+1
     read(9,*,end=90) Exp_wave(i), Exp_inten(i)
     goto 91
90	continue
      Exp_num=i-1

! read the ion fraction
	i=0
101	i=i+1
     read(10,*,end=100) Te_value(i),Ne_value(i),Ion4_fraction(i),Ion5_fraction(i),Ion6_fraction(i),&
     Ion7_fraction(i),Ion8_fraction(i),Ion9_fraction(i),Ion10_fraction(i),Ion11_fraction(i),&
	 Ion12_fraction(i), Ion13_fraction(i)
     goto 101
100	continue
     fraction_num=i-1

    do k=1, fraction_num, 1
! read the Ion_4 calculated data
	i=0
111	i=i+1
     read(11,*,end=110) lowenergy(i),&
     wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
     goto 30
     goto 111
110	continue
      Ion4_num=i-1


! calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
      do j=2,trans
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

      do k=1,trans
         wfwhm(k)=fwhmgauss
      enddo

      mm=0
      do 100 i=1,trans
      if(abs(wavelength(i)).gt. emin .and. abs(wavelength(i)).lt.emax)then
		      mm=mm+1
			wavenew(mm)=abs(1239.85/(1239.85/wavelength(i)-delta))
			intensitynew(mm)=abs(gf(i))
            fwhmnew(mm)=wfwhm(i)*2
	        flowk(mm)=lowk(i)
	        fhighk(mm)=highk(i)
	        newfjtp(mm)=fjtp(i)
            newfjt(mm)=fjt(i)
            fwspecen(mm)=wspecen(i)
            flowenergy(mm)=lowenergy(i)            
      write(16,"(f9.3,f9.4)")1239.85/wavelength(i),gf(i)
	  end if
100   continue

         do i=1,mm
       if(flowenergy(i)>fwspecen(i))then
           newwspecen(i)=flowenergy(i)
           ffjtp(mm)=newfjt(i)
	   else
           newwspecen(i)=fwspecen(i)
           ffjtp(mm)=newfjtp(i)
	   end if
           population(i)=(2*ffjtp(i)+1)*exp(-abs(newwspecen(i)-minimum)*0.124/Te_value(k))/(2*lj+1)  !---------Te_value(k)
	    write(15,"(f10.4)") population(i)
         enddo

       do i=1,Exp_num
           wave(i)=Exp_wave(i)  
           crossb_4(i)=0                            !--------------- set initial value
       enddo

      do i=1,Exp_num
	     uu=0
      do j=1,mm
         uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	   &fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	   &fwhmnew(j)**2/4)))
      enddo
         crossb_4(i)=uu*Ion4_fraction(k)            !--------------- set initial value
	  enddo
      do i=1,Exp_num
         if(crossb_4(i).lt.0.00001 ) then
         crossb_4(i)=0                              !--------------- set the small value to Zero
	     endif
	     write(19,*) 1239.85/wave(i), crossb_4(i)   !--------- output the emission of nth ion considering the population and ion-fraction
	enddo






















	i=0
30	i=i+1
     read(11,*,end=40) lowenergy(i),&
     wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
     goto 30
40	continue
      trans=i-1



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
            fwhmnew(mm)=wfwhm(i)*2
	        flowk(mm)=lowk(i)
	        fhighk(mm)=highk(i)
	        newfjtp(mm)=fjtp(i)
            newfjt(mm)=fjt(i)
            fwspecen(mm)=wspecen(i)
            flowenergy(mm)=lowenergy(i)            
        write(16,"(f9.3,f9.4)")1239.85/wavelength(i),gf(i)
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
	     uu=0
      do j=1,mm
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
           if(crossb(i).lt.0.00001 ) then
           crossb(i)=0
	     endif
	     write(19,*) 1239.85/wave(i), crossb(i)
	enddo

      end
