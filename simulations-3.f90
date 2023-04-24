
    program Simulations


	parameter(m=150000,n=2000,delta=-0,dd=100000)
    real wavelength(m),gf(m),specen(m),fwhm(m),wfwhm(m),lj,&
         aversum,averintensity,average,wspecen(m),fwspecen(m),lowk(m),highk(m),&
         wave(n),intensity(n),fwhmnew(m),wavenew(m),newwspecen(m),&
         intensitynew(m),tt,emin,emax,cross(n),ss,flowk(m),fwhmgauss,&
         fhighk(m),fjt(m),fjtp(m),ffjtp(m),ffjt(m),lowenergy(m),&
         flowenergy(m),minimum,crossb(n),population(m),uu,ffwhm(m),newfjt(m),newfjtp(m),&
		 Exp_inten(dd),Exp_wave(dd),Theo_wave(10,3000),Theo_inten(10,3000),&
		 Ion4_fraction(dd),Ion5_fraction(dd),Ion6_fraction(dd),&
         Ion7_fraction(dd),Ion8_fraction(dd),Ion9_fraction(dd),Ion10_fraction(dd),Ion11_fraction(dd),&
	     Ion12_fraction(dd),Ion13_fraction(dd),Te_value(dd),Ne_value(dd),crossb_4(m)
   integer Exp_num,fraction_num,Ion4_num,k,t,mm


	emin=49
!   49 eV=25.3 nm
	emax=250
!   250 eV=4.96 nm
!    Te=25
    fwhmgauss=0.27


open(9,file="EXP.dat",status="unknown")
open(10,file="Ion-fraction.dat",status="unknown") 


  

open(11,file="spectra-Ge4.dat",status="unknown")
open(12,file="spectra-Ge5.dat",status="unknown")
open(13,file="spectra-Ge6.dat",status="unknown")
open(14,file="spectra-Ge7.dat",status="unknown")
open(15,file="spectra-Ge8.dat",status="unknown")
open(16,file="spectra-Ge9.dat",status="unknown")
open(17,file="spectra-Ge10.dat",status="unknown")
open(18,file="spectra-Ge11.dat",status="unknown")
open(19,file="spectra-Ge12.dat",status="unknown")
open(20,file="spectra-Ge13.dat",status="unknown")

!open(11,file="spectra-temp.dat",status="unknown")

open(110,file="cross-P.dat",status="unknown") 

!     read EXP data
	  i=0
91	  i=i+1
      read(9,*,end=90) Exp_wave(i), Exp_inten(i)
      goto 91
90	  continue
      Exp_num=i-1
!	  write(6,*) Exp_num

!     read the ion fraction
	  i=0
101	  i=i+1
      read(10,*,end=100) Te_value(i),Ne_value(i),Ion4_fraction(i),Ion5_fraction(i),Ion6_fraction(i),&
            Ion7_fraction(i),Ion8_fraction(i),Ion9_fraction(i),Ion10_fraction(i),Ion11_fraction(i),&
	        Ion12_fraction(i),Ion13_fraction(i)
      goto 101
100	  continue
            fraction_num=i-1
      write(6,*) fraction_num

!     do i=1,fraction_num,1
!	  write(110,"(f9.3, e14.5,10f12.6)") Te_value(i),Ne_value(i),Ion4_fraction(i),Ion5_fraction(i),Ion6_fraction(i),&
!     Ion7_fraction(i),Ion8_fraction(i),Ion9_fraction(i),Ion10_fraction(i),Ion11_fraction(i),&
!	  Ion12_fraction(i), Ion13_fraction(i)
!     enddo


      do k=60, 60,1
!	  write(6,"(f9.3,f12.6)") Te_value(k),Ion4_fraction(k)
!     read the Ion_4 calculated data
	  i=0
111	  i=i+1
      read(11,*,end=110) lowenergy(i),&
           wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 111
110	  continue
!-----number of Ion4 file 
      Ion4_num=i-1                               
!      write(6,*) Ion4_num

!	  do i=1,Ion4_num
!	  write(110,"(8f14.3)") lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
!     enddo


!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!-----number of Ion4 file  
      do j=2,Ion4_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo
!     write(6,*) minimum,lj

      do t=1,Ion4_num                            
         wfwhm(t)=fwhmgauss
!        write(6,*) t,wfwhm(t)
      enddo

      mm=0
!-----number of Ion4 file  
      do 50 i=1,Ion4_num                         
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
!     write(110,"(i5,2f9.4)") mm, 1239.85/wavelength(mm),gf(mm)
	  end if
50    continue

       do i=1,mm
           if(flowenergy(i) .gt. fwspecen(i))then
           newwspecen(i)=flowenergy(i)
           ffjtp(mm)=newfjt(i)
	       else
           newwspecen(i)=fwspecen(i)
           ffjtp(mm)=newfjtp(i)
	       end if
!---------Te_value(k)
           population(i)=(2*ffjtp(i)+1)*exp(-abs(newwspecen(i)-minimum)*0.124/Te_value(k))/(2*lj+1)  
     	   write(110,"(i4,2f10.4)")i,population(i),Te_value(k)
       enddo

       do i=1,Exp_num
           wave(i)=Exp_wave(i)  
!----------set initial value
           crossb_4(i)=0                                 
!		   write(110,*) 1239.85/wave(i), crossb_4(i)
       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
!             write(110,*) i,uu
           enddo
           crossb_4(i)=uu
!----------obtain the calculated value
!          crossb_4(i)=uu*Ion4_fraction(k)               
!		   write(110,*) 1239.85/wave(i),uu
           write(6,*) mm
	   enddo

       do i=1,Exp_num
!          if(crossb_4(i).lt.0.00001 ) then
!----------set the small value to Zero
!          crossb_4(i)=0                                 
!	       endif
!--------- output the emission of nth ion considering the population and ion-fraction
!	       write(110,*) i,1239.85/wave(i), crossb_4(i)   
	   enddo

     enddo
     end
