 
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
	     Ion12_fraction(dd),Ion13_fraction(dd),Te_value(dd),Ne_value(dd),crossb_4(m),crossb_5(m),&
         crossb_6(m),crossb_7(m),crossb_8(m),crossb_9(m),crossb_10(m),crossb_11(m),crossb_12(m),&
		 crossb_13(m)
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

open(999,file="cross-P.dat",status="unknown") 

!     read EXP data
	  i=0
91	  i=i+1
!-----Exp_wave(i):  nm as unite
      read(9,*,end=90) Exp_wave(i), Exp_inten(i)
!	  write(999,*)  Exp_wave(i), Exp_inten(i)
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

!      write(6,*) fraction_num

!     do i=1,fraction_num,1
!	  write(999,"(f9.3, e14.5,10f12.6)") Te_value(i),Ne_value(i),Ion4_fraction(i),Ion5_fraction(i),Ion6_fraction(i),&
!     Ion7_fraction(i),Ion8_fraction(i),Ion9_fraction(i),Ion10_fraction(i),Ion11_fraction(i),&
!	  Ion12_fraction(i), Ion13_fraction(i)
!     enddo

      do k=60,60,1
!	  write(6,"(f9.3,f12.6)") Te_value(k),Ion4_fraction(k)

!####################################################################   Ion_4    ##########################

!     read the Ion_4 calculated data
	  i=0
111	  i=i+1
      read(11,*,end=110) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
!      write(999,"(8f9.4)") lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 111
110	  continue
!-----number of Ion4 file 
      Ion4_num=i-1 
!      write(6,*) Ion4_num


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
      do i=1,Ion4_num                         
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
	        end if
      enddo
!      write(6,*) mm


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
!     	   write(999,"(i4,2f10.4)")i,population(i),Te_value(k)
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------set initial value
           crossb_4(i)=0                                 
!		   write(999,*) i,wave(i), crossb_4(i)
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
          crossb_4(i)=uu*Ion4_fraction(k)               
!		   write(999,*) 1239.85/wave(i),uu
!           write(6,*) mm
	   enddo

       do i=1,Exp_num
          if(crossb_4(i).lt.0.00001 ) then
!----------set the small value to Zero
          crossb_4(i)=0                                 
	       endif
!--------- output the emission of nth ion considering the population and ion-fraction
!	       write(999,*) 1239.85/wave(i), crossb_4(i)   
	   enddo


!####################################################################   Ion_4_end    ##########################

!####################################################################   Ion_5_start  ##########################

	  i=0
121	  i=i+1
      read(12,*,end=120) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 121
120	  continue
!-------->
      Ion5_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion5_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion5_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion5_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_5(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->---------->
          crossb_5(i)=uu*Ion5_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_5(i).lt.0.00001 ) then
!---------------->
          crossb_5(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_5(i)   
	   enddo
!####################################################################   Ion_5_end    #########################

!####################################################################   Ion_6_start  ##########################  13-Ge6+

	  i=0
131	  i=i+1
      read(13,*,end=130) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 131
130	  continue
!-------->
      Ion6_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion6_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion6_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion6_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_6(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->---------->
          crossb_6(i)=uu*Ion6_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_6(i).lt.0.00001 ) then
!---------------->
          crossb_6(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_6(i)   
	   enddo
!####################################################################   Ion_6_end    #########################

!####################################################################   Ion_7_start  ##########################  14-Ge7+

	  i=0
141	  i=i+1
      read(14,*,end=140) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 141
140	  continue
!-------->
      Ion7_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion7_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion7_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion7_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_7(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->---------->
          crossb_7(i)=uu*Ion7_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_7(i).lt.0.00001 ) then
!---------------->
          crossb_7(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_7(i)   
	   enddo
!####################################################################   Ion_7_end    #########################

!####################################################################   Ion_8_start  ##########################  15-Ge8+

	  i=0
151	  i=i+1
      read(15,*,end=150) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 151
150	  continue
!-------->
      Ion8_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion8_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion8_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion8_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_8(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->---------->
          crossb_8(i)=uu*Ion8_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_8(i).lt.0.00001 ) then
!---------------->
          crossb_8(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_8(i)   
	   enddo
!####################################################################   Ion_8_end    #########################


!####################################################################   Ion_9_start  ##########################  16-Ge9+

	  i=0
161	  i=i+1
      read(16,*,end=160) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 161
160	  continue
!-------->
      Ion9_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion9_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion9_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion9_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_9(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->---------->
          crossb_9(i)=uu*Ion9_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_9(i).lt.0.00001 ) then
!---------------->
          crossb_9(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_9(i)   
	   enddo
!####################################################################   Ion_9_end    #########################

!####################################################################   Ion_10_start  ##########################  17-Ge10+

	  i=0
171	  i=i+1
      read(17,*,end=170) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 171
170	  continue
!-------->
      Ion10_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion10_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion10_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion10_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_10(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->----------->
          crossb_10(i)=uu*Ion10_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_10(i).lt.0.00001 ) then
!---------------->
          crossb_10(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_10(i)   
	   enddo
!####################################################################   Ion_10_end    #########################

!####################################################################   Ion_11_start  ##########################  18-Ge11+

	  i=0
181	  i=i+1
      read(18,*,end=180) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 181
180	  continue
!-------->
      Ion11_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion11_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion11_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion11_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_11(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->----------->
          crossb_11(i)=uu*Ion11_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_11(i).lt.0.00001 ) then
!---------------->
          crossb_11(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_11(i)   
	   enddo
!####################################################################   Ion_11_end    #########################

!####################################################################   Ion_12_start  ##########################  19-Ge12+

	  i=0
191	  i=i+1
      read(19,*,end=190) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 191
190	  continue
!-------->
      Ion12_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion12_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion12_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion12_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_12(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->----------->
          crossb_12(i)=uu*Ion12_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_12(i).lt.0.00001 ) then
!---------------->
          crossb_12(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_12(i)   
	   enddo
!####################################################################   Ion_12_end    #########################

!####################################################################   Ion_13_start  ##########################  20-Ge13+

	  i=0
201	  i=i+1
      read(20,*,end=200) lowenergy(i),wspecen(i),wavelength(i),gf(i),lowk(i),highk(i),fjt(i),fjtp(i)
      goto 201
200	  continue
!-------->
      Ion13_num=i-1 

!     calculate the simulated spcetra of Ion4 ion
      minimum=lowenergy(1)
!--------------->  
      do j=2,Ion13_num                            
	  if (minimum .gt. lowenergy(j)) then
         minimum=lowenergy(j)
	     lj=fjt(i)
	  endif
      enddo

!--------------->
      do t=1,Ion13_num                            
         wfwhm(t)=fwhmgauss
      enddo

      mm=0
!---------------> 
      do i=1,Ion13_num                         
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
	        end if
      enddo

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
       enddo

       do i=1,Exp_num
           wave(i)=1239.85/Exp_wave(i)  
!----------------->
           crossb_13(i)=0                                 

       enddo

       do i=1,Exp_num
	       uu=0
           do j=1,mm
           uu=uu+(intensitynew(j)*population(j)/(2*ffjtp(j)+1))*&
     	      fwhmnew(j)/(2*3.1415928*(((wavenew(j)-wave(i))**2+&
     	      fwhmnew(j)**2/4)))
           enddo
!-------------- ->----------->
          crossb_13(i)=uu*Ion13_fraction(k)               
	   enddo

       do i=1,Exp_num
!------------------->
          if(crossb_13(i).lt.0.00001 ) then
!---------------->
          crossb_12(i)=0                                 
	       endif
!----------------------------------------------->
	       write(999,*) 1239.85/wave(i), crossb_13(i)   
	   enddo
!####################################################################   Ion_13_end    #########################

!     SUM
      do i=1,,Exp_num
      sim_spec(i)=0
      enddo

      do i=1,Exp_num
         sim_spec(i)=crossb_4(i)+crossb_5(i)+crossb_6(i)+crossb_7(i)+crossb_8(i)+crossb_9(i)+crossb_10(i)+crossb_11(i)+crossb_12(i)+crossb_13(i)
      enddo


     enddo


     end
