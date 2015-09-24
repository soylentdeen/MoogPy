
      subroutine synspec ()
c******************************************************************************
c     This routine does synthetic spectra                                
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Linex.com'
      include 'Factor.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      real*8 dd(5000)

cf2py intent(callback, hide) recorder
      external recorder

c*****initialize the synthesis
c      write (nf1out,1101)
c      write (nf2out,1002) moditle(1:73)
c      if (iunits .eq. 1) then
c         write (nf2out,1103) oldstart,oldstop,oldstep,olddelta
c      else
c         write (nf2out,1102) start,sstop,step,delta
c      endif
c      if (iraf .eq. 1) then
c         npoint = (sstop-start)/step
c         write (nf4out,1104) npoint,wave,wave,step,step
c         write (nf4out,1105)
c         write (nf4out,1106) moditle
c         do j=1,93
c            if (pec(j) .gt. 0 ) then
c               dummy1(j) = dlog10(xabund(j)) + 12.0
c               write (nf4out,1107) names(j),dummy1(j)
c            endif
c         enddo
c         write (nf4out,1108) vturb(1)
c         write (nf4out,1109)
c      endif
      n = 1           
      num = 0
      nsteps = 1
      if (mode .ne. 4) then 
         lim1line = 0
         lim2line = 0
         lim1obs = 0
         lim2obs = 0
         lim1 = 0
         lim2 = 0
      endif


c*****calculate continuum quantities at the spectrum wavelength
      wave = start
      wavl = 0.
30    if (dabs(wave-wavl)/wave .ge. 0.001) then
         wavl = wave   
         call opacit (2,wave)    
c         if (modprintopt .ge. 2) 
c     .       write (nf1out,1001) wave,(kaplam(i),i=1,ntau)
         call cdcalc (1)  
         first = 0.4343*cd(1)
         flux = rinteg(xref,cd,dummy1,ntau,first)
c         if (iunits .eq. 1) then
c            write (nf1out,1003) 1.d-4*wave,flux
c         else
c            write (nf1out,1004) wave,flux
c         endif
      endif


c*****find the appropriate set of lines for this wavelength, reading 
c     in a new set if needed
      if (mode .eq. 3) then
20       call linlimit
         if (lim2line .lt. 0) then
            call inlines (2)
            call nearly (1)
            go to 20
         endif
         lim1 = lim1line
         lim2 = lim2line
      endif


c*****compute a spectrum depth at this point
      call taukap   
      call cdcalc (2)
      first = 0.4343*cd(1)
      d(n) = rinteg(xref,cd,dummy1,ntau,first)
      call recorder(wave, d(n))
c      if (mod(n,10) .eq. 0) then
c         if (iraf .eq. 1) then
c            do j=1,10
c               dd(num+j) = 1. - d(num+j)
c            enddo
c            write (nf4out,1110) (dd(num+j),j=1,10)
c         endif
c         if (iunits .eq. 1) then
c            wave3 = 1.d-4*(wave - 9.0*step)
c            write (nf1out,1112) wave3,(d(num+j),j=1,10)
c         else
c            wave3 = wave - 9.0*step
c            write (nf1out,1111) wave3,(d(num+j),j=1,10)
c         endif
c         if (nf2out .gt. 0) write (nf2out,1110) (d(num+j),j=1,10)
c         num = num + 10
c      endif


c*****step in wavelength and try again 
      wave = oldstart + step*nsteps
c      write (*,*) wave, step, nsteps, sstop
      if (wave .le. sstop) then
         n = n + 1        
         nsteps = nsteps + 1
         if (n .gt. 5000) then
            n = 1                                      
            num = 0
         endif
         go to 30                   


c*****finish the synthesis
      else
c         nn = mod(n,10)
c         if (nn .ne. 0) then
c            if (iraf .eq. 1) then
c               do j=1,nn
c                  dd(num+j) = 1. - d(num+j)
c               enddo
c               write (nf4out,1110) (dd(num+j),j=1,nn)
c            endif
c            if (iunits .eq. 1) then
c               wave3 = 1.d-4*(wave - 9.0*step)
c               write (nf1out,1112) wave3,(d(num+j),j=1,nn)
c            else
c               wave3 = wave - 9.0*step
c               write (nf1out,1111) wave3,(d(num+j),j=1,nn)
c            endif
c            if (nf2out .gt. 0) write (nf2out,1110) (d(num+j),j=1,nn)
c         endif
c         if (iunits .eq. 1) then
c            write (nf1out,1113) 1.d-4*wave
c         else
c            write (nf1out,1114) wave
c         endif
c         do j =1,30
c             write (*,*) "blah - ", j
c         enddo
         return 
      endif


c*****format statements
1001  format ('  kaplam from 1 to ntau at wavelength',f10.2/
     .        (6(1pd12.4)))
1002  format ('MODEL: ',a73)
1003  format ('AT WAVELENGTH/FREQUENCY =',f11.7,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1004  format ('AT WAVELENGTH/FREQUENCY =',f11.3,
     .        '  CONTINUUM FLUX/INTENSITY =',1p,d12.5)
1101  format (/'SPECTRUM DEPTHS')
1102  format (4f11.3)
1103  format (4f10.7)
1104  format ('SIMPLE  =    t'/'NAXIS   =     1'/'NAXIS1  = ',i10,/
     .        'W0      =',f10.4/'CRVAL1  =',f10.4/'WPC     =',f10.4/
     .        'CDELT1  =',f10.4)
1105  format (16HORIGIN  = 'moog'/21HDATA-TYP= 'synthetic'/
     .        18HCTYPE1  = 'lambda'/21HCUNIT1  = 'angstroms')
1106  format (11HTITLE   = ',A65,1H')
1107  format ('ATOM    = ',1H',7x,a2,1H',/,'ABUND   = ',f10.2)
1108  format ('VTURB   = ',d10.4,'     /  cm/sec  ')
1109  format ('END')
1110  format (10f7.4)
1111  format (f10.3,': depths=',10f6.3)
1112  format (f10.7,': depths=',10f6.3)
1113  format ('FINAL WAVELENGTH/FREQUENCY =',f10.7/)
1114  format ('FINAL WAVELENGTH/FREQUENCY =',f10.3/)


      end                                




