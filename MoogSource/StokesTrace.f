
      subroutine stokestrace
c******************************************************************************
c     This program can synthesize the full stokes parameters for multiple
c     sections of spectra for multiple input model atmospheres
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Stokes.com'
      include 'Angles.com'

      zeros(1,1) = 0.0
      zeros(1,2) = 0.0
      zeros(1,3) = 0.0
      zeros(2,1) = 0.0
      zeros(2,2) = 0.0
      zeros(2,3) = 0.0
      zeros(3,1) = 0.0
      zeros(3,2) = 0.0
      zeros(3,3) = 0.0

c*****examine the parameter file
      call params
      linfileopt = 2
      linprintopt = linprintalt
      
c*****open and read the model atmosphere file
      nfmodel = 30
c      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel

c*****open the line list file and the strong line list file
      nflines = 31
c      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
c         lscreen = lscreen + 2
         array = 'THE STRONG LINE LIST'
         nchars = 20
         call infile ('input  ',nfslines,'formatted  ',0,nchars,
     .                 fslines,lscreen)
      endif

      B_sph(1) = 1.0
      B_sph(2) = 0.0
      B_sph(3) = 0.0
      
      B_xyz(1) = 0.0
      B_xyz(2) = 0.0
      B_xyz(3) = 0.0

      start = oldstart
      isynth = 1
      call wavegrid
c*****Read in the line list and calculate the equilibria
      call inlines (1)
      call eqlib
      call nearly (1)

c***** Calculate zdepth, the physical depth scale
      call spl_def(ntau, xref, kapref, kref_knots, n_kref_knots,
     .             kref_coeffs)
      nz = 1
      dtau = 0.05
      do ztau = xref(1), xref(ntau), dtau
          if (nz.eq.1) then
              zdepth(nz) = 0.0
          else
              h1=1.0/spl_ev(kref_knots, n_kref_knots, kref_coeffs,
     .                      ztau)
              h2=1.0/spl_ev(kref_knots, n_kref_knots, kref_coeffs,
     .                      ztau-dtau)
              dt = (10.0**ztau-10.0**(ztau-dtau))
              zdepth(nz) = zdepth(nz-1)+(h1+h2)/2.0*dt
          endif
          taus(nz) = ztau
          nz=nz+1
      enddo
      taus(nz) = xref(ntau)
      h2=1.0/kapref(ntau)
      dt = 10.0**taus(nz)-10.0**taus(nz-1)
      zdepth(nz) = zdepth(nz-1)+(h2+h1)/2.0*dt
      call spl_def(nz, taus, zdepth, z_knots, n_z_knots, z_coeffs)

      curr_strong = 1
      curr_weak = 1
c*****Perform the Synthesis
      wavl = 0.
      mode = 3
      wave = dipstickwave
30    if (dabs(wave-wavl)/wave .ge. 0.001) then
         wavl = wave
         call opacit (2,wave)
      endif
20    call linlimit
      if (lim2line .lt. 0) then
          call inlines(2)
          call nearly (1)
          go to 20
      endif
      lim1 = lim1line
      lim2 = lim2line
      call calcopacities
      call stokesdipstick(dble(0.0), dble(0.0), dble(1.0))
      

c*****finish
      call finish (0)

c1001  format (a80)
c12345 format (f10.4,5e15.5)
12346 format (i5,8e16.5)
12347 format (2i5,2e16.5)
12348 format (i5)
12349 format (i5,2e16.5)
6520  format (f10.4)
6521  format (e16.8)
      end 
