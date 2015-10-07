
      subroutine synstokes
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
      integer n_cells, icell, curr_strong, curr_weak
      real*8 az_start, az_stop, daz, az, long, dlong
      real*8 ring_area, cell_area, cell_a, phi_ang, chi_ang, mu

cf2py intent(callback, hide) stokesrecorder
      external stokesrecorder

cf2py intent(callback, hide) beachball
      external beachball

cf2py intent(callback, hide) diskoball
      external diskoball

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
      
c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra
c      nf1out = 20     
c      lscreen = 4
c      array = 'STANDARD OUTPUT'
c      nchars = 15
c      call infile ('output ',nf1out,'formatted  ',0,nchars,
c     .             f1out,lscreen)
c      nf2out = 21               
cc      lscreen = lscreen + 2
c      array = 'RAW SYNTHESIS OUTPUT'
c      nchars = 20
c      call infile ('output ',nf2out,'formatted  ',0,nchars,
c     .             f2out,lscreen)
c
c*****Open the output files for the Stokes Output
c      nfAngles = 60               
cc      lscreen = lscreen + 2
c      array = 'ANGLES OUTPUT'
c      nchars = 20
c      call infile ('output ',nfAngles,'formatted  ',0,nchars,
c     .             fAngles,lscreen)
c      nfStokesI = 61               
cc      lscreen = lscreen + 2
c      array = 'Stokes I OUTPUT'
c      nchars = 20
c      call infile ('output ',nfStokesI,'formatted  ',0,nchars,
c     .             fStokesI,lscreen)
c      nfStokesQ = 62               
cc      lscreen = lscreen + 2
c      array = 'Stokes Q OUTPUT'
c      nchars = 20
c      call infile ('output ',nfStokesQ,'formatted  ',0,nchars,
c     .             fStokesQ,lscreen)
c      nfStokesU = 63               
cc      lscreen = lscreen + 2
c      array = 'Stokes U OUTPUT'
c      nchars = 20
c      call infile ('output ',nfStokesU,'formatted  ',0,nchars,
c     .             fStokesU,lscreen)
c      nfStokesV = 64               
cc      lscreen = lscreen + 2
c      array = 'Stokes V OUTPUT'
c      nchars = 20
c      call infile ('output ',nfStokesV,'formatted  ',0,nchars,
c     .             fStokesV,lscreen)
c      nfContinuum = 65               
cc      lscreen = lscreen + 2
c      array = 'Continuum OUTPUT'
c      nchars = 20
c      call infile ('output ',nfContinuum,'formatted  ',0,nchars,
c     .             fContinuum,lscreen)
c
c*****open and read the model atmosphere file
      nfmodel = 30
c      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel

      write(*,*) fmodel


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

      if (diskflag .eq. 0) then
c          write (nfAngles, 12347) ncells, nrings, inclination,
c     .                               position_angle
          cell_area = 4.0*3.14159262/ncells
          radtodeg = 180.0/3.1459262
          icell = 1
          do i=1,nrings
             az_start = (i-1)*3.14159262/nrings - 3.14159262/2.0
             az_stop = i*3.14159262/nrings - 3.14159262/2.0
             az = (az_start + az_stop)/2.0
             daz= sin(az_stop) - sin(az_start)
             ring_area = 2.0*3.14159262 * daz ! total sterrad in circular ring
             n_cells = nint(ring_area/cell_area) ! # of cells in ring
             cell_a = ring_area/float(n_cells)
             dlong = 2.0*3.14159262/n_cells
             do j=1,n_cells
                long = -3.14159262+(j-0.5)*dlong
                call calcGeom(az, long, phi_ang, chi_ang, mu)
                if (mu .ge. 0.001) THEN
                   phi_angle(icell) = phi_ang
                   chi_angle(icell) = chi_ang
                   azimuth(icell) = az
                   longitude(icell) = long
                   mus(icell) = mu
c                   write (nfAngles, 12346)icell, az, az_start, az_stop,
c     .                    long, dlong, phi_ang, chi_ang, mu
                   icell = icell + 1
                endif
             enddo
          enddo
          icell = icell-1
          call diskoball
      else
c*****  Now need to make the simple Stokes I, V disk sampling routine.
          phi_angle(1) = 0.270640041
          chi_angle(1) = 0.0
          mus(1) = 0.93659998
          phi_angle(2) = 0.481286
          chi_angle(2) = 0.0
          mus(2) = 0.88634
          phi_angle(3) = 0.640495
          chi_angle(3) = 0.0
          mus(3) = 0.8018
          phi_angle(4) = 0.785389
          chi_angle(4) = 0.0
          mus(4) = 0.7071
          phi_angle(5) = 0.929793
          chi_angle(5) = 0.0
          mus(5) = 0.5979998
          phi_angle(6) = 1.089532
          chi_angle(6) = 0.0
          mus(6) = 0.4629
          phi_angle(7) = 1.300206
          chi_angle(7) = 0.0
          mus(7) = 0.2673
          icell = 7
          
c          write(nfAngles, 12348) 7
c          do i=1,7
c              write (nfAngles, 12349) i, phi_angle(i), mus(i)
c          enddo
          call beachball
      endif

      start = oldstart
      isynth = 1
      call wavegrid
c      write (*,*) "Number of Strong Lines ", ns_lines
c      write (*,*) "Number of weak Lines ", nw_lines
c      read (*,*)
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
      wave = oldstart
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
      if (testflag .eq. 1) then
         call traceStokes(dble(0.69813), dble(0.0), dble(1.0))
         call stokesrecorder(1, wave, Stokes, continuum)
c         write (*,*) wave, Stokes, continuum
      else
         do i = 1, icell
            call traceStokes(phi_angle(i), chi_angle(i), mus(i))
            call stokesrecorder(i, wave, Stokes, continuum)
         enddo
      endif
      
c****      Calculate the distances to the closest strong/weak line
      if (curr_strong .eq. 1) then
          strong_blue_distance = 1000.0
      else
          strong_blue_distance = dabs(wave - strong(curr_strong-1))
      endif
      if (curr_weak .eq. 1) then
          weak_blue_distance = 1000.0
      else
          weak_blue_distance = dabs(wave - weak(curr_weak-1))
      endif

      if (curr_strong .eq. ns_lines) then
          strong_red_distance = 1000.0
      else
          strong_red_distance = dabs(strong(curr_strong) - wave)
      endif
          
      if (curr_weak .eq. nw_lines) then
          weak_red_distance = 1000.0
      else
          weak_red_distance = dabs(weak(curr_weak) - wave)
      endif

65    if ((wave .ge. strong(curr_strong)).and.
     .          (curr_strong .lt. ns_lines)) then
          curr_strong = curr_strong + 1
          goto 65
      endif
66    if ((wave .ge. weak(curr_weak)).and.
     .          (curr_weak .lt. nw_lines)) then
          curr_weak = curr_weak + 1
          goto 66
      endif

c****      Calculate the next wavelength step
      strong_distance = MIN(strong_red_distance,
     .           strong_blue_distance)
      weak_distance = MIN(weak_red_distance, 
     .           weak_blue_distance)

      strong_step = 0.15 - 0.14/(1.0+exp(beta_strong*
     .          (strong_distance/R_strong-1.0)))
      weak_step = 0.15 - 0.14/(1.0+exp(beta_weak*
     .          (weak_distance/R_weak-1.0)))
      
c      write (*,*) wave, strong_step, weak_step
      write (*,*) wave
      wave = wave + MIN(strong_step, weak_step)

      if (wave .le. sstop) then
          go to 30
      endif

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
