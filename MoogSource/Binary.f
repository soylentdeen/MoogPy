
      subroutine binary
c******************************************************************************
c     This program synthesizes a section of a spectroscopic binary star
c     spectrum, using one model atmosphere for the primary and one model
c     for the secondary, and then compares it to an observed spectrum.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Equivs.com'
      include 'Multistar.com'
      integer ncall


c*****examine the parameter file
      begin = 0
      ncall = 1
1     call params
      linprintopt = linprintalt
      if (begin .eq. 0) then
         if (numpecatom .gt. 0) then
            do i=3,95
               binpec(syncount,i) = pec(i)
               do j=1,numatomsyn
                  binpecabund(syncount,i,j) = pecabund(i,j)
               enddo
            enddo
         endif
      else
         if (numpecatom .gt. 0) then
            do i=3,95
               pec(i) = binpec(syncount,i)
               do j=1,numatomsyn
                  pecabund(i,j) = binpecabund(syncount,i,j)
               enddo
            enddo
         endif
      endif


c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra
      istat = ivcleof(5,1)
      nf1out = 20     
      lscreen = 4
      array = 'STANDARD OUTPUT'
      nchars = 15
      call infile ('output ',nf1out,'formatted  ',0,nchars,
     .             f1out,lscreen)
      nf2out = 21               
      lscreen = lscreen + 2
      array = 'RAW SYNTHESIS OUTPUT'
      nchars = 20
      call infile ('output ',nf2out,'formatted  ',0,nchars,
     .             f2out,lscreen)
      if (syncount .eq. 1) then
         f7out = f2out
      else
         f8out = f2out
      endif
      if (iraf .ne. 0) then
         nf4out = 23               
         lscreen = lscreen + 2
         array = 'IRAF ("rtext") OUTPUT'
         nchars = 24
         call infile ('output ',nf4out,'formatted  ',0,nchars,
     .                f4out,lscreen)
      endif


c*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel


c*****open the line list file and the strong line list file
      nflines = 31
      lscreen = lscreen + 2
      array = 'THE LINE LIST'
      nchars = 13
      call infile ('input  ',nflines,'formatted  ',0,nchars,
     .              flines,lscreen)
      if (dostrong .gt. 0) then
         nfslines = 32
         lscreen = lscreen + 2
         array = 'THE STRONG LINE LIST'
         nchars = 20
         call infile ('input  ',nfslines,'formatted  ',0,nchars,
     .                 fslines,lscreen)
      endif
      

c*****do the syntheses
      if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isorun = 1
         nlines = 0
         mode = 3
         call inlines (1)
         call eqlib
         call nearly (1)
         call synspec
      else
         do n=1,numatomsyn
            isynth = n
            isorun = 1
            start = oldstart
            sstop = oldstop
            mode = 3
            call inlines (1)
            call eqlib
            call nearly (1)
            call synspec
            linprintopt = 0
         enddo
      endif
      if (syncount .eq. 1) then
         fluxprimary = flux
      else
         fluxsecondary = flux
      endif
         

c*****finish the syntheses
      call finish (1)
      istat = ivcleof(4,1)
      if (control .ne. 'gridend') go to 1


c*****combine the synthetic spectra for plotting


c*****now plot the spectrum, maybe iterating abundances, and end the program
      call finish (0)

      end






