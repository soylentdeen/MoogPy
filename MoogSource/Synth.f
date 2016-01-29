
      subroutine synth
c******************************************************************************
c     This program synthesizes a section of spectrum and compares it
c     to an observation file.
c******************************************************************************

      implicit real*8 (a-h,o-z)
      include 'Atmos.com'
      include 'Factor.com'
      include 'Mol.com'
      include 'Linex.com'
      include 'Pstuff.com'
      include 'Dummy.com'
      include 'Stokes.com'

cf2py intent(callback, hide) tennisball
      external tennisball

c*****examine the parameter file
      call params


c*****open the files for: standard output, raw spectrum depths, smoothed 
c     spectra, and (if desired) IRAF-style smoothed spectra

c     Commenting out the output files.
c
c      nf1out = 20     
c      lscreen = 4
c      array = 'STANDARD OUTPUT'
c      nchars = 15
c      call infile ('output ',nf1out,'formatted  ',0,nchars,
c     .             f1out,lscreen)
c      nf2out = 21               
c      lscreen = lscreen + 2
c      array = 'RAW SYNTHESIS OUTPUT'
c      nchars = 20
c      call infile ('output ',nf2out,'formatted  ',0,nchars,
c     .             f2out,lscreen)
c      if (plotopt .ne. 0) then
c         nf3out = 22               
c         lscreen = lscreen + 2
c         array = 'SMOOTHED SYNTHESES OUTPUT'
c         nchars = 25
c         call infile ('output ',nf3out,'formatted  ',0,nchars,
c     .                f3out,lscreen)
c         if (f5out .ne. 'optional_output_file') then
c            nf5out = 26
c            lscreen = lscreen + 2
c            array = 'POSTSCRIPT PLOT OUTPUT'
c            nchars = 22
c            call infile ('output ',nf5out,'formatted  ',0,nchars,
c     .                   f5out,lscreen)
c         endif
c      endif
c      if (iraf .ne. 0) then
c         nf4out = 23               
c         lscreen = lscreen + 2
c         array = 'IRAF ("rtext") OUTPUT'
c         nchars = 24
c         call infile ('output ',nf4out,'formatted  ',0,nchars,
c     .                f4out,lscreen)
c      endif


c*****open and read the model atmosphere file
      nfmodel = 30
      lscreen = lscreen + 2
      array = 'THE MODEL ATMOSPHERE'
      nchars = 20
      call infile ('input  ',nfmodel,'formatted  ',0,nchars,
     .             fmodel,lscreen)
      call inmodel
c      write (*,*) "Read in the Model"


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
      ncall = 1
      call tennisball
10    if (numpecatom .eq. 0 .or. numatomsyn .eq. 0) then
         isynth = 1
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
            isorun = n
            start = oldstart
            sstop = oldstop
            mode = 3
            call inlines (1)
              molopt = 2
            call eqlib
            call nearly (1)
            call synspec
            linprintopt = 0
         enddo
      endif
         

c*****now plot the spectrum


c*****if the syntheses need to be redone: first rewind the output files,
c     then close/reopen line list(s), then rewrite model atmosphere output
c      if (choice .eq. 'n') then
c         call chabund
cc         if (choice .eq. 'x') go to 20
c         rewind nf1out
c         rewind nf2out
c         if (nflines .ne. 0) then
c            close (unit=nflines)
c            open (unit=nflines,file=flines,access='sequential',
c     .            form='formatted',blank='null',status='old',
c     .            iostat=jstat,err=10)
c         endif
c         if (nfslines .ne. 0) then
c            close (unit=nfslines)
c            open (unit=nfslines,file=fslines,access='sequential',
c     .            form='formatted',blank='null',status='old',
c     .            iostat=jstat,err=10)
c         endif
c         if (plotopt .ne. 0) then
c            rewind nf3out
c         endif
c         write (nf1out,1002) modtype
c         if (modprintopt .ge. 1) then
c            if (modtype .eq. 'begn      ' .or.
c     .          modtype .eq. 'BEGN      ') write (nf1out,1003)
c            write (nf1out,1102) moditle
c            do i=1,ntau
c               dummy1(i) = dlog10(pgas(i))
c               dummy2(i) = dlog10(ne(i)*1.38054d-16*t(i))
c            enddo
c            write (nf1out,1103) wavref,(i,xref(i),tauref(i),t(i),
c     .                          dummy1(i), pgas(i),dummy2(i),ne(i),
c     .                          vturb(i),i=1,ntau)
c            write (nf1out,1104)
c            do i=1,95
c               dummy1(i) = dlog10(xabund(i)) + 12.0
c            enddo
c            write (nf1out,1105) (names(i),i,dummy1(i),i=1,95)
c            write (nf1out,1106) modprintopt, molopt, linprintopt, 
c     .                          fluxintopt
c            write (nf1out,1107) (kapref(i),i=1,ntau)
c         endif
c         linprintopt = linprintalt
c         ncall = 2
c         choice = '1'
c         go to 10


c*****otherwise end the code gracefully
c      else
c         call finish (0)
c      endif

      call finish (0)

c*****format statements
1002  format (13('-'),'MOOG OUTPUT FILE',10('-'),
     .        '(MOOG version from 23/04/07)',13('-')//
     .        'THE MODEL TYPE: ',a10)
1003  format ('   The Rosseland opacities and optical depths have ',
     .        'been read in')
1102  format (/'MODEL ATMOSPHERE HEADER:'/a80/)
1103  format ('INPUT ATMOSPHERE QUANTITIES',10x,
     .        '(reference wavelength =',f10.2,')'/3x,'i',2x,'xref',3x,
     .        'tauref',7x,'T',6x,'logPg',4x,'Pgas',6x,'logPe',
     .        5x,'Ne',9x,'Vturb'/
     .        (i4,0pf6.2,1pd11.4,0pf9.1,f8.3,1pd11.4,0pf8.3,
     .        1pd11.4,d11.2))
1104  format (/'INPUT ABUNDANCES: (log10 number densities, log H=12)'/
     .       '      Default solar abundances: Anders and Grevesse 1989')
1105  format (5(3x,a2,'(',i2,')=',f5.2))
1106  format (/'OPTIONS: atmosphere = ',i1,5x,'molecules  = ',i1/
     .        '         lines      = ',i1,5x,'flux/int   = ',i1)
1107  format (/'KAPREF ARRAY:'/(6(1pd12.4)))



      end 





