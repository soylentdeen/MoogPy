      real*8 function expint(x, n)
      INTEGER n, MAXIT
      REAL*8 x, EPS, FPMIN, EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)

      INTEGER i, ii, nm1
      REAL a, b, c, d, del, fact, h, psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
          write(*,*) 'bad arguments in expint', n, x
      else if(n.eq.0)then            ! Special Case
          expint=exp(-x)/x
      else if(x.eq.0.)then           ! Another special case
          expint=1./nm1
      else if(x.gt.1.) then          ! Lentz's algorithm
          b=x+n
          c=1./FPMIN
          d=1./b
          h=d
          do i=1,MAXIT
              a=-1*(nm1+1)
              b=b+2
              d=1./(a*d+b)
              c=b+a/c
              del=c*d
              h=h*del
              if(abs(del-1.).lt.EPS)then
                  expint=h*exp(-x)
                  return
              endif
          enddo
          write(*,*) 'Continued Fraction failed in expint'
      else
          if(nm1.ne.0)then
              expint=1./nm1
          else
              expint=-log(x)-EULER
          endif
          fact=1.
          do i=1,MAXIT
              fact=-fact*x/i
              if(i.ne.nm1)then
                  del=-fact/(i-nm1)
              else
                  psi=-EULER
                  do ii=1,nm1
                      psi=psi+1./ii
                  enddo
              endif
              expint=expint+del
              if(abs(del).lt.abs(expint)*EPS) return
          enddo
          write(*,*) 'series failed in expint'
      endif
      return
      END
