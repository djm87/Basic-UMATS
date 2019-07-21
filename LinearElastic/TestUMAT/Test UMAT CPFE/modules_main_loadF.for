c      
c***********************************************************************
c *** dimensions *******************************************************
c***********************************************************************
      module dimensions
      
      INTEGER, PARAMETER :: NTENS=6,NSTATV=700,NPROPS=2
      
      end module dimensions 
c      
c***********************************************************************
c *** miscellaneous_sub ************************************************
c***********************************************************************
      module miscellaneous_sub_main
      !jacobi :
      !eigsrt :
      !det :Calculates determinant of 3x3 matrix      
      !lubksbc :Solves the set of N linear equations
      !ludcmpc :Performs the LU decomposition of a matrix 
      !tmismatch :Calculates scalar difference between tensors
      !tnorm :Calculates norm of a tensor
      !ran2 :Generates random number given seed
      !euler :Calc. transform matrix for each grain 
      !reorient_grain :Calc. rotation matrix from incrmental spin
      !rodrigues :Calc. rotation matrix from incrmental spin
      contains
      
      SUBROUTINE jacobi(a,n,np,d,v,nrot,ier)

      !Copied from VPSC6 by JN, 8-7-07
c **********************************************************************

      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do ip=1,n
        do iq=1,n
          v(ip,iq)=0.
        enddo
        v(ip,ip)=1.
      enddo
      do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
      enddo
      nrot=0
      do i=1,50
        sm=0.
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))

          enddo
        enddo
c
        if(sm.eq.0.)then
        ier=0
        return
        endif
c
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     #        g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)

                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
        enddo
      enddo
c      pause 'too many iterations in jacobi'
c
      ier=1
c
      return
      END SUBROUTINE
      
      FUNCTION tnorm(v,nrows,ncols)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION v(nrows*ncols)

      tnorm=0.0
      do i=1,nrows*ncols
        tnorm=tnorm+v(i)*v(i)
      enddo
      tnorm=sqrt(tnorm)

c      call IEEE('tnorm',iOutExcept)
      RETURN
      END FUNCTION

      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500) !Maximum anticipated value of n.
c     USES pythag
c     Given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes its
c     singular value decomposition, A = U W V^T. The matrix U replaces a on output. The
c     diagonal matrix of singular values W is output as a vector w(1:n). The matrix V (not the
c     transpose V T ) is output as v(1:n,1:n).
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) !,pythag
      g=0.0 !Householder reduction to bidiagonal form.
      scale=0.0
      anorm=0.0
      do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do k=i,m
            scale=scale+abs(a(k,i))
          enddo 
          if(scale.ne.0.0)then
            do  k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
            enddo 
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do  j=l,n
              s=0.0
              do  k=i,m
                s=s+a(k,i)*a(k,j)
              enddo 
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              enddo 
      
      
            enddo 
            do k=i,m
              a(k,i)=scale*a(k,i)
            enddo 
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do k=l,n
            scale=scale+abs(a(i,k))
          enddo 
          if(scale.ne.0.0)then
            do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
            enddo 
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do  k=l,n
              rv1(k)=a(i,k)/h
            enddo 
            do  j=l,m
              s=0.0
              do  k=l,n
               s=s+a(j,k)*a(i,k)
              enddo 
              do  k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              enddo 
            enddo 
            do  k=l,n
              a(i,k)=scale*a(i,k)
            enddo 
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      enddo 
      do  i=n,1,-1 !Accumulation of right-hand transformations.
        if(i.lt.n)then
          if(g.ne.0.0)then
            do  j=l,n !Double division to avoid possible underflow.
              v(j,i)=(a(i,j)/a(i,l))/g
            enddo 
            do  j=l,n
              s=0.0
              do  k=l,n
                s=s+a(i,k)*v(k,j)
              enddo 
              do  k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              enddo 
            enddo 
          endif
          do  j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
          enddo 
        endif
        v(i,i)=1.0
      
      
        g=rv1(i)
        l=i
      enddo 
      do  i=min(m,n),1,-1 !Accumulation of left-hand transformations.
        l=i+1
        g=w(i)
        do  j=l,n
          a(i,j)=0.0
        enddo 
        if(g.ne.0.0)then
          g=1.0/g
          do  j=l,n
            s=0.0
            do  k=l,m
              s=s+a(k,i)*a(k,j)
            enddo 
            f=(s/a(i,i))*g
            do  k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            enddo 
          enddo 
          do  j=i,m
            a(j,i)=a(j,i)*g
          enddo 
        else
          do  j= i,m
            a(j,i)=0.0
          enddo 
        endif
        a(i,i)=a(i,i)+1.0
      enddo 
      do  k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
        do its=1,30 
          do  l=k,1,-1 !Test for splitting.
            nm=l-1 !Note that rv1(1) is always zero.
            if((abs(rv1(l))+anorm).eq.anorm) goto 2
            if((abs(w(nm))+anorm).eq.anorm) goto 1
          enddo 
    1     c=0.0 !Cancellation of rv1(l), if l > 1.
          s=1.0
          do  i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do  j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            enddo 
          enddo 
    2     z=w(k)
          if(l.eq.k)then !Convergence.
            if(z.lt.0.0)then !Singular value is made nonnegative.
              w(k)=-z
              do  j=1,n
                v(j,k)=-v(j,k)
              enddo       
          
          
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l) !Shift from bottom 2-by-2 minor.
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0 !Next QR transformation:
          s=1.0
          do  j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do  jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            enddo 
            z=pythag(f,h)
            w(j)=z !Rotation can be arbitrary if z = 0.
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do  jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            enddo 
          enddo 
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
        enddo 
    3   continue
      enddo 
c      call IEEE('svdcmp',iOutExcept)
      return
      END SUBROUTINE
      
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500) !Maximum anticipated value of n.
c      Solves A X = B for a vector X, where A is specied by the arrays u, w, v as returned by
c      svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices. mp
c      and np are the physical dimensions of a. b(1:m) is the input right-hand side. x(1:n) is
c      the output solution vector. No input quantities are destroyed, so the routine may be called
c      sequentially with dierent b's.
      INTEGER i,j,jj
      REAL s,tmp(NMAX)
      do 12 j=1,n !Calculate UTB.
          s=0.
          if(w(j).ne.0.)then !Nonzero result only if wj is nonzero.
              do 11 i=1,m
                  s=s+u(i,j)*b(i)
11            enddo 
              s=s/w(j) !This is the divide by wj .
          endif
          tmp(j)=s
12    enddo 
      do 14 j=1,n !Matrix multiply by V to get answer.
          s=0.
          do 13 jj=1,n
              s=s+v(j,jj)*tmp(jj)
13        enddo
          x(j)=s
14    enddo 
c      call IEEE('svbksb',iOutExcept)
      return
      END SUBROUTINE
      
      FUNCTION pythag(a,b)
      REAL a,b,pythag
c     Computes (a2 + b2)1=2 without destructive underflow or overflow.      
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
        pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
c      call IEEE('pythag',iOutExcept)
      return
      END FUNCTION
      
      subroutine solve_proc(NOEL,tolerance,ibc,bc,DROT,DFGRD0,TSTART
     #    ,TEND,time_pro,dtime,STATEV1,STRESS1,STRAN,time,TEMP)
      
      use dimensions
      
      integer, intent(in) :: ibc(NTENS),NOEL
      real, intent(in) :: tolerance,bc(NTENS),TSTART,TEND,time_pro
      
      real, intent(inout) :: STATEV1(NSTATV),dtime
      
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      character(len=150) :: cmname
      
      real :: bc_t_dt(NTENS),bc_r(6),aux66(NTENS,NTENS)
      real :: aux67(NTENS,NTENS),w(NTENS),v(NTENS,NTENS)
      real :: CORRECTION(NTENS),aux6(NTENS),aux7(NTENS)
      real :: STRAN_START(NTENS),STRESS_START(NTENS),STRESS1(NTENS)
      real :: DDFGRD(3,3),aux33(3,3),aux34(3,3),aux3(3),DUMMY(20)
      
      integer NPT
      
      NPT = NOEL
      
      !define initial guess for dstran
      do i=1,NTENS
          DSTRAN(i)=0.0
      enddo
      TEMP=TSTART
      write(222+NOEL,'(10000e32.16)') STRAN,STRESS1,TIME(2)
     #        ,STATEV1(NSTATV-200:NSTATV)
c     #        ,STATEV(NSTATV-28-40+1:NSTATV-28)
c     #        ,STATEV(180:203),STATEV(205:228)
      !save state at beginning of process 
      STRAN_START=STRAN
      STRESS_START=STRESS1
      !find size of problem
      isize=NTENS-sum(ibc)                  
      time(1)=0.0
      NSTEPS=116
      OPEN(10,FILE='OneElComp_defgrad.txt')      
      
      DO II=1,NSTEPS    
          !write
          write(*,*) 'Total time: ', time(2), 'NOEL: ', NOEL
          !find bc at t+dt    
          DTEMP=0

        !restore state variables
        STATEV=STATEV1
        STRESS=STRESS1
        !define dfgrd
        READ(10,*)(DUMMY(I),I=1,20)
        DTIME=DUMMY(1)
        TIME(1)=DUMMY(2)
        TIME(2)=DUMMY(2)
        DFGRD0(1,1)= DUMMY(3)
        DFGRD0(1,2)= DUMMY(4)
        DFGRD0(1,3)= DUMMY(5)
        DFGRD0(2,1)= DUMMY(6)
        DFGRD0(2,2)= DUMMY(7)
        DFGRD0(2,3)= DUMMY(8)
        DFGRD0(3,1)= DUMMY(9)
        DFGRD0(3,2)= DUMMY(10)
        DFGRD0(3,3)= DUMMY(11)
        DFGRD1(1,1)= DUMMY(12)
        DFGRD1(1,2)= DUMMY(13)
        DFGRD1(1,3)= DUMMY(14)
        DFGRD1(2,1)= DUMMY(15)
        DFGRD1(2,2)= DUMMY(16)
        DFGRD1(2,3)= DUMMY(17)
        DFGRD1(3,1)= DUMMY(18)
        DFGRD1(3,2)= DUMMY(19)
        DFGRD1(3,3)= DUMMY(20)              
              
              
        PNEWDT = 1.0
        NDI=3
        !find STRESS and DDSDDE
            CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1            DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2            PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,
     3            NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
     4            NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

          !store state variables
          STATEV1=STATEV
          !store stress
          STRESS1=STRESS
          !update variables
C          time(1)=time(1)+dtime
C          TIME(2)=TIME(2)+dtime
          !find stran
          STRAN=STRAN+DSTRAN
          !find temperature
          TEMP=TEMP+DTEMP
          !update def grad
C          DFGRD0=DFGRD1
          !write output
          write(222+NOEL,'(10000e32.16)') TIME(2),STRESS
c     #            ,STATEV(NSTATV-28-40+1:NSTATV-28)
c     #            ,STATEV(180:203),STATEV(205:228)
          enddo    
          CLOSE(10)
      
      end subroutine solve_proc
      
      subroutine statev_inc
      
      use dimensions
      
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      character(len=150) :: cmname
      
      real :: bc_t_dt(NTENS),bc_r(6),aux66(NTENS,NTENS)
      real :: aux67(NTENS,NTENS),w(NTENS),v(NTENS,NTENS)
      real :: CORRECTION(NTENS),aux6(NTENS),aux7(NTENS)
      real :: STRAN_START(NTENS),STRESS_START(NTENS),STRESS1(NTENS)
      real :: DDFGRD(3,3),aux33(3,3),aux34(3,3),aux3(3)
      
      do i=1,3
          DROT(i,i)=1.0
          DFGRD0(i,i)=1.0
      enddo
      
      !call umat (1st call to sub)
      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1    DRPLDE,DRPLDT,STRAN,DSTRAN*0.0,TIME,DTIME,TEMP,DTEMP,
     2    PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,
     3    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
     4    NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      !read initial state
      open(333,file='statev.in')
      read(333,*)
      read(333,*)
      read(333,*) (DSTRAN(i),i=1,NTENS),((DROT(i,j),i=1,3),j=1,3)
     #    ,(STRAN(i),i=1,NTENS),(STRESS(i),i=1,NTENS)
      read(333,*) 
      do ns=1,NSTATV
          read(333,*) STATEV(ns)
      enddo
      
      do i=1,3
          DROT(i,i)=1.0
          DFGRD0(i,i)=1.0
      enddo
      
      time=1.
      
      !call umat
      CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,
     1    DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2    PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,
     3    NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
     4    NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      
      
      end subroutine statev_inc
      
      end module miscellaneous_sub_main