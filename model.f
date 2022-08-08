!------------------------------------------------------------------------
        program main
!------------------------------------------------------------------------
! 08/02/22 toy model to recover EOM from dynamics 
        implicit real*8 (a-h,o-z)
        real*8 x0(2), p0(2)
        parameter(fs_to_au=41.341375d0)
        parameter(bohr_to_A=0.529177d0)
        common/pot_par/s(6,500) !cations: x,y,sigma,eta
! position, cation size r0=sigma*1.122,hardness parameter

        open(5,file = 'IN')
        read(5,*) dt, kmax, kout !time step, N large steps, and substeps
        read(5,*) x0, p0, am !initial position, mom, particle mass
        read(5,*) ncat,dcat,ycel !number of cations&spacing,pore size
        read(5,*) sig,eta,cf  ! particle size, hardness, coulomb coeff

        dt = dt*fs_to_au ! in a.u.
!        dcat = dcat/bohr_to_A ! in bohr
!        ycel = ycel/bohr_to_A ! in bohr
        write(*,*) 'spacing and ycel in bohr=',dcat,ycel 

        call define_pot(ncat,dcat,ycel,sig,eta,cf)

        call prop(dt,kmax,kout,x0,p0,am,ncat)

        stop
        end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
        subroutine prop(dt,kmax,kout,x,p,am,ncat)
! Trotter propagator for a classical trajectory 
        implicit real*8 (a-h,o-z)
        real*8 x(2), p(2),dv(2)
        parameter(fs_to_au=41.341375d0,ps_to_au=41341.375d0)
        parameter(bohr_to_A=0.529177d0)

        open(11,file='traj.dat')
        write(11,*)'# t [ps],x(1),x(2),p(1),p(2),E [a.u.]'

        dt = dt/kout
        dt2 = dt/2d0
        tau = dt/am

        t = 0d0

        call pot(ncat,x,V,dv)

        write(11,*)t/ps_to_au,x,p,energy(p,V,am)

        do k = 1,kmax
         do kk = 1,kout
         
         do i = 1,2
          p(i) = p(i) - dv(i)*dt2
         enddo

         do i = 1,2
          x(i) = x(i) + p(i)*tau
         enddo

         call pot(ncat,x,V,dv)

         do i = 1,2
          p(i) = p(i) - dv(i)*dt2
         enddo

         enddo

         t = t+dt*kout

         write(11,*)t/ps_to_au,x,p,energy(p,V,am)

        enddo

        close(11)

        return
        end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
        function energy(p,V,am)
! computes classical energy of a trajectory
        implicit real*8 (a-h,o-z)
        real*8 p(2)

        energy = V+(p(1)**2+p(2)**2)/2d0/am

        return
        end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
        subroutine define_pot(ncat,dcat,ycel,sig,eta,cf)
! put cations on a grid
        implicit real*8 (a-h,o-z)
        common/pot_par/s(6,500) !cations: x,y,sigma,eta
        parameter(bohr_to_A=0.529177d0)
! position, cation size r0=sigma*1.122,hardness parameter

        xmn = -dcat*(ncat/2-1)/2d0
!        xmn = -dcat*(ncat/4-1)/2d0

        ymn = -ycel/2d0

        ymx =  ycel/2d0

        ic = 0

!        do n = 1,ncat,4 ! cation positions
        do n = 1,ncat,2 ! cation positions

          ic = ic+1

          s(1,ic) = xmn+dcat*(n-1)/2d0

          s(2,ic) = ymn

          ic = ic+1

          s(1,ic) = s(1,n)+dcat/2d0

         s(2,ic) = ymx

!          ic = ic+1

!          s(1,ic) = s(1,n)+dcat

!          s(2,ic) = ymn-ycel

!          ic = ic+1

!          s(1,ic) = s(1,n)+dcat*1.5d0

!          s(2,ic) = ymx+ycel

         enddo                

         ncat = ic

         write(*,*) 'number of cations =',ic 
         write(*,*) 'spacing and cell size =',dcat,ycel
!         write(*,*) 'size and harndess =',sig,eta
         write(*,*) 'r0, eta, D for Morse =',sig,eta,cf

         do n = 1,ncat ! cation size and hardness

          s(3,n) = sig ! r0  Morse V=D*(exp(-eta*(r-r0))-1)^2+D
          s(4,n) = eta ! eta 
          s(5,n) = cf  ! D

         enddo

! write cation params 
         open(10,file = 'cation.dat')
         write(10,*) '#cations:  x[A], y[A], r0[a0], eta[1/a0], D[Eh]'
         do n = 1,ncat 
!          write(10,*) s(1,n)*bohr_to_A,s(2,n)*bohr_to_A,(s(k,n),k=3,5)
          write(10,*) (s(k,n),k=1,5) !in bohr
         enddo
         close(10)

         return
         end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
        subroutine pot(ncat,x,v,dv)
! computes LJ + Coulomb potential                             
        implicit real*8 (a-h,o-z)
        real*8 x(2),dv(2),dx(2)
        parameter(Cf = 1d0) ! Coulomb strength
        common/pot_par/s(6,500)

!        ncat = 1 ! debug
!        write(*,*) 'anion',x
!        write(*,*) 'cat',(s(k,1),k=1,5)

        V = 0d0

        dv(1) = 0d0

        dv(2) = 0d0

        do n = 1,ncat ! loop over cations

         r = dist(x,s(1,n),dx) !|r_cation-r_anion|

         y = exp(-s(4,n)*(r-s(3,n)))

         V = V + s(5,n)*(y*y-2d0*y)

         do k=1,2

          dv(k) = dv(k)-2d0*s(4,n)*s(5,n)*(y-1d0)*y*dx(k)

         enddo

!LJ         rs = s(3,n)/r ! sigma/r 
!LJ         V = V+ s(4,n)*(rs**12-rs**6)+s(5,n)/r !LJ+Coulomb
!LJ         do k = 1,2 ! loop over dimensions
!LJ          dv(k) = dv(k)+s(4,n)*(-12d0*rs**12+6d0*rs**6)*dx(k)/r
!LJ     &                 -s(5,n)/r**2*dx(k)
!LJ
!LJ          enddo

         enddo

!         write(*,*) r,dx
!         write(*,*) v,dv
!         stop

         return
         end
!------------------------------------------------------------------------


!------------------------------------------------------------------------
        function dist(r1,r2,dx)
! r=|r1-r2| and grad_r1(r)
        implicit real*8 (a-h,o-z)
        real*8 r1(2),r2(2),dx(2)

        d1 = r1(1)-r2(1)

        d2 = r1(2)-r2(2)

        dist = dsqrt(d1*d1+d2*d2)

        dx(1) = d1/dist

        dx(2) = d2/dist

        return
        end







