!******************************************************
! Copyright (c) 2016 Shunsuke A. Sato                 *
! Released under the MIT license                      *
! http://opensource.org/licenses/mit-license.php      *
!******************************************************
module gloval_variables

! parameter
  integer,parameter :: dp = kind(0d0), zp = kind((1d0,1d0))
  complex(zp),parameter :: zI = (0d0,1d0)
  real(dp),parameter :: pi=3.14159265359d0

! mesh
  integer :: Nx
  real(dp) :: length_x,dx
  real(dp),allocatable :: xn(:)

! GS: wave-function, density, potential and so on
  integer :: N_orbit, Nscf, Ncg
  real(dp),allocatable :: wfn(:,:), rho(:), v_ext(:), esp(:), esp_res(:)
  
  
end module gloval_variables

program main
  use gloval_variables
  implicit none

  write(*,'(A)')'Start qm1d'

  call input
  call mesh

  call preparation_GS

  call ground_state  
  write(*,'(A)')'Complete qm1d'

end program main
!=======10========20========30========40========50========60========70========80========90=======100
subroutine input
 use gloval_variables
 implicit none

! == input parameter == !                                                                                                                       
  write(*,'(A)')'===== Input parameter ============================================================='
  write(*,'(A)')

  N_orbit = 4
  Nx = 1200
  length_x = 10d0*pi
  Nscf = 100
  Ncg = 200

  write(*,'(A,2x,I4)')'N_orbit =',N_orbit
  write(*,'(A,2x,I4)')'Nx =',Nx
  write(*,'(A,2x,e26.16e3)')'length_x =',length_x

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete Input parameter ==================================================='

  return
end subroutine input
!=======10========20========30========40========50========60========70========80========90=======100
subroutine mesh
  use gloval_variables
  implicit none
  integer :: ix 

  write(*,'(A)')'===== Making mesh ================================================================'
  write(*,'(A)')
  write(*,'(A)')

  allocate(xn(0:Nx))
  dx = length_x/dble(Nx)
  
  do ix = 0,Nx
     xn(ix) = dx*dble(ix) - 0.5d0*length_x
  end do

  write(*,'(A)')'===== Complete Making mesh ========================================================'
  return
end subroutine mesh
!=======10========20========30========40========50========60========70========80========90=======100
subroutine preparation_GS
  use gloval_variables
  implicit none
  integer :: ix,iorb
  real(dp) :: tmp

  write(*,'(A)')'===== preparatin_GS =============================================================='
  write(*,'(A)')
  write(*,'(A)')

  allocate(wfn(0:Nx,N_orbit), rho(0:Nx), v_ext(0:Nx), esp(N_orbit), esp_res(N_orbit))

  write(*,'(A)')'=== preparing initial wave-function ===='
! wfn(0,:) == wfn(Nx,:) == 0 
  wfn(:,:) = 0d0
  do iorb=1,N_orbit
     do ix=1,Nx-1
        call random_number(tmp)
        tmp = tmp-0.5d0
        wfn(ix,iorb) = tmp
     end do

     tmp = sum(wfn(:,iorb)**2)*dx
     wfn(:,iorb)=wfn(:,iorb)/sqrt(tmp)

  end do

  write(*,'(A)')'=== preparing external potential ===='
 v_ext(:)=0d0
!  do ix=1,Nx-1
!     v_ext(ix) = 0.5d0*xn(ix)**2
!  end do

  write(*,'(A)')
  write(*,'(A)')
  write(*,'(A)')'===== Complete preparatin_GS ====================================================='

  return
end subroutine preparation_GS
!=======10========20========30========40========50========60========70========80========90=======100
subroutine hpsi(f,hf)
  use gloval_variables
  implicit none
  real(dp) :: f(0:Nx), hf(0:Nx)
  integer :: ix
  real(dp) :: c0,c1

! three-points formula  
  c0=-0.5d0*(-2d0/dx**2)
  c1=-0.5d0*(1d0/dx**2)

  do ix=1,Nx-1
     hf(ix)=c1*f(ix+1)+c0*f(ix)+c1*f(ix-1)
  end do

  hf(:)=hf(:)+v_ext(:)*f(:)
  hf(0)=0d0; hf(Nx)=0d0 ! Boundary condition  

  return
end subroutine hpsi
!=======10========20========30========40========50========60========70========80========90=======100
subroutine zhpsi(zf,zhf)
  use gloval_variables
  implicit none
  complex(zp) :: zf(0:Nx), zhf(0:Nx)
  integer :: ix
  complex(zp) :: c0,c1

! three-points formula  
  c0=-0.5d0*(-2d0/dx**2)
  c1=-0.5d0*(1d0/dx**2)

  do ix=1,Nx-1
     zhf(ix)=c1*zf(ix+1)+c0*zf(ix)+c1*zf(ix-1)
  end do

  zhf(:)=zhf(:)+v_ext(:)*zf(:)
  zhf(0)=0d0; zhf(Nx)=0d0 ! Boundary condition  

  return
end subroutine zhpsi
!=======10========20========30========40========50========60========70========80========90=======100
subroutine ground_state
  use gloval_variables
  implicit none
  integer :: iter_scf,iorb

  write(*,'(A)')'===== Ground state calculation ===================================================='
  write(*,'(A)')
  write(*,'(A)')

  do iter_scf=1,Nscf
     call conjugate_gradient_method
     write(*,'(A,2x,I5)')'iter_scf =',iter_scf
     write(*,'(A)')'iorb,     esp,     esp_res'
     do iorb=1,N_orbit
        write(*,'(I4,3x,e16.6e3,3x,e16.6e3)')iorb,esp(iorb),esp_res(iorb)
     end do
  end do

  return
end subroutine ground_state
!=======10========20========30========40========50========60========70========80========90=======100
subroutine conjugate_gradient_method
  use gloval_variables
  implicit none
  real(dp) :: xvec(0:Nx),pvec(0:Nx),rvec(0:Nx),hxvec(0:Nx),gvec(0:Nx),hpvec(0:Nx)
  real(dp) :: xx,pp,xp,xhx,php,xhp,esp_t,esp_res_t,gg,gg0
  real(dp) :: ss,lambda,alpha,beta,aa,bb,cc
  integer :: iorb,iorb_t,ix,iter_cg

  do iorb=1,N_orbit
     xvec(:)=wfn(:,iorb)
     do iorb_t=1,iorb-1
        ss=sum(wfn(:,iorb_t)*xvec(:))*dx
        xvec(:)=xvec(:)-ss*wfn(:,iorb_t)
     end do
        
     call hpsi(xvec,hxvec)
     do iorb_t=1,iorb-1
        ss=sum(wfn(:,iorb_t)*hxvec(:))*dx
        hxvec(:)=hxvec(:)-ss*wfn(:,iorb_t)
     end do

     xx=sum(xvec(:)**2)*dx
     xhx=sum(xvec(:)*hxvec(:))*dx
     lambda=xhx/xx

     do iter_cg=1,Ncg
        gvec(:)=2d0*(hxvec(:)-lambda*xvec(:))/xx
        do iorb_t=1,iorb-1
           ss=sum(wfn(:,iorb_t)*gvec(:))*dx
           gvec(:)=gvec(:)-ss*wfn(:,iorb_t)
        end do
        
        gg0=sum(gvec(:)**2)*dx
        select case(iter_cg)
        case(1)
           pvec(:)=-gvec(:)
        case default
           beta=gg0/gg
           pvec(:)=-gvec(:)+beta*pvec(:)
           do iorb_t=1,iorb-1
              ss=sum(wfn(:,iorb_t)*pvec(:))*dx
              pvec(:)=pvec(:)-ss*wfn(:,iorb_t)
           end do
        end select
        gg=gg0

        call hpsi(pvec,hpvec)
        do iorb_t=1,iorb-1
           ss=sum(wfn(:,iorb_t)*hpvec(:))*dx
           hpvec(:)=hpvec(:)-ss*wfn(:,iorb_t)
        end do

        pp=sum(pvec(:)**2)*dx
        php=sum(pvec(:)*hpvec(:))*dx
        xp=sum(xvec(:)*pvec(:))*dx
        xhp=sum(hxvec(:)*pvec(:))*dx

        aa=php*xp-xhp*pp
        bb=php*xx-xhx*pp
        cc=xhp*xx-xhx*xp
        ss=bb**2-4d0*aa*cc
        if(ss > 0d0)then
           alpha=(-bb+sqrt(ss))/(2d0*aa)
        else
           exit
        end if
        
        xvec(:)=xvec(:)+alpha*pvec(:)

        call hpsi(xvec,hxvec)
        do iorb_t=1,iorb-1
           ss=sum(wfn(:,iorb_t)*hxvec(:))*dx
           hxvec(:)=hxvec(:)-ss*wfn(:,iorb_t)
        end do

        xx=sum(xvec(:)**2)*dx
        xhx=sum(xvec(:)*hxvec(:))*dx
        lambda=xhx/xx

     end do
     xvec(:)=xvec(:)/sqrt(xx)
     call hpsi(xvec,hxvec)
     esp(iorb)=sum(xvec(:)*hxvec(:))*dx
     esp_res(iorb)=sum((hxvec(:)-esp(iorb)*xvec(:))**2)*dx
     wfn(:,iorb)=xvec

  end do

  return
end subroutine conjugate_gradient_method
!=======10========20========30========40========50========60========70========80========90=======100
!=======10========20========30========40========50========60========70========80========90=======100
