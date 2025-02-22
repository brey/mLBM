!======================================
! Lattice Kinetic Code
! 3D Module
! MPI Port
!======================================
      module lkcmod

      use mpi

      implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer status(MPI_STATUS_SIZE),ierr
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      type bn
        integer::xn,yn,kn,ln
        real(8)::dn
        type(bn),pointer::next
      end type

      type(bn),pointer::ibfirst,iblast

      real(8) cm(0:6,3),wm(0:6),b_m(3),w(0:18)

      real(8),dimension(:,:,:,:),allocatable:: fa,df,u,b,fm
      real(8),dimension(:,:,:,:,:),allocatable:: ba,db
      real(8),dimension(:,:,:),allocatable:: dens,ctot,vtot
      real(8),dimension(:,:),allocatable::ms,invms,c
      integer,allocatable::e(:,:),lb(:)
      real(8),allocatable::c_(:),s(:),phi(:),meq(:)

      integer ks,kend,ksize,msize,mtype
      integer ni,nj,nk,nt,kl,kf,n1,l,l1
      integer i,j,k,inode,is,ir
      integer m,mbs,mbe,nframes,opt
      integer numtasks,taskid,ci
      integer::it=0,ncount=0,ik
      integer:: nbit,nn

      integer em(0:6,3)

      real(8) Ha,kel,h,visc,time,length
      real(8) bx,by,bz,x1,x2,x3,ts,pi,a
      real(8) tu,tm,relax,u0,r0,crit
      real(8) dens0,ddenmax,dden,cs02,ddenmaxt
      real(8) u_tot,b_tot,u_c,b_c
      real(8) dmax,den
      real(8) dd,div,divmax
      real(8) u2
      real(8) wxx,we,wej

      character(2) problem,magn
      character(4) bgk_option
      character(3) mode

!------------------------
      data (em(0,j),j=1,3)/0,0,0/
      data (em(1,j),j=1,3)/1,0,0/
      data (em(2,j),j=1,3)/0,0,-1/
      data (em(3,j),j=1,3)/-1,0,0/
      data (em(4,j),j=1,3)/0,0,1/
      data (em(5,j),j=1,3)/0,-1,0/
      data (em(6,j),j=1,3)/0,1,0/

!------------------------
!     data (phi(j),j=1,12)/0.,0.,1.,0.,1.5,0.5,0.25,1.25,1.75,0.5,0.5,0./
!------------------------

       contains

!___________________________
      subroutine cur
!___________________________

     real(8) curx,cury,curz,vortx,vorty,vortz,et,et1,et2
     real(8) emagnetic,ekinetic,st1,st2,cor,cors,cpp(numtasks),vpp(numtasks)
     real(8) emt,ekt,et1t,et2t,divmaxt,cmax,vmax,cmaxt,vmaxt
     integer cmpos(3),vmpos(3),cp(3,numtasks),vp(3,numtasks)

     integer lks,lkend

       nn=(ni-1)*(nj-1)*(nk-1)

       divmax=0.d0
        st1=0.d0
        st2=0.d0
       ekinetic=0.d0
       emagnetic=0.d0
       cor=0.d0

         do k=ks,kend
        do i=1,ni
        do j=1,nj

!--------------------------
!  Compute current and divB
!--------------------------

       if(divb(i,j,k).gt.divmax) divmax=divb(i,j,k)

!--------------------------
!  Compute vorticity
!  &
!  Compute total enetgy dissipation
!--------------------------
          vortx=dfi(3,2)-dfi(2,3)
          vorty=dfi(1,3)-dfi(3,1)
          vortz=dfi(2,1)-dfi(1,2)

!----------------------------------
! Alternative current computation
!----------------------------------

!         curx=dbi(3,2)-dbi(2,3)
!         cury=dbi(1,3)-dbi(3,1)
!         curz=dbi(2,1)-dbi(1,2)

          curx=curf(2,3)
          cury=curf(3,1)
          curz=curf(1,2)

!----------------------------------
        vtot(i,j,k)=dsqrt(vortx**2+vorty**2+vortz**2)
          call vorf
        vtot(i,j,k)=vtot(i,j,k)
        ctot(i,j,k)=dsqrt(curx**2+cury**2+curz**2)

!       if(ncount.eq.nframes) then

!       write(10,101) vtot(i,j,k)
!       write(9,101) ctot(i,j,k)

!       endif


!---------------------
! Compute Lorenz force
!---------------------

      fm(i,j,k,1)=b(3,i,j,k)*cury-b(2,i,j,k)*curz
      fm(i,j,k,2)=-b(3,i,j,k)*curx+b(1,i,j,k)*curz
      fm(i,j,k,3)=b(2,i,j,k)*curx-b(1,i,j,k)*cury

      if(magn=='eq'.and.mode=='bgk') fm(i,j,k,:)=0.d0

         enddo
        enddo
        enddo

!-------------------------
! Integrable quantities
!-------------------------

         lks=ks
         lkend=kend

       if(kend==nk) lkend=lkend-1

         do k=lks,lkend
        do i=1,ni-1
        do j=1,nj-1

!---------------------
! Compute enstrophies
!---------------------
!       st1=st1+vtot(i,j,k)**2/2.d0
        st1=st1+vtot(i,j,k)/2.d0
        st2=st2+ctot(i,j,k)**2/2.d0
!---------------------
! Compute energies
!---------------------
      u_tot=dsqrt(u(1,i,j,k)**2+u(2,i,j,k)**2+u(3,i,j,k)**2)
      b_tot=dsqrt(b(1,i,j,k)**2+b(2,i,j,k)**2+b(3,i,j,k)**2)

      cor=cor+dot_product(u(:,i,j,k),b(:,i,j,k))
      ekinetic=ekinetic+u_tot**2/2.d0
      emagnetic=emagnetic+b_tot**2/2.d0

         enddo
        enddo
        enddo

      et1=st1*length/real(nj-1)
!     et1=st1
      et2=st2*length/real(nj-1)


      ekinetic=ekinetic/real(nn)
      emagnetic=emagnetic/real(nn)
      cor=cor/real(nn)

      cmax=maxval(ctot)
      vmax=maxval(vtot)

      cmpos=maxloc(ctot)
      vmpos=maxloc(vtot)

       call MPI_Reduce(divmax,divmaxt,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

       call MPI_Reduce(emagnetic,emt,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       call MPI_Reduce(ekinetic,ekt,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       call MPI_Reduce(et1,et1t,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       call MPI_Reduce(et2,et2t,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

! Maximum current
       call MPI_Reduce(cmax,cmaxt,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

! Maximum Vorticity
       call MPI_Reduce(vmax,vmaxt,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

       call MPI_Reduce(cor,cors,1,&
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

! Position of maximum current
       call MPI_Gather(cmpos,3,&
        MPI_INTEGER,cp,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       call MPI_Gather(cmax,1,&
        MPI_DOUBLE_PRECISION,cpp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

! Position of maximum Vorticity
       call MPI_Gather(vmpos,3,&
        MPI_INTEGER,vp,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       call MPI_Gather(vmax,1,&
        MPI_DOUBLE_PRECISION,vpp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      if(taskid.eq.0) then

      write(30,104) real(it)*ts,et1t/u0**2/length**3,et2t/u0**2/length**3,(et1t+et2t)/u0**2/length**3
      write(31,104) real(it)*ts,ekt/u0**2,emt/u0**2
!     write(31,104) real(it),dsqrt(2.*ekt), emt/u0**2, (ekt+emt)/u0**2
      write(11,102) real(it)*ts,divmax
      write(12,103) real(it)*ts,cmaxt,vmaxt
      write(15,102) real(it)*ts,cors/(ekt+emt)
      write(13,111) cp(:,maxloc(cpp))
      write(14,111) vp(:,maxloc(vpp))

!      print*,dmax
!      print*,divmaxt

      endif

 101  format(e12.6)
 111  format(3(1x,i3))
 102  format(2(2x,e12.6))
 103  format(3(2x,e12.6))
 104  format(4(1x,e12.6))

!--------------------------
!   END SUBROUTINE CUR
!--------------------------
      end subroutine

!---------------------------
!   DISCRIBE INNER BOUNDARY
!---------------------------
      subroutine innerb

      real(8) y1,y2,y3,r,rr,dr

         r=dsqrt(x1**2+x2**2)

         do l=1,nbit-1
         y1=x1+e(l,1)*dd
         y2=x2+e(l,2)*dd
         y3=x3+e(l,3)*dd
         rr=dsqrt(y1**2+y2**2)

         if(r<r0.and.rr>r0) then

         dr=-x1*dcos(phi(l)*pi)-x2*dsin(phi(l)*pi)-&
     dsqrt(r0**2-x1**2-x2**2+(x1*dcos(phi(l)*pi)+x2*dsin(phi(l)*pi))**2)
         if(dr<0.) then
         dr=-x1*dcos(phi(l)*pi)-x2*dsin(phi(l)*pi)+&
     dsqrt(r0**2-x1**2-x2**2+(x1*dcos(phi(l)*pi)+x2*dsin(phi(l)*pi))**2)
         endif

         if(associated(ibfirst)) then

         allocate(iblast%next)

         iblast=>iblast%next

         iblast%next=>null()
         iblast%xn=i
         iblast%yn=j
         iblast%kn=k
         iblast%ln=l
         iblast%dn=dr/dsqrt((x1-y1)**2+(x2-y2)**2)

         else

         allocate(ibfirst)

         iblast=>ibfirst

         iblast%next=>null()
         iblast%xn=i
         iblast%yn=j
         iblast%kn=k
         iblast%ln=l
         iblast%dn=dr/dsqrt((x1-y1)**2+(x2-y2)**2)

         endif


         endif

         enddo

      end subroutine
!------------------------
!   BOUNDARY CONDITIONS
!------------------------

!------------------------
!  D.Yu et al (2003)
!------------------------
      subroutine bound

      real(8) fwil,ff,fff,y1,y2,y3,ddd
      real(8), pointer::d
      integer, pointer::ix,iy,iz,il
!------------------------
!   INLET TOROIDAL BOUNDARY
!------------------------

      do while(associated(ibfirst))

         ix=>ibfirst%xn
         iy=>ibfirst%yn
         iz=>ibfirst%kn
         il=>ibfirst%ln
         d=>ibfirst%dn

            x1=-0.5d0+(ix-1)*dd
            x2=-0.5d0+(iy-1)*dd
            x3=0.d0+(iz-1)*dd
         y1=x1+e(il,1)*dd
         y2=x2+e(il,2)*dd
         y3=x3+e(il,3)*dd

       ddd=dsqrt((x1-y1)**2+(x2-y2)**2)

       write(20,*) x1+d*ddd*dcos(phi(il)*pi),x2+d*ddd*dsin(phi(il)*pi)


         ff=fa(lb(il),ix-e(il,1),iy-e(il,2),iz-e(il,3))
         fff=fa(lb(il),ix-2*e(il,1),iy-2*e(il,2),iz-2*e(il,3))

!----------------------------------------
!  fa
!  No-Slip Boundary conditions
!----------------------------------------
         fwil=fa(il,ix,iy,iz)+d*(fa(il,ix+e(il,1),iy+e(il,2),iz+e(il,3))-fa(il,ix,iy,iz))
         fa(il,ix+e(il,1),iy+e(il,2),iz+e(il,3))=df(lb(il),ix+e(il,1),iy+e(il,2),iz+e(il,3))
!----------------------------------------
!  linear
!----------------------------------------
         fa(lb(il),ix,iy,iz)=fwil+d/(1.d0-d)*(ff-fwil)

!----------------------------------------
!  2nd order
!----------------------------------------
!        fa(lb(il),ix,iy,iz)=fwil*2.d0/((1.d0+d)*(2.d0+d))+ff*2.d0*d/(1.d0+d)-fff*d/(2.d0+d)

!----------------------------------------
!  ba
!  Dirichlet Conditions
!----------------------------------------

       ibfirst=>ibfirst%next

       enddo

        end subroutine

!-----------------------------
!   compute beq with new u,b
!-----------------------------
      real(8) function beq(l,m)
      integer l,m

      select case(m)

      case(1)

      beq=wm(l)*(b(1,i,j,k)&
         +4.d0*cm(l,2)*(u(2,i,j,k)*b(1,i,j,k)-b(2,i,j,k)*u(1,i,j,k))&
         +4.d0*cm(l,3)*(u(3,i,j,k)*b(1,i,j,k)-b(3,i,j,k)*u(1,i,j,k)))

      case(2)

      beq=wm(l)*(b(2,i,j,k)&
         +4.d0*cm(l,1)*(u(1,i,j,k)*b(2,i,j,k)-b(1,i,j,k)*u(2,i,j,k))&
         +4.d0*cm(l,3)*(u(3,i,j,k)*b(2,i,j,k)-b(3,i,j,k)*u(2,i,j,k)))

      case(3)

      beq=wm(l)*(b(3,i,j,k)&
         +4.d0*cm(l,1)*(u(1,i,j,k)*b(3,i,j,k)-b(1,i,j,k)*u(3,i,j,k))&
         +4.d0*cm(l,2)*(u(2,i,j,k)*b(3,i,j,k)-b(2,i,j,k)*u(3,i,j,k)))

      end select

       end function

!------------------------------------
!   compute velocity derivatives with stencil
!------------------------------------
      real(8) function dfi(ld,m)
      integer ld,m
      integer ip,jp,kp,im,jm,km

          jp=j+1
          jm=j-1
          kp=k+1
          km=k-1
          ip=i+1
          im=i-1

          if(jp>nj) jp=2
          if(jm<1) jm=nj-1
          if(ip>ni) ip=2
          if(im<1) im=ni-1

       if(numtasks.eq.1) then
          if(kp>nk) kp=2
          if(km<1) km=nk-1
       endif

      select case(nbit)


      case(13)

      select case(m)

      case(1)

          dfi=1.d0/8.d0*(u(ld,ip,jp,k)-u(ld,im,jp,k)&
                   +u(ld,ip,jm,k)-u(ld,im,jm,k)&
                   +u(ld,ip,j,kp)-u(ld,im,j,kp)&
                   +u(ld,ip,j,km)-u(ld,im,j,km))

      case(2)

          dfi=1.d0/8.d0*(u(ld,i,jp,kp)-u(ld,i,jm,kp)&
                   +u(ld,i,jp,km)-u(ld,i,jm,km)&
                   +u(ld,ip,jp,k)-u(ld,ip,jm,k)&
                   +u(ld,im,jp,k)-u(ld,im,jm,k))

      case(3)

          dfi=1.d0/8.d0*(u(ld,ip,j,kp)-u(ld,ip,j,km)&
                   +u(ld,im,j,kp)-u(ld,im,j,km)&
                   +u(ld,i,jp,kp)-u(ld,i,jp,km)&
                   +u(ld,i,jm,kp)-u(ld,i,jm,km))

      end select

      case(19)

      select case(m)

      case(1)

          dfi=u(ld,ip,j,k)-u(ld,im,j,k)-1.d0/8.d0*(u(ld,ip,jp,k)-u(ld,im,jp,k)&
                   +u(ld,ip,jm,k)-u(ld,im,jm,k)&
                   +u(ld,ip,j,kp)-u(ld,im,j,kp)&
                   +u(ld,ip,j,km)-u(ld,im,j,km))

          dfi=(u(ld,ip,j,k)-u(ld,im,j,k))/2.d0

      case(2)

          dfi=u(ld,i,jp,k)-u(ld,i,jm,k)-1.d0/8.d0*(u(ld,i,jp,kp)-u(ld,i,jm,kp)&
                   +u(ld,i,jp,km)-u(ld,i,jm,km)&
                   +u(ld,ip,jp,k)-u(ld,ip,jm,k)&
                   +u(ld,im,jp,k)-u(ld,im,jm,k))

          dfi=(u(ld,i,jp,k)-u(ld,i,jm,k))/2.d0

      case(3)

          dfi=u(ld,i,j,kp)-u(ld,i,j,km)-1.d0/8.d0*(u(ld,ip,j,kp)-u(ld,ip,j,km)&
                   +u(ld,im,j,kp)-u(ld,im,j,km)&
                   +u(ld,i,jp,kp)-u(ld,i,jp,km)&
                   +u(ld,i,jm,kp)-u(ld,i,jm,km))

          dfi=(u(ld,i,j,kp)-u(ld,i,j,km))/2.d0

      end select

      end select

       end function

!------------------------------------
!   compute magnetic derivatives with stencil
!------------------------------------
      real(8) function dbi(ld,m)
      integer ld,m
      integer ip,jp,kp,im,jm,km

          jp=j+1
          jm=j-1
          kp=k+1
          km=k-1
          ip=i+1
          im=i-1

          if(jp>nj) jp=2
          if(jm<1) jm=nj-1
          if(ip>ni) ip=2
          if(im<1) im=ni-1

       if(numtasks.eq.1) then
          if(kp>nk) kp=2
          if(km<1) km=nk-1
       endif

      select case(m)

      case(1)

          dbi=1.d0/8.d0*(b(ld,ip,jp,k)-b(ld,im,jp,k)&
                   +b(ld,ip,jm,k)-b(ld,im,jm,k)&
                   +b(ld,ip,j,kp)-b(ld,im,j,kp)&
                   +b(ld,ip,j,km)-b(ld,im,j,km))

      case(2)

          dbi=1.d0/8.d0*(b(ld,i,jp,kp)-b(ld,i,jm,kp)&
                   +b(ld,i,jp,km)-b(ld,i,jm,km)&
                   +b(ld,ip,jp,k)-b(ld,ip,jm,k)&
                   +b(ld,im,jp,k)-b(ld,im,jm,k))

      case(3)

          dbi=1.d0/8.d0*(b(ld,ip,j,kp)-b(ld,ip,j,km)&
                   +b(ld,im,j,kp)-b(ld,im,j,km)&
                   +b(ld,i,jp,kp)-b(ld,i,jp,km)&
                   +b(ld,i,jm,kp)-b(ld,i,jm,km))

      end select

       end function

!--------------------------
!   OUTPUT DATA
!--------------------------
      subroutine output

      if(taskid.eq.0) then


       do k=ks,kend
        do i=1,ni

!           x1=0.d0+(i-1)*dd
!           x2=0.d0+(j-1)*dd

!        if(dsqrt(x1**2+x2**2)<r0) u(:,i,j,k)=0.


        write(2,103) (u(1,i,j,k),u(2,i,j,k),u(3,i,j,k),j=1,nj)
        write(7,103) (b(1,i,j,k),b(2,i,j,k),b(3,i,j,k),j=1,nj)
        write(3,101) (dens(i,j,k),j=1,nj)
        write(9,101) (ctot(i,j,k),j=1,nj)
        write(10,101) (vtot(i,j,k),j=1,nj)


!       write(2,103) u(1,i,j,k),u(2,i,j,k),u(3,i,j,k)
!       write(7,103) b(1,i,j,k),b(2,i,j,k),b(3,i,j,k)
!       write(3,101) dens(i,j,k)
!       write(9,101) ctot(i,j,k)
!       write(10,101) vtot(i,j,k)

        enddo
       enddo

      endif


      do inode=1,numtasks-1


!-----------
      if(taskid.eq.inode) then

       call MPI_SEND(u(1,0,0,ks),3*n1*n1*ksize,&
        MPI_DOUBLE_PRECISION,0,mtype,MPI_COMM_WORLD,ierr)

      endif


      if(taskid.eq.0) then

       call MPI_RECV(u(1,0,0,ks),3*n1*n1*ksize,&
        MPI_DOUBLE_PRECISION,inode,mtype,MPI_COMM_WORLD,status,ierr)

       endif
!-----------
      if(taskid.eq.inode) then

       call MPI_SEND(b(1,0,0,ks),3*n1*n1*ksize,&
        MPI_DOUBLE_PRECISION,0,mtype,MPI_COMM_WORLD,ierr)

      endif


      if(taskid.eq.0) then

       call MPI_RECV(b(1,0,0,ks),3*n1*n1*ksize,&
        MPI_DOUBLE_PRECISION,inode,mtype,MPI_COMM_WORLD,status,ierr)

      endif
!-----------
      if(taskid.eq.inode) then

       call MPI_SEND(dens(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,0,mtype,MPI_COMM_WORLD,ierr)

      endif


      if(taskid.eq.0) then

       call MPI_RECV(dens(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,inode,mtype,MPI_COMM_WORLD,status,ierr)

      endif
!-----------
      if(taskid.eq.inode) then

       call MPI_SEND(vtot(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,0,mtype,MPI_COMM_WORLD,ierr)

      endif


      if(taskid.eq.0) then

       call MPI_RECV(vtot(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,inode,mtype,MPI_COMM_WORLD,status,ierr)

      endif
!-----------
      if(taskid.eq.inode) then

       call MPI_SEND(ctot(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,0,mtype,MPI_COMM_WORLD,ierr)

      endif


      if(taskid.eq.0) then

       call MPI_RECV(ctot(1,1,ks),ni*nj*ksize,&
        MPI_DOUBLE_PRECISION,inode,mtype,MPI_COMM_WORLD,status,ierr)

      endif
!-----------
      if(taskid.eq.0) then

       do k=ks,kend
        do i=1,ni

        write(2,103) (u(1,i,j,k),u(2,i,j,k),u(3,i,j,k),j=1,nj)
        write(7,103) (b(1,i,j,k),b(2,i,j,k),b(3,i,j,k),j=1,nj)
        write(3,101) (dens(i,j,k),j=1,nj)
        write(9,101) (ctot(i,j,k),j=1,nj)
        write(10,101) (vtot(i,j,k),j=1,nj)

        enddo
       enddo


       endif


       enddo


!--------------------------
  101 format(e12.6)
  103 format(3(1x,e12.6))
!--------------------------
      end subroutine

!-----------------------------
!   compute meq with new u,b
!-----------------------------
      subroutine mrt

      real(8) cf(0:nbit-1)


        do i=1,ni
         do j=1,nj

          do k=ks,kend


      u2=u(1,i,j,k)**2+u(2,i,j,k)**2+u(3,i,j,k)**2

           if(nbit==19) call meq19
           if(nbit==13) call meq13

               cf=matmul(ms,fa(:,i,j,k))

!  forcing term


               cf(3)=cf(3)+fm(i,j,k,1)/2.d0
               cf(5)=cf(5)+fm(i,j,k,2)/2.d0
               cf(7)=cf(7)+fm(i,j,k,3)/2.d0

!  macroscopic
!              dens(i,j)=cf(0)
!              u(1,i,j,k)=cf(1)
!              u(2,i,j,k)=cf(2)
!              u(3,i,j,k)=cf(3)


               cf=cf-s*(cf-meq)

!  forcing term

               cf(3)=cf(3)+fm(i,j,k,1)/2.d0
               cf(5)=cf(5)+fm(i,j,k,2)/2.d0
               cf(7)=cf(7)+fm(i,j,k,3)/2.d0


               df(:,i,j,k)=matmul(invms,cf)

         enddo

         enddo

         enddo

      end subroutine

!---------------------------------
!  compute collision term with bgk
!---------------------------------

      subroutine collf

       do i=1,ni
       do j=1,nj
         do k=ks,kend

            do l=0,nbit-1

        df(l,i,j,k)=fa(l,i,j,k)-1.d0/(tu+0.5d0)*(fa(l,i,j,k)-feq(l))&
                              +tu/(tu+.5d0)*fr(l)

            enddo

         enddo
       enddo
       enddo

      end subroutine

!-----------------------------
!   compute feq with new u,b
!-----------------------------
       real(8) function feq(l)
       integer l

        u_c=c(l,1)*u(1,i,j,k)+c(l,2)*u(2,i,j,k)+c(l,3)*u(3,i,j,k)
        b_c=c(l,1)*b(1,i,j,k)+c(l,2)*b(2,i,j,k)+c(l,3)*b(3,i,j,k)

        u_tot=dsqrt(u(1,i,j,k)**2+u(2,i,j,k)**2+u(3,i,j,k)**2)
        b_tot=dsqrt(b(1,i,j,k)**2+b(2,i,j,k)**2+b(3,i,j,k)**2)

        select case (bgk_option)

       case('comp')

       feq=dens(i,j,k)*w(l)*(&
                                 1.d0+3.d0*u_c&
                                +9.d0/2.d0*u_c**2&
                                -3.d0/2.d0*u_tot**2&
                                )

       case('inco')

       feq=w(l)*(dens(i,j,k)&
                                +3.d0*u_c&
                                +9.d0/2.d0*u_c**2&
                                -3.d0/2.d0*u_tot**2&
                                )

       end select

       if(magn=='eq') then

       feq=feq&
       +w(l)*9.d0/2.d0*(1.d0/3.d0*b_tot**2*c_(l)**2-b_c**2)

       endif

       end function

!-----------------------------
!   compute forcing term
!-----------------------------
      real(8) function fr(l)
      integer l

        u_c=c(l,1)*u(1,i,j,k)+c(l,2)*u(2,i,j,k)+c(l,3)*u(3,i,j,k)

      fr=3.d0*w(l)*dens(i,j,k)*(&
         fm(i,j,k,1)*((c(l,1)-u(1,i,j,k))+3.d0*u_c*c(l,1))&
        +fm(i,j,k,2)*((c(l,2)-u(2,i,j,k))+3.d0*u_c*c(l,2))&
        +fm(i,j,k,3)*((c(l,3)-u(3,i,j,k))+3.d0*u_c*c(l,3))&
             )

       end function
!---------------------------------
!  initialize local quantities
!---------------------------------

      subroutine init


       do i=1,ni
       do j=1,nj
         do k=ks,kend


        select case (mode)


        case('mrt')

        if(nbit==19)  call meq19
        if(nbit==13)  call meq13

          do l=0,nbit-1
               fa(l,i,j,k)=dot_product(invms(l,:),meq(:))
          enddo

        case ('bgk')


         do l=0,nbit-1

             fa(l,i,j,k)=feq(l)

         enddo

        end select

          do l=0,6
             ba(1,l,i,j,k)=beq(l,1)
             ba(2,l,i,j,k)=beq(l,2)
             ba(3,l,i,j,k)=beq(l,3)
          enddo

        enddo
      enddo
      enddo

!     call thermal

      call cur

      ddenmaxt=1.d0
      ddenmax=1.d0
      ik=1

      do while (ddenmaxt>crit)
!--------------------------------
! Initialize pressure
!--------------------------------
      ddenmax=0.d0

       do i=1,ni
       do j=1,nj
         do k=ks,kend

        select case (mode)

        case('mrt')
!  mrt
        den=dot_product(ms(0,:),fa(:,i,j,k))

        case('bgk')

        den=sum(fa(:,i,j,k))

        end select

         dden=dabs(dens(i,j,k)-den)

         if(dabs(dden/dens(i,j,k)).gt.ddenmax) ddenmax=dabs(dden/dens(i,j,k))

         dens(i,j,k)=den


        enddo
         enddo
         enddo


! Maximum ddenmax
       call MPI_Reduce(ddenmax,ddenmaxt,1,&
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

       call MPI_BCAST(ddenmaxt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        if(taskid.eq.0)  print*, ddenmaxt


       if (ik==1) then
        ddenmaxt=1.d0
        ik=0
       endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------
! Compute collision term
!----------------
          if(mode=='mrt') call mrt
          if(mode=='bgk') call collf

          call collb

!----------------
! Triple Periodic Condition
!----------------
!      if(opt.eq.0) then
           do k=ks,kend

           do j=1,nj

                df(:,0,j,k)=df(:,ni-1,j,k)
                df(:,ni+1,j,k)=df(:,2,j,k)

                db(:,:,0,j,k)=db(:,:,ni-1,j,k)
                db(:,:,ni+1,j,k)=db(:,:,2,j,k)
           enddo

           do i=0,ni+1

                df(:,i,0,k)=df(:,i,nj-1,k)
                df(:,i,nj+1,k)=df(:,i,2,k)

                db(:,:,i,0,k)=db(:,:,i,nj-1,k)
                db(:,:,i,nj+1,k)=db(:,:,i,2,k)

           enddo

           enddo
!       endif

!-------------------------------------------------
! Each node sends their boundary collision terms
! for the streaming process
!-------------------------------------------------
      if(numtasks.gt.1)  call mpicom
!-------------------------------------------------
!-----------------
! for serial run
!-----------------
      if(numtasks.eq.1.) then

        df(:,:,:,0)=df(:,:,:,nk-1)
        df(:,:,:,nk+1)=df(:,:,:,2)

        db(:,:,:,:,0)=db(:,:,:,:,nk-1)
        db(:,:,:,:,nk+1)=db(:,:,:,:,2)

      endif
!------------------------
!  streaming fa, ba
!------------------------

      call stream

!     call no_slip

      enddo

      end subroutine

!     subroutine thermal
!

!     real(8) dthx,dthy,laplth,divj

!       temp(:,js)=u0
!       temp(:,jend)=-u0

!     do i=2,ni-1
!      do j=js+1,jend-1

!       dthx=temp(i+1,j)-temp(i-1,j)-1./4.*(temp(i+1,j+1)-temp(i-1,j+1)+temp(i+1,j-1)-temp(i-1,j-1))
!       dthy=temp(i,j+1)-temp(i,j-1)-1./4.*(temp(i+1,j+1)-temp(i+1,j-1)+temp(i-1,j+1)-temp(i-1,j-1))
!       laplth=2.*(temp(i+1,j)+temp(i-1,j)+temp(i,j+1)+temp(i,j-1))&
!              -1./2.*(temp(i+1,j+1)+temp(i-1,j+1)+temp(i-1,j-1)+temp(i+1,j-1))&
!              -6.*temp(i,j)

!       divj=u(1,i+1,j)-u(1,i-1,j)-1./4.*(u(1,i+1,j+1)-u(1,i-1,j+1)+u(1,i+1,j-1)-u(1,i-1,j-1))&
!           +u(2,i,j+1)-u(2,i,j-1)-1./4.*(u(2,i+1,j+1)-u(2,i+1,j-1)+u(2,i-1,j+1)-u(2,i-1,j-1))

!       temp(i,j)=temp(i,j)&
!                 -u(1,i,j)*dthx&
!                 -u(2,i,j)*dthy&
!                 +kth*laplth&
!                 +(gama-1.)*cs02*divj

!      enddo
!     enddo

!       temp(1,:)=temp(2,:)
!       temp(ni,:)=temp(ni-1,:)

!     end subroutine


!     subroutine no_slip

!

!     end subroutine

!----------------------------
!  mpi communication
!----------------------------

      subroutine mpicom

!-------------------------------------------------

      if(mod(taskid,2)==0) then
!-----------------------------------
! EACH NODE SENDS THE RIGHT BOUNDARY
!-----------------------------------
      is=taskid-1
      ir=taskid+1

      if(ir.eq.numtasks) ir=0
      if(is.eq.-1) is=numtasks-1

      kl=kend
      if(ir.eq.0) kl=kend-1

!-----------------
! SEND/RECEIVE
!-----------------

!     mtype=11

!      call MPI_SENDRECV(df(0,0,0,kl),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,df(0,0,0,ks-1),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(df(0,0,0,kl),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       call MPI_RECV(df(0,0,0,ks-1),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)
!-----------
! SEND/RECEIVE
!-----------

!     mtype=12

!      call MPI_SENDRECV(db(1,0,0,0,kl),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,db(1,0,0,0,ks-1),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(db(1,0,0,0,kl),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       call MPI_RECV(db(1,0,0,0,ks-1),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

      endif
      if(mod(taskid,2)/=0) then
!-----------------------------------
! EACH NODE SENDS THE RIGHT BOUNDARY
!-----------------------------------
      is=taskid-1
      ir=taskid+1

      if(ir.eq.numtasks) ir=0
      if(is.eq.-1) is=numtasks-1

      kl=kend
      if(ir.eq.0) kl=kend-1

!-----------------
! SEND/RECEIVE
!-----------------

!     mtype=11

!      call MPI_SENDRECV(df(0,0,0,kl),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,df(0,0,0,ks-1),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(df(0,0,0,ks-1),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(df(0,0,0,kl),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

!-----------
! SEND/RECEIVE
!-----------

!     mtype=12

!      call MPI_SENDRECV(db(1,0,0,0,kl),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,db(1,0,0,0,ks-1),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(db(1,0,0,0,ks-1),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(db(1,0,0,0,kl),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

      endif
      if(mod(taskid,2)==0) then
!-----------------------------------
! EACH NODE SENDS THE LEFT COLUMN
!-----------------------------------

      is=taskid+1
      ir=taskid-1

      if(ir.eq.-1) ir=numtasks-1
      if(is.eq.numtasks) is=0

      kf=ks
      if(ir.eq.(numtasks-1)) kf=ks+1
!-----------
! SEND/RECEIVE
!-----------

!     mtype=13

!      call MPI_SENDRECV(df(0,0,0,kf),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,df(0,0,0,kend+1),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(df(0,0,0,kf),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)


       call MPI_RECV(df(0,0,0,kend+1),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)
!-----------
! SEND/RECEIVE
!-----------

!     mtype=14

!      call MPI_SENDRECV(db(1,0,0,0,kf),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,db(1,0,0,0,kend+1),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(db(1,0,0,0,kf),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)


       call MPI_RECV(db(1,0,0,0,kend+1),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       endif
      if(mod(taskid,2)/=0) then
!-----------------------------------
! EACH NODE SENDS THE LEFT COLUMN
!-----------------------------------

      is=taskid+1
      ir=taskid-1

      if(ir.eq.-1) ir=numtasks-1
      if(is.eq.numtasks) is=0

      kf=ks
      if(ir.eq.(numtasks-1)) kf=ks+1
!-----------
! SEND/RECEIVE
!-----------

!     mtype=13

!      call MPI_SENDRECV(df(0,0,0,kf),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,df(0,0,0,kend+1),nbit*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(df(0,0,0,kend+1),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(df(0,0,0,kf),nbit*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

!-----------
! SEND/RECEIVE
!-----------

!     mtype=14

!      call MPI_SENDRECV(db(1,0,0,0,kf),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,db(1,0,0,0,kend+1),3*7*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(db(1,0,0,0,kend+1),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(db(1,0,0,0,kf),3*7*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       endif


      end subroutine


      subroutine stream


        do i=1,ni
         do j=1,nj
           do k=ks,kend

      do l=0,nbit-1
      fa(l,i,j,k)=df(l,i-e(l,1),j-e(l,2),k-e(l,3))
      enddo

      do l=0,6
      ba(:,l,i,j,k)=db(:,l,i-em(l,1),j-em(l,2),k-em(l,3))
      enddo


         enddo
        enddo
       enddo


      end subroutine


      subroutine update

      if(mod(taskid,2)==0) then
!-----------------------------------
! EACH NODE SENDS THE LAST COLUMN
!-----------------------------------

      is=taskid-1
      ir=taskid+1

      if(ir.eq.numtasks) ir=0
      if(is.eq.-1) is=numtasks-1

      kl=kend
      if(ir.eq.0) kl=kend-1
!-----------
! SEND/RECEIVE
!-----------

!     type=15

!      call MPI_SENDRECV(u(1,0,0,kl),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,u(1,0,0,ks-1),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(u(1,0,0,kl),3*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)


       call MPI_RECV(u(1,0,0,ks-1),3*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)


       endif
      if(mod(taskid,2)/=0) then
!-----------------------------------
! EACH NODE SENDS THE LAST COLUMN
!-----------------------------------

      is=taskid-1
      ir=taskid+1

      if(ir.eq.numtasks) ir=0
      if(is.eq.-1) is=numtasks-1

      kl=kend
      if(ir.eq.0) kl=kend-1
!-----------
! SEND/RECEIVE
!-----------

!     type=15

!      call MPI_SENDRECV(u(1,0,0,kl),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,u(1,0,0,ks-1),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(u(1,0,0,ks-1),3*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(u(1,0,0,kl),3*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       endif

      if(mod(taskid,2)==0) then
!-----------------------------------
! EACH NODE SENDS THE LEFT COLUMN
!-----------------------------------
      is=taskid+1
      ir=taskid-1

      if(ir.eq.-1) ir=numtasks-1
      if(is.eq.numtasks) is=0

      kf=ks
      if(ir.eq.(numtasks-1)) kf=ks+1
!-----------
! SEND/RECEIVE
!-----------

!     mtype=16

!      call MPI_SENDRECV(u(1,0,0,kf),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,u(1,0,0,kend+1),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(u(1,0,0,kf),3*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       call MPI_RECV(u(1,0,0,kend+1),3*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       endif

      if(mod(taskid,2)/=0) then
!-----------------------------------
! EACH NODE SENDS THE LEFT COLUMN
!-----------------------------------
      is=taskid+1
      ir=taskid-1

      if(ir.eq.-1) ir=numtasks-1
      if(is.eq.numtasks) is=0

      kf=ks
      if(ir.eq.(numtasks-1)) kf=ks+1
!-----------
! SEND/RECEIVE
!-----------

!     mtype=16

!      call MPI_SENDRECV(u(1,0,0,kf),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,u(1,0,0,kend+1),3*n1*n1,&
!       MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_RECV(u(1,0,0,kend+1),3*n1*n1,&
        MPI_DOUBLE_PRECISION,is,mtype,MPI_COMM_WORLD,status,ierr)

       call MPI_SEND(u(1,0,0,kf),3*n1*n1,&
        MPI_DOUBLE_PRECISION,ir,mtype,MPI_COMM_WORLD,ierr)

       endif

      end subroutine


      subroutine macro

       dmax=0.d0
!--------------------------------
! Update macroscopic variables
!--------------------------------
       select case (mode)

       case('mrt')

       do i=1,ni
        do j=1,nj
         do k=ks,kend

         dens(i,j,k)=dot_product(ms(0,:),fa(:,i,j,k))
         u(1,i,j,k)=dot_product(ms(3,:),fa(:,i,j,k))
         u(2,i,j,k)=dot_product(ms(5,:),fa(:,i,j,k))
         u(3,i,j,k)=dot_product(ms(7,:),fa(:,i,j,k))

       enddo
        enddo
         enddo

       case('bgk')

           if(bgk_option=='comp') then

       do i=1,ni
        do j=1,nj
         do k=ks,kend

          dens(i,j,k)=sum(fa(:,i,j,k))

          u(1,i,j,k)=dot_product(fa(:,i,j,k),c(:,1))/dens(i,j,k)&
           +dens(i,j,k)*fm(i,j,k,1)/2.d0
          u(2,i,j,k)=dot_product(fa(:,i,j,k),c(:,2))/dens(i,j,k)&
           +dens(i,j,k)*fm(i,j,k,2)/2.d0
          u(3,i,j,k)=dot_product(fa(:,i,j,k),c(:,3))/dens(i,j,k)&
           +dens(i,j,k)*fm(i,j,k,3)/2.d0

       enddo
        enddo
         enddo


           elseif(bgk_option=='inco') then

       do i=1,ni
        do j=1,nj
         do k=ks,kend

          dens(i,j,k)=sum(fa(:,i,j,k))

          u(1,i,j,k)=dot_product(fa(:,i,j,k),c(:,1))&
           +dens(i,j,k)*fm(i,j,k,1)/2.d0
          u(2,i,j,k)=dot_product(fa(:,i,j,k),c(:,2))&
           +dens(i,j,k)*fm(i,j,k,2)/2.d0
          u(3,i,j,k)=dot_product(fa(:,i,j,k),c(:,3))&
           +dens(i,j,k)*fm(i,j,k,3)/2.d0

       enddo
        enddo
         enddo

        endif

        end select

       do i=1,ni
        do j=1,nj
         do k=ks,kend

         b_m(:)=0.d0

         do l=0,6

          b_m(:)=b_m(:)+ba(:,l,i,j,k)

         enddo

         b(:,i,j,k)=b_m(:)

       enddo
        enddo
         enddo


        u(:,0,:,ks:kend)=u(:,ni-1,:,ks:kend)
        u(:,:,0,ks:kend)=u(:,:,nj-1,ks:kend)
        u(:,ni+1,:,ks:kend)=u(:,2,:,ks:kend)
        u(:,:,nj+1,ks:kend)=u(:,:,2,ks:kend)

        b(:,0,:,ks:kend)=b(:,ni-1,:,ks:kend)
        b(:,:,0,ks:kend)=b(:,:,nj-1,ks:kend)
        b(:,ni+1,:,ks:kend)=b(:,2,:,ks:kend)
        b(:,:,nj+1,ks:kend)=b(:,:,2,ks:kend)

      end subroutine

!-----------------------------
!   compute collision term
!-----------------------------
      subroutine collb

       do i=1,ni
       do j=1,nj
         do k=ks,kend

       do m=1,3
        do l=0,6
         db(m,l,i,j,k)=ba(m,l,i,j,k)&
                    -1.d0/(tm+0.5d0)*(ba(m,l,i,j,k)-beq(l,m))
        enddo
       enddo

!====================================================

         enddo   ! next k
        enddo   ! next j
        enddo   ! next i

       end subroutine

!-----------------------------

      real(8) function curf(l,m)
      integer l,m,ib
      real(8) clm,cml,bbq(0:6,2)

      do ib=0,6
      bbq(ib,1)=beq(ib,l)
      bbq(ib,2)=beq(ib,m)
      enddo

!   bbq => l==1, m==2

      clm=(dot_product(cm(:,l),ba(m,:,i,j,k))+1.d0/(2.d0*tm)*dot_product(cm(:,l),bbq(:,2)))/(1.d0+1.d0/(2.d0*tm))
      cml=(dot_product(cm(:,m),ba(l,:,i,j,k))+1.d0/(2.d0*tm)*dot_product(cm(:,m),bbq(:,1)))/(1.d0+1.d0/(2.d0*tm))

      curf=-4.d0*(clm-dot_product(cm(:,l),bbq(:,2)))/tm+4.d0*(cml-dot_product(cm(:,m),bbq(:,1)))/tm

      end function

!-----------------------------

      real(8) function divb(i,j,k)
      integer ib,id,i,j,k
      real(8) dbq(0:6,3),d(3)

      do ib=0,6
      dbq(ib,1)=beq(ib,1)
      dbq(ib,2)=beq(ib,2)
      dbq(ib,3)=beq(ib,3)
      enddo

      do id=1,3
      d(id)=(dot_product(cm(:,id),ba(id,:,i,j,k))+1.d0/(2.d0*tm)*dot_product(cm(:,id),dbq(:,id)))/(1.d0+1.d0/(2.d0*tm))
      enddo

       divb=dabs(sum(d))*4.d0/tm

      end function

!--------------------------
!   TIME INTEGRATION
!--------------------------
      subroutine timeint

      do it=1,nt

        if(taskid.eq.0) then

          print*
          print*,'>>>>>>>>'
          print*,'it=',it,'time=',it*ts

        endif

      ncount=ncount+1

!-------------------------------
!   compute collision terms
!-------------------------------

      if(mode=='mrt') call mrt
      if(mode=='bgk') call collf

      call collb

!----------------
! Triple Periodic Condition
!----------------
       if(opt.eq.0) then
           do k=ks,kend

           do j=1,nj

                df(:,0,j,k)=df(:,ni-1,j,k)
                df(:,ni+1,j,k)=df(:,2,j,k)

                db(:,:,0,j,k)=db(:,:,ni-1,j,k)
                db(:,:,ni+1,j,k)=db(:,:,2,j,k)
           enddo

           do i=0,ni+1

                df(:,i,0,k)=df(:,i,nj-1,k)
                df(:,i,nj+1,k)=df(:,i,2,k)

                db(:,:,i,0,k)=db(:,:,i,nj-1,k)
                db(:,:,i,nj+1,k)=db(:,:,i,2,k)

           enddo

           enddo
         endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------
! Each node sends their boundary collision terms
! for the streaming process
!-------------------------------------------------
      if(numtasks.gt.1) call mpicom
!-----------------
! for serial run
!-----------------
      if(numtasks.eq.1.) then

        df(:,:,:,0)=df(:,:,:,nk-1)
        df(:,:,:,nk+1)=df(:,:,:,2)

        db(:,:,:,:,0)=db(:,:,:,:,nk-1)
        db(:,:,:,:,nk+1)=db(:,:,:,:,2)

      endif
!------------------------
!  streaming fa, ba
!------------------------

      call stream

!=========================================================
       if (opt.ne.0) call bound
!=========================================================

       call macro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------
!   UPDATE CURRENT/VORTICITY
!-----------------------------
!  Get the neighboring u
!--------------------------
      if(numtasks.gt.1) call update

       call cur

!=================
       if(ncount.eq.nframes) then
        call output
        ncount=0
       endif
!=================

!     if(dmax.lt.crit) exit

!--------------------------
!   END TIME INTEGRATION
!--------------------------
      enddo   ! next it

      end subroutine

!--------------------------------------
!  Taylor Green Vortex initialization
!--------------------------------------
      subroutine tg


      do i=1,ni
       do j=1,nj
        do k=ks,kend

        x1=0.d0+real(i-1)*dd
        x2=0.d0+real(j-1)*dd
        x3=0.d0+real(k-1)*dd

!      call innerb
!------------------------
!  initialize macroscopic quantities
!------------------------
        u(1,i,j,k)=u0*dsin(x1)&
                    *dcos(x2)*dcos(x3)
        u(2,i,j,k)=-u0*dcos(x1)&
                    *dsin(x2)*dcos(x3)
        u(3,i,j,k)=0.d0

        b(1,i,j,k)=0.d0
        b(2,i,j,k)=by
        b(3,i,j,k)=0.d0

        b(1,i,j,k)=a*u0*dsin(x1)&
                    *dsin(x2)*dcos(x3)
        b(2,i,j,k)=a*u0*dcos(x1)&
                    *dcos(x2)*dcos(x3)
        b(3,i,j,k)=0.d0



        dens(i,j,k)=dens0


        enddo
        enddo
        enddo

        u(:,0,:,ks:kend)=u(:,ni-1,:,ks:kend)
        u(:,:,0,ks:kend)=u(:,:,nj-1,ks:kend)
        u(:,ni+1,:,ks:kend)=u(:,2,:,ks:kend)
        u(:,:,nj+1,ks:kend)=u(:,:,2,ks:kend)

        b(:,0,:,ks:kend)=b(:,ni-1,:,ks:kend)
        b(:,:,0,ks:kend)=b(:,:,nj-1,ks:kend)
        b(:,ni+1,:,ks:kend)=b(:,2,:,ks:kend)
        b(:,:,nj+1,ks:kend)=b(:,:,2,ks:kend)


        end subroutine

!--------------------------------------
!  Orszag Tang Vortex initialization
!--------------------------------------
        subroutine ot

      do i=1,ni
       do j=1,nj
        do k=ks,kend

        x1=0.d0+real(i-1)*dd
        x2=0.d0+real(j-1)*dd
        x3=0.d0+real(k-1)*dd

!      call innerb
!------------------------
!  initialize macroscopic quantities
!------------------------

        u(1,i,j,k)=-2.d0*u0*dsin(x2)
        u(2,i,j,k)=2.d0*u0*dsin(x1)
        u(3,i,j,k)=0.d0

        b(1,i,j,k)=-2.d0*a*u0*dsin(2.d0*x2)+a*u0*dsin(x3)
        b(2,i,j,k)=2.d0*a*u0*dsin(x1)+a*u0*dsin(x3)
        b(3,i,j,k)=a*u0*dsin(x1)+a*u0*dsin(x2)

!       if(dsqrt(x1**2+x2**2)<=r0) u(3,i,j,k)=u0

        dens(i,j,k)=dens0


        enddo
        enddo
        enddo

        u(:,0,:,ks:kend)=u(:,ni-1,:,ks:kend)
        u(:,:,0,ks:kend)=u(:,:,nj-1,ks:kend)
        u(:,ni+1,:,ks:kend)=u(:,2,:,ks:kend)
        u(:,:,nj+1,ks:kend)=u(:,:,2,ks:kend)

        b(:,0,:,ks:kend)=b(:,ni-1,:,ks:kend)
        b(:,:,0,ks:kend)=b(:,:,nj-1,ks:kend)
        b(:,ni+1,:,ks:kend)=b(:,2,:,ks:kend)
        b(:,:,nj+1,ks:kend)=b(:,:,2,ks:kend)

        end subroutine

        subroutine setup

!-----------------------
! DOMAIN DECOMPOSITION
!-----------------------

      ksize=nk/numtasks
      ks=ksize*taskid+1
      Kend=ksize*(taskid+1)

!     print*, taskid,ks,kend
!-----------------------
! MEMORY ALLOCATION
!-----------------------
      allocate (fa(0:nbit-1,1:ni,1:nj,ks:kend))
      allocate(ba(1:3,0:6,1:ni,1:nj,ks:kend))
      allocate(db(1:3,0:6,0:(ni+1),0:(nj+1),(ks-1):(kend+1)))
      allocate(df(0:nbit-1,0:(ni+1),0:(nj+1),(ks-1):(kend+1)))
      allocate(u(1:3,0:(ni+1),0:(nj+1),(ks-1):(kend+1)))
      allocate(b(1:3,0:(ni+1),0:(nj+1),(ks-1):(kend+1)))
      allocate(dens(1:ni,1:nj,ks:kend))
      allocate(ctot(1:ni,1:nj,ks:kend))
      allocate(vtot(1:ni,1:nj,ks:kend))
      allocate(fm(1:ni,1:nj,ks:kend,3))

      allocate(ms(0:nbit-1,0:nbit-1),e(0:nbit-1,3),lb(1:nbit-1),meq(0:nbit-1))
      allocate(invms(0:nbit-1,0:nbit-1),c(0:nbit-1,3),c_(0:nbit-1))
      allocate(phi(1:12),s(0:nbit-1))

! Lattice

      if(nbit==19)  call bit19
      if(nbit==13)  call bit13

!--------------------------
!   SET VARIABLES
!--------------------------
       c=real(e)
       cm=real(em)

!------------------------
!  coefficients
!------------------------
       w(0)=1.d0/3.d0
       w(1:6)=1.d0/18.d0
       w(7:18)=1.d0/36.d0

       wm(0)=1.d0/4.d0
       wm(1:6)=1.d0/8.d0

!------------------------
!  Value of |c|
!------------------------
      c_(:)=dsqrt(c(:,1)**2+c(:,2)**2+c(:,3)**2)


      end subroutine


      subroutine bit19

!----------------------------
!  directional vectors 19bit
!----------------------------

      e(0,:)=(/0,0,0/)
      e(1,:)=(/1,0,0/)
      e(2,:)=(/-1,0,0/)
      e(3,:)=(/0,1,0/)
      e(4,:)=(/0,-1,0/)
      e(5,:)=(/0,0,1/)
      e(6,:)=(/0,0,-1/)
      e(7,:)=(/1,1,0/)
      e(8,:)=(/-1,1,0/)
      e(9,:)=(/1,-1,0/)
      e(10,:)=(/-1,-1,0/)
      e(11,:)=(/1,0,1/)
      e(12,:)=(/-1,0,1/)
      e(13,:)=(/1,0,-1/)
      e(14,:)=(/-1,0,-1/)
      e(15,:)=(/0,1,1/)
      e(16,:)=(/0,-1,1/)
      e(17,:)=(/0,1,-1/)
      e(18,:)=(/0,-1,-1/)
!------------------------
      lb(:)=(/2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15/)
!------------------------
! Moment transformation matrix
!-----------------------------

      ms(0,:)=(/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
      ms(1,:)=(/-30.d0,-11.d0,-11.d0,-11.d0,-11.d0,-11.d0,-11.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0,8.d0/)
      ms(2,:)=(/12.d0,-4.d0,-4.d0,-4.d0,-4.d0,-4.d0,-4.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
      ms(3,:)=(/0.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(4,:)=(/0.d0,-4.d0,4.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(5,:)=(/0.d0,0.d0,0.d0,1.d0,-1.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0/)
      ms(6,:)=(/0.d0,0.d0,0.d0,-4.d0,4.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0/)
      ms(7,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0/)
      ms(8,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,-4.d0,4.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0/)
      ms(9,:)=(/0.d0,2.d0,2.d0,-1.d0,-1.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,-2.d0,-2.d0,-2.d0,-2.d0/)
      ms(10,:)=(/0.d0,-4.d0,-4.d0,2.d0,2.d0,2.d0,2.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,-2.d0,-2.d0,-2.d0,-2.d0/)
      ms(11,:)=(/0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,1.d0,-1.d0,-1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(12,:)=(/0.d0,0.d0,0.d0,-2.d0,-2.d0,2.d0,2.d0,1.d0,1.d0,1.d0,1.d0,-1.d0,-1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(13,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(14,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0/)
      ms(15,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(16,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(17,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,-1.d0,1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,1.d0,-1.d0/)
      ms(18,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,-1.d0,-1.d0,1.d0,1.d0/)

!==========================
! collision operators
!==========================
! lbgk
      we=3.d0
      wej=-11.d0/2.d0
      wxx=-1.d0/2.d0
! optimum
      we=0.d0
      wej=-475.d0/63.d0
      wxx=0.d0

!==========================
      s(0)=0.d0

!     s(1)=1.d0/(9.d0/2.d0*bvisc+0.5d0)
      s(1)=1.19d0
      s(2)=1.4d0

      s(3)=0.d0
      s(4)=1.2d0
      s(5)=0.d0
      s(6)=1.2d0
      s(7)=0.d0
      s(8)=1.2d0

      s(9)=1.d0/(3.d0*visc+0.5d0)
      s(10)=1.4d0
      s(11)=1.d0/(3.d0*visc+0.5d0)
      s(12)=1.4d0

      s(13)=1.d0/(3.d0*visc+0.5d0)
      s(14)=1.d0/(3.d0*visc+0.5d0)
      s(15)=1.d0/(3.d0*visc+0.5d0)

      s(16)=1.98d0
      s(17)=1.98d0
      s(18)=1.98d0

!     s(:)=1.d0/(visc+0.5d0)


      end subroutine


      subroutine bit13

!----------------------------
!  directional vectors 13bit
!----------------------------

      e(0,:)=(/0,0,0/)
      e(1,:)=(/1,1,0/)
      e(2,:)=(/-1,1,0/)
      e(3,:)=(/1,-1,0/)
      e(4,:)=(/-1,-1,0/)
      e(5,:)=(/0,1,1/)
      e(6,:)=(/0,-1,1/)
      e(7,:)=(/0,1,-1/)
      e(8,:)=(/0,-1,-1/)
      e(9,:)=(/1,0,1/)
      e(10,:)=(/1,0,-1/)
      e(11,:)=(/-1,0,1/)
      e(12,:)=(/-1,0,-1/)
!------------------------
      lb(:)=(/4,3,2,1,8,7,6,5,12,11,10,9/)
!-----------------------------
! Moment transformation matrix
!-----------------------------

      ms(0,:)=(/1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
      ms(1,:)=(/0.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0/)
      ms(2,:)=(/0.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(3,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0/)
      ms(4,:)=(/-12.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
      ms(5,:)=(/0.d0,1.d0,1.d0,1.d0,1.d0,-2.d0,-2.d0,-2.d0,-2.d0,1.d0,1.d0,1.d0,1.d0/)
      ms(6,:)=(/0.d0,1.d0,1.d0,1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,-1.d0,-1.d0,-1.d0/)
      ms(7,:)=(/0.d0,1.d0,-1.d0,-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(8,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(9,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0/)
      ms(10,:)=(/0.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,-1.d0,1.d0,1.d0/)
      ms(11,:)=(/0.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,-1.d0,1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0/)
      ms(12,:)=(/0.d0,0.d0,0.d0,0.d0,0.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,-1.d0,1.d0,-1.d0/)

!     a2=0.d0

      s(0)=0.d0
      s(1)=0.d0
      s(2)=0.d0
      s(3)=0.d0
      s(4)=1.5d0
      s(5)=1.d0/(4.d0*visc+.5d0)
      s(6)=1.d0/(4.d0*visc+.5d0)
      s(7)=1.d0/(2.d0*visc+.5d0)
      s(8)=1.d0/(2.d0*visc+.5d0)
      s(9)=1.d0/(2.d0*visc+.5d0)
      s(10)=1.8d0
      s(11)=1.8d0
      s(12)=1.8d0

      end subroutine

      subroutine meq19


          meq(0)=dens(i,j,k)
          meq(1)=-11.d0*dens(i,j,k)+19.d0*u2
          meq(2)=we*dens(i,j,k)+wej*u2

          meq(3)=u(1,i,j,k)
          meq(4)=-2.d0/3.d0*u(1,i,j,k)
          meq(5)=u(2,i,j,k)
          meq(6)=-2.d0/3.d0*u(2,i,j,k)
          meq(7)=u(3,i,j,k)
          meq(8)=-2.d0/3.d0*u(3,i,j,k)

          meq(9)=1.d0/3.d0*(3.d0*u(1,i,j,k)**2-u2)
          meq(10)=wxx*meq(9)
          meq(11)=u(2,i,j,k)**2-u(3,i,j,k)**2
          meq(12)=wxx*meq(11)

          meq(13)=u(1,i,j,k)*u(2,i,j,k)
          meq(14)=u(2,i,j,k)*u(3,i,j,k)
          meq(15)=u(1,i,j,k)*u(3,i,j,k)

          meq(16)=0.d0
          meq(17)=0.d0
          meq(18)=0.d0


      end subroutine

!-----------------------------
!   compute meq with new u,b
!-----------------------------
      subroutine meq13

          meq(0)=dens(i,j,k)
          meq(1)=u(1,i,j,k)
          meq(2)=u(2,i,j,k)
          meq(3)=u(3,i,j,k)
          meq(4)=3.d0/2.d0*(13.*cs02-8.)*dens(i,j,k)+13.d0/2.d0*u2
!         meq(4)=3.d0/2.d0*(13.d0*cs02-8.d0)*dens(i,j,k)+13.d0/(2.d0*dens(i,j,k))*u2
          meq(5)=3.d0*u(1,i,j,k)**2-u2
!       meq(5)=meq(5)/dens(i,j,k)
          meq(6)=u(2,i,j,k)**2-u(3,i,j,k)**2
!       meq(6)=meq(6)/dens(i,j,k)
          meq(7)=u(1,i,j,k)*u(2,i,j,k)
!       meq(7)=meq(7)/dens(i,j,k)
          meq(8)=u(2,i,j,k)*u(3,i,j,k)
!       meq(8)=meq(8)/dens(i,j,k)
          meq(9)=u(1,i,j,k)*u(3,i,j,k)
!       meq(9)=meq(9)/dens(i,j,k)
          meq(10)=0.
          meq(11)=0.
          meq(12)=0.

      end subroutine

!-----------------------------

      subroutine vorf
      integer l1,l2,lk
      real(8) p(3,3),ffb(0:nbit-1)

      do lk=0,nbit-1
       ffb(lk)=feq(lk)
      enddo

!   bbq => l==1, m==2
      do l1=1,3
      do l2=1,3

      p(l1,l2)=(sum(c(:,l1)*c(:,l2)*fa(:,i,j,k))-sum(c(:,l1)*c(:,l2)*ffb(:)))/(1.d0+1.d0/(2.d0*tu))

      enddo
      enddo

      p=-p*3.d0/(dens(i,j,k)*tu)

!     write(55,*) p(1,1)+p(2,2)+p(3,3)



      vtot(i,j,k)=sum(p**2)/2.d0
!         -(p(1,1)+p(2,2)+p(3,3))/2.d0

      end subroutine

!--------------------------
!   END OF MODULE
!--------------------------
      end
