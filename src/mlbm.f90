!===================================================================
! 3D Lattice Kinetic Code utilizing MRT approximation
! Main
! MPI Port
! Author: George Breyiannis, Dec 2004 
!===================================================================
      program lkc

      use lkcmod

      implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------
!  MASTER TASKS
!-----------------------------

      if(taskid.eq.0) then
      open(1,file='input')

      open(12,file='stat')
      open(13,file='statc')
      open(14,file='statv')
      open(15,file='corel')
      open(30,file='disr')
      open(31,file='ene')

      open(2,file='uv')
      open(3,file='de')
      open(7,file='bi')
      open(9,file='c')
      open(10,file='v')
      open(11,file='d')

!--------------------------
!   READ DATA    
!--------------------------
      read(1,*) ni,nj,nk,time
      read(1,*) visc,kel
      read(1,*) dens0,u0
      read(1,*) Ha,a
      read(1,*) crit
      read(1,*) relax
      read(1,*) nframes
      read(1,*) mode
      read(1,*) nbit
      read(1,*) problem
      read(1,*) bgk_option
      read(1,*) magn
! triple periodic
      read(1,*) opt
      read(1,*) r0

       endif

      call MPI_BCAST(ni,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(opt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(visc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kel,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Ha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(u0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dens0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(crit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(nframes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(mode,3,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(problem,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(magn,2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bgk_option,4,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(relax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(r0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!------------------------
! Various properties
!------------------------
      pi=4.d0*datan(1.d0)

      length=2.d0*pi

      cs02=1.d0/3.d0

!     u0=u0*dsqrt(cs02)
      if (taskid.eq.0) then
            print*, 'Ma=',u0/dsqrt(cs02)
      endif

 
      dd=length/real(nj-1)
      h=real(nj-1)/2.d0

      tm=kel*4.d0
      tu=visc*3.d0

      ts=length/(nj-1)*u0

      nt=time/ts

      if (taskid.eq.0) then
            print*, 'nt=',nt
      endif


      if(nframes==0) then
      nframes=nt+1
      else
      nframes=nt/nframes
      endif


      n1=nj+2

! Allocate/Setip

      call setup

!mrt

      invms=matmul(ms,transpose(ms))
      
      forall(i=0:nbit-1,j=0:nbit-1,i==j)
       invms(i,j)=1.d0/invms(i,j)
      endforall

      invms=matmul(transpose(ms),invms)

!------------------------
!  set problem     
!------------------------


      if(problem=='tg') call tg
      if(problem=='ot') call ot


!---------------------------------
!  initialize local quantities 
!---------------------------------
      if(numtasks.gt.1) call update

      call init

      call cur

!--------------------------
!   TIME INTEGRATION
!--------------------------
      call timeint

      if(nt.eq.0.and.nframes<nt) call output

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call MPI_FINALIZE(ierr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------
!   END OF PROGRAM   
!--------------------------
      end
