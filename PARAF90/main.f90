
module basicmod
  implicit none
  integer::nhy
  integer,parameter::nhymax=2000
  integer,parameter::ngridx=90
  integer,parameter::ngridy=180
  integer,parameter::ngridz=1
  integer,parameter::in=ngridx &
 &                  ,jn=ngridy &
 &                  ,kn=ngridz
  integer,parameter::is=1 &
 &                  ,js=1 &
 &                  ,ks=1
  integer,parameter::ie=ngridx &
 &                  ,je=ngridy &
 &                  ,ke=ngridz

  real(8),dimension(in,jn,kn):: A,B,S
  
end module basicmod

subroutine setmpi
  use mpimod, only: ntiles
  implicit none
  ntiles(1) = 4
  ntiles(2) = 2
  ntiles(3) = 1
end subroutine setmpi

program main
  use basicmod
  use omp_lib
  use mpimod
  implicit none
  real(8)::time_begin,time_end
  integer::threadsnum
  logical,parameter::nooutput=.true.
  logical,parameter::debug=.false.
  
  call InitializeMPI
  if(myid_w == 0) print *, "grid size for x y",ngridx*ntiles(1),ngridy*ntiles(2)

  threadsnum = omp_get_max_threads()
  if(myid_w == 0) print *, "threads=",threadsnum

  call InitialzeVariables
  
  if(myid_w == 0) print *, "entering main loop"
! main loop
  time_begin = omp_get_wtime()
  do nhy=1,nhymax
     call ParaOpenmp
  enddo
  time_end = omp_get_wtime()
  
  if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
  if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngridx*ntiles(1)*ngridy*ntiles(2))/nhymax

  time_begin = omp_get_wtime()
  do nhy=1,nhymax
     call ParaConcurent
  enddo
  time_end = omp_get_wtime()
  
  if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
  if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngridx*ntiles(1)*ngridy*ntiles(2))/nhymax

  call FinalizeMPI
  if(myid_w == 0) print *, "program has been finished"
end program main

subroutine InitialzeVariables
  use mpimod, only:myid
  use basicmod
  implicit none
  integer i,j,k
  do k=ks,ke
  do j=js,je
  do i=is,ie
!     S(i,j,k) =  myid*in*jn*kn + in*jn*(k-1)+ in*(j-1) + i 
     A(i,j,k) =  1.0d0
     B(i,j,k) =  2.0d0
     S(i,j,k) =  0.0d0
  enddo
  enddo
  enddo
end subroutine InitialzeVariables


subroutine ParaOpenmp
  use mpimod
  use basicmod
  use omp_lib
  implicit none
  integer i,j,k
  
!$omp parallel
!$omp do collapse(3)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     S(i,j,k) = A(i,j,k) + B(i,j,k)  
  enddo
  enddo
  enddo
!$omp end do
!$omp end parallel
end subroutine ParaOpenmp


subroutine ParaConcurent
  use mpimod
  use basicmod
  use omp_lib
  implicit none
  integer i,j,k
  
  do concurrent (i=is:ie, j=js:je,k=ks:ke)
     S(i,j,k) = A(i,j,k) + B(i,j,k)  
  enddo
end subroutine ParaConcurent
