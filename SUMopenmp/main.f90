
module basicmod
  implicit none
  integer::nhy
  integer,parameter::nhymax=1
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

  real(8),dimension(in,jn,kn):: S
  real(8):: sumall
  
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
  threadsnum = omp_get_max_threads()
  if(myid_w == 0) print *, "threads=",threadsnum
  if(myid_w == 0) print *, "grid size for x y",ngridx*ntiles(1),ngridy*ntiles(2)

  call InitialzeVariables
  
  if(myid_w == 0) print *, "entering main loop"
! main loop
  time_begin = omp_get_wtime()
  do nhy=1,nhymax
     call SumVariables
  enddo
  
  if(myid_w == 0) print *, "sum",sumall
  time_end = omp_get_wtime()
  
  if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
  if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngridx*ngridy)/nhymax
  
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
     S(i,j,k) =  1.0d0
  enddo
  enddo
  enddo
end subroutine InitialzeVariables


subroutine SumVariables
  use mpimod
  use basicmod
  use omp_lib
  implicit none
  integer i,j,k
  real(8)::sum
  real(8)::sumsend, sumrecv
  
  sum = 0.0d0  
!$omp parallel reduction(+:sum)
!$omp do collapse(3)
  do k=ks,ke
  do j=js,je
  do i=is,ie
     sum = sum + S(i,j,k) 
  enddo
  enddo
  enddo
!$omp end do
  print *,"p1",myid,sum
!$omp end parallel
  sumsend = sum

  print *,"p2",myid,sumsend
  
  call MPI_ALLREDUCE( sumsend, sumrecv, 1  &
  &                 , MPI_DOUBLE_PRECISION &
  &                 , MPI_SUM, comm3d, ierr)

  sumall = sumrecv
  
  if(myid==0) print *,"p3",sumall
end subroutine SumVariables
