! compile in fortran 
!  ifort f90_master.f90 -I. libff-mmap-semaphore.o -o f90_master
! launch
!   ./f90_master
! warning path of freefem++ command and of ffslave.edp  is in hard. 
   
program main

  integer*8            :: sem_ff, sem_fo, shd, status, ret
  double precision     :: cff, rff
  integer :: i

  call ffsem_init(sem_ff,'ff-slave1'//achar(0),1)
  call ffsem_init(sem_fo,'ff-master1'//achar(0),1)
  call ffmmap_init(shd,'shared-data'//achar(0),1024)

  status=1
  call system('../../src/nw/FreeFem++ ../../examples/plugin/ffslave.edp -nw -ns  &')
! here warning path of freefem++ command and of ffslave.edp  is in hard. 

  call ffmmap_write(shd, status, 8, 8,ret)
  call ffmmap_msync(shd,0,32,ret)
  call ffsem_wait(sem_ff,ret)
 
  do i=1,10
     cff=10+i
     call ffmmap_write(shd, cff, 8, 0,ret)
     write(*,*) 'bb ffsem_post'
     call ffsem_post(sem_fo,ret)
     write(*,*) 'bb ffsem_wait'
     call ffsem_wait(sem_ff,ret)
     call ffmmap_read(shd, rff, 8, 16,ret)
     write(*,*) 'iter',i,rff
     
  enddo

  status =0
  call ffmmap_write(shd, status, 8, 8,ret)
  call ffsem_post(sem_fo,ret)
  write(*,*) 'Fin Master'
  call ffsem_wait(sem_ff,ret)
  write(*,*) 'Fin Freefem dans master'
  call ffsem_del(sem_ff)
  call ffsem_del(sem_fo)
  call ffmmap_del(shd)
     


end
