program Mandelbrot
      implicit none

      integer :: nx = 1000
      integer :: ny = 1000
      integer :: maxiter = 2000
      integer,dimension(:),allocatable :: MSet
      integer :: xmin = -3
      integer :: ymin = -2
      integer :: xmax = 1
      integer :: ymax = 2
      real(8) :: threshold = 1
      real(8) :: dist = 0
      integer :: ix, iy
      real(8) :: cx, cy
      integer :: iter, i
      real(8) :: x,y,x2,y2
      real(8) :: temp = 0.0
      real(8) :: xder = 0.0
      real(8) :: yder = 0.0
      real(8),dimension(:),allocatable :: xorbit
      real(8),dimension(:),allocatable :: yorbit
      integer :: hugenum = 100000
      logical :: flag = .false.
      integer,parameter :: overflow = 2147483647
      real(8) :: delta 
      delta = (threshold*(xmax-xmin))/real((nx-1))

      allocate(MSet(nx*ny))
      allocate(xorbit(maxiter))
      allocate(yorbit(maxiter))

      do iy=1,(ny-1)
        cy = (ymin+iy*(ymax-ymin))/real((ny-1))
        do ix=1,(nx-1)
            iter = 0
            i = 0
            x = 0
            y = 0
            x2 = 0
            y2 = 0
            temp = 0
            xder = 0
            yder = 0
            dist = 0
            cx = (xmin + ix*(xmax-xmin))/real(ny-1)
            do iter = 0,maxiter
              temp = x2-y2 + cx
              y = 2.0*x*y+cy
              x = temp
              x2 = x**2
              y2 = y**2
              xorbit(iter+1) = x
              yorbit(iter+1) = y
              if (x2+y2>hugenum) then 
                exit
              endif
            enddo

            if (x2+y2>=hugenum) then
              xder = 0
              yder = 0
              i = 0
              flag = .false.

              do i=0, iter
                
                if (flag == .true.) then
                  exit
                endif
                
                temp = 2*(xorbit(i)*xder-yorbit(i)*yder)+1
                yder = 2*(yorbit(i)*xder+xorbit(i)*yder)
                xder = temp
                flag = max(abs(xder),abs(yder)) > overflow
              enddo

              if (flag == .false.) then
                dist = (log(x2+y2)*sqrt(x2+y2))/sqrt(xder**2+yder**2)
              endif
            endif

            if (dist < delta) then
              MSet(iy*ny+ix) = 1
            else
              MSet(iy*ny+ix) = 0
            endif
        enddo  
      enddo

      do i = 1, (nx*ny)
        print*, "MSet[",i,"]: ",MSet(i)
      enddo

      ! Call my C function for calc_pixel_value(nx,ny,MSet,maxiter)

      deallocate(MSet)
end program Mandelbrot
