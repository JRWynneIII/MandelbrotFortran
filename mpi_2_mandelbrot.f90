
!*****************************************************************************80
!
!! MPI_2_Mandelbrot produces an image with the Mandelbrot pattern.
!
!  Discussion:
!
!   This project creates an image of a Mandelbrot fractal 
!   (http://en.wikipedia.org/wiki/Mandelbrot_set). This code is a hybrid 
!   OpenMP and MPI version.
!
!   This code uses the Distance Estimator 
!   (http://mrob.com/pub/muency/distanceestimator.html) method to calculate 
!   the Mandelbrot set. Each pixel requires its own independent unit of work. 
!   Because these units of work are independent, it is a prime candidate 
!   for parallelizaion. 
!
!   The code will perform an iterative calculation at each point that will 
!   determine weather the point escapes the set 
!   (http://en.wikipedia.org/wiki/Mandelbrot_set#Formal_definition). If the 
!   point doesn't escape, then it is considered part of the Mandelbrot set 
!   and the pixel is painted black.
!
!   The Distance Estimator comes into play if the point escapes very slowly 
!   (takes many iterations before it escapes). For visualization purposes, we 
!   consider these points to be part of the set. By doing this, it will 
!   "reveal" more of the set than the standard method. 
!
!   At the beginning of the code, a 2 dimensional grid is allocated 
!   in memory. This grid will hold either a 1 or 0 in each element to 
!   represent whether a point is in the set or not. 
!
!   The code then iterates over each point in the grid and performs the 
!   Distance Estimator calculation. When the calculation determines 
!   if the point escapes, it breaks out of the loop and writes a "0" in the 
!   grid for that point. Else, it will infer that the point does not escape 
!   and will write a "1" in the grid. 
! 
!   After iterating over each point, the grid is then passed to a function 
!   that will write out a pgm image to the current working directory. The 
!   function is code by John Burkardt and can be obtained at 
!   http://people.sc.fsu.edu/~jburkardt/f_src/pgma_io/pgma_io.html
!
!   This version uses a dynamic decomposition of the domain. Small portions
!   of the entire domain are dispatched to MPI processes to compute. When
!   an process finishes its computation, it requests another portion of the 
!   the domain to compute.
!
!   To compile the code use your fortran compiler with MPI wrapper, 
!   for example gfortran or OpenMPI/mpich/mvapich and gfortran
!
!   mpif90 -fopenmp mpi_2_mandelbrot.f90 -o mpi_2_mandelbrot
!
!   The code can then be run
!
!   mpirun mpi_2_mandelbrot
!
!   and once it completes, should put a picture file in the running directory.
!   The picture can be opened in many viewers, such as GIMP 
!   (https://www.gimp.org/) or Gwenview (https://userbase.kde.org/Gwenview)
!
!  Licensing:
!
!     
!
!  Modified:
!
!    3 January 2016
!
!  Author:
!
!    James Riley Wynne III
!
!  Last Modified By:
!
!    Benson Muite
!
! IN PROGRESS

program Mandelbrot
      use mpi
      implicit none

      ! x and y resolution
      integer(kind=4), parameter :: nx = 1000
      integer(kind=4), parameter :: ny = 1000
      ! Maximum number of iterations
      integer(kind=4), parameter :: maxiter = 2000
      ! Number of points to process in each chunk given to an MPI task
      ! this should be smaller than nx * ny
      integer(kind=4), parameter :: processchunk = 100
      ! Use a two dimensional array to hold data for the global image
      integer(kind=4), dimension(:,:), allocatable :: MSet
      ! Use a one dimensional array to hold data for the global image
      integer(kind=4), dimension(:), allocatable :: myMSet
      ! Change maximum and minimum values of the window that the
      ! appear in the image. Increase the values and you zoom out
      ! decrease the values and you zoom in
      integer(kind=4), parameter :: xmin = -2
      integer(kind=4), parameter :: ymin = -2
      integer(kind=4), parameter :: xmax = 2
      integer(kind=4), parameter :: ymax = 2
      integer(kind=4) :: ix, iy, ii, proc
      integer(kind=4) :: iter, i, ierror, allocatestatus
      ! name of file to write image out to
      character ( len = 80 ) :: file_name = 'mandelbrot.ascii.pgm'
      ! MPI variables
      integer(kind=4) :: ierr, myid, numprocs, mystart, myend
      integer(kind=4) :: itemcount, procstart, procend
      integer(kind=4) :: mpistat(MPI_STATUS_SIZE) 
      delta = (threshold*(xmax-xmin))/real(nx-1)
      ! initialisation of MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      if (myid.eq.0) then
        print *,"Program to calculate Mandelbrot set started"
      end if
      mystart=1+floor(myid*Nx*Ny/real(numprocs))
      myend=min(floor((myid+1.0d0)*Nx*Ny/real(numprocs)),Nx*Ny)
      allocate(MSet(1:nx,1:ny),myMSet(1:processchunk),stat=allocatestatus)
      if (allocatestatus .ne. 0) stop
      MSet(1:nx,1:ny)=0
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (myid.eq.0) then
        print *,"Memory allocated"
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (myid.eq.0) then
        print *,"Mandelbrot set calculated on each process"
      end if

      ! Get Process 0 to collect data and write out a file. This is scalable
      ! to a moderate number of processes, limited by memory on process 0.
      if (myid.eq.0) then

        do ii=1,nx*ny,processchunk


        end do
        do ii=mystart,myend
          iy=1+floor((ii-1)/real(Nx))
          ix=ii-(iy-1)*Nx
          MSet(ix,iy)=myMSet(ix+Nx*(iy-1))
        end do

        if (numprocs.gt.1) then
          do proc=1,numprocs-1
            procstart=1+floor(proc*Nx*Ny/real(numprocs))
            procend=min(floor((proc+1.0d0)*Nx*Ny/real(numprocs)),Nx*Ny)
            itemcount=1+procend-procstart
            call MPI_RECV(myMSet,itemcount,MPI_INTEGER,proc,proc, &
                          MPI_COMM_WORLD,mpistat,ierr)
            ! copy data into main array
            do ii=procstart,procend
              iy=1+floor((ii-1)/real(Nx))
              ix=ii-(iy-1)*Nx
              MSet(ix,iy)=myMSet(mystart+ix+Nx*(iy-1)-procstart)
            end do
          end do
        end if
      else
         ! other processes send data to process 0
         itemcount=1+myend-mystart
         call MPI_SEND(myMSet,itemcount,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (myid.eq.0) then
        print *,"Mandelbrot set collected on one process"
      end if

      if (myid.eq.0) then
        ! Call Burkardt's Fortran code to write Mandelbrot picture to disk
        call pgma_write ( file_name, nx, ny, Mset, ierror )
        print *,"Saved image file to disk"
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      deallocate(MSet,myMSet,stat=allocatestatus)
      if (allocatestatus .ne. 0) stop
      call MPI_FINALIZE(ierr)
end program Mandelbrot


subroutine calcset(nx,ny,mystart,myend,processchunk,xmax,xmin,ymax,ymin,myMSet)

      ! We use a nested loop here to effectively traverse over each part of the
      ! grid (pixel of the image) in sequence. First, the complex values of the
      ! points are determined and then used as the basis of the computaion. 
      ! Effectively, it will loop over each point (pixel) and according on how
      ! many iterations it takes for the value that the mathematical function 
      ! returns on each iteration it will determine whether or not the point 
      ! "escapes" to infinity (or an arbitrarily large number.) or not. If it 
      ! takes few iterations to escape then it will decide that this point is
      ! NOT part of the Mandelbrot set and will put a 0 in that point's index
      ! in MSet. If it takes nearly all or all of the iterations to escape, 
      ! then it will decide that the point/pixel is part of the Mandelbrot set
      ! and instead put a 1 in its place in myMSet.
      !
      ! The use of the OpenMP pragma here will divide up the iterations between
      ! threads and execute them in parallel.
      !
      ! This region is VERY easily parallelized because there is NO data shared
      ! between the loop iterations.


      use omp_lib
      implicit none

      integer(kind=4), intent(in) :: nx, ny, maxiter, mystart, myend 
      integer(kind=4), intent(in) :: processchunk 
      real(kind=8), intent(in) :: xmax, xmin, ymax, ymin
      integer(kind=8), intent(out), dimension(1:processchunk) :: myMset
      real(kind=8) :: threshold = 1
      real(kind=8) :: dist = 0
      integer(kind=4) :: ix, iy, ii, proc
      real(kind=8) :: cx, cy
      integer(kind=4) :: iter, i, ierror, allocatestatus
      real(kind=8) :: x,y,x2,y2
      real(kind=8) :: temp = 0.0
      real(kind=8) :: xder = 0.0
      real(kind=8) :: yder = 0.0
      real(kind=8), dimension(1:maxiter) :: xorbit
      real(kind=8), dimension(1:maxiter) :: yorbit
      integer(kind=4), parameter :: hugenum = 100000
      logical :: flag = .false.
      integer(kind=4), parameter :: overflow = 2147483647
      real(kind=8) :: delta 
      ! chunk size for dynamic scheduling of shared memory loop
      integer(kind=4), parameter :: chunk = 1 

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ii,ix,iy,cx,cy,iter,i,x,y,x2,y2) &
      !$OMP PRIVATE(temp,xder,yder,dist,xorbit,yorbit,flag) &
      !$OMP ORDERED SHARED(myMSet) SCHEDULE(DYNAMIC,chunk)
      do ii=mystart,myend
        iy=1+floor((ii-1)/real(nx))
        ix=ii-(iy-1)*nx
        cy = ymin+(iy-1)*(ymax-ymin)/real(ny-1) 
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
        cx = xmin + (ix-1)*(xmax-xmin)/real(nx-1) 
        ! This is the main loop that determins whether or not the point escapes
        ! or not. It breaks out of the loop when it escapes
        do iter = 0,maxiter-1
          temp = x2-y2 + cx
          y = 2.0*x*y+cy
          x = temp
          x2 = x**2
          y2 = y**2
          xorbit(iter+1) = x
          yorbit(iter+1) = y
          ! if point escapes then break to next loop
          if (x2+y2 .gt. hugenum) exit
        enddo
        ! if the point escapes, find the distance from the set, just in case 
        ! its close to the set. if it is, it will make it part of the set.

        if (x2+y2 .ge. hugenum) then
          xder = 0
          yder = 0
          i = 0
          flag = .false.

          do i=1, iter
                
            if (flag .eqv. .true.) exit
                
            temp = 2*(xorbit(i)*xder-yorbit(i)*yder)+1
            yder = 2*(yorbit(i)*xder+xorbit(i)*yder)
            xder = temp
            flag = max(abs(xder),abs(yder)) > overflow
          enddo

          if (flag .eqv. .false.) then
            dist = (log(x2+y2)*sqrt(x2+y2))/sqrt(xder**2+yder**2)
          endif
        endif
        ! Assign the appropriate values to MSet in the place relating to 
        ! the point in question
        if (dist .lt. delta) then
          myMSet(ix+Nx*(iy-1)-mystart+1) = 1
        else
          myMSet(ix+Nx*(iy-1)-mystart+1) = 0
        endif
      enddo  
      !$OMP END PARALLEL DO
end subroutine calcset

subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end

subroutine pgma_write ( file_out_name, row_num, col_num, g, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE writes an ASCII PGM file.
!
!  Example:
!
!    P2
!    # feep.pgm
!    24 7
!    15
!    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
!    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
!    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
!    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
!    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
!    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
!    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    2 January 2017
!
!  Original Source Code:
!
!     http://people.sc.fsu.edu/%7Ejburkardt/f_src/pgma_io/pgma_io.f90
!
!  Author:
!
!    John Burkardt
!
!  Last Modified by:
! 
!    Benson Muite
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), the gray value of each 
!    pixel.  These should be nonnegative.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: col_num
  integer ( kind = 4 ), intent(in) :: row_num
  character ( len = * ), intent(in) ::  file_out_name
  integer ( kind = 4 ), dimension(1:row_num,1:col_num), intent(in) :: g
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) g_max

  ierror = 0
!
!  Compute the maximum color value.
!
  g_max = maxval ( g(1:row_num,1:col_num) )
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    stop
  end if
!
!  Write the header.
!
  call pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
    g_max, ierror )
!
!  Write the data.
!
  call pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE - Note:'
    write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
    write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
    write ( *, '(a,i8)' ) '  Maximum gray value G_MAX =       ', g_max
  end if

  return
end
subroutine pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE_DATA writes the data of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    2 January 2017
!
!  Original Source Code:
!
!     http://people.sc.fsu.edu/%7Ejburkardt/f_src/pgma_io/pgma_io.f90
!
!  Author:
!
!    John Burkardt
!
!  Last Modified by:
! 
!    Benson Muite
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), the gray value of each 
!    pixel.  These should be nonnegative.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ), intent(in) :: col_num
  integer ( kind = 4 ), intent(in) :: row_num
  integer ( kind = 4 ), dimension(1:row_num,1:col_num), intent(in) ::  g
  integer ( kind = 4 ), intent(in) :: file_out_unit
  integer ( kind = 4 ), intent(out) :: ierror
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo

  ierror = 0

  do i = 1, row_num
    do jlo = 1, col_num, 12
      jhi = min ( jlo + 11, col_num )
      write ( file_out_unit, '(12i5)' ) g(i,jlo:jhi)
    end do
  end do

  return
end
subroutine pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
  g_max, ierror )

!*****************************************************************************80
!
!! PGMA_WRITE_HEADER writes the header of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!
!  Modified:
!
!    2 January 2017
!
!  Original Source Code:
!
!     http://people.sc.fsu.edu/%7Ejburkardt/f_src/pgma_io/pgma_io.f90
!
!  Author:
!
!    John Burkardt
!
!  Last Modified by:
! 
!    Benson Muite
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer ( kind = 4 ) G_MAX, the maximum gray value.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * ), intent(in) :: file_out_name
  integer ( kind = 4 ), intent(in) :: file_out_unit
  integer ( kind = 4 ), intent(in) :: g_max
  integer ( kind = 4 ), intent(in) :: col_num
  integer ( kind = 4 ), intent(in) :: row_num
  integer ( kind = 4 ), intent(out) :: ierror
  character ( len = 2 ) :: magic = 'P2'

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PGMA_IO::PGMA_WRITE.'
  write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num
  write ( file_out_unit, '(i8)' ) g_max

  return
end
