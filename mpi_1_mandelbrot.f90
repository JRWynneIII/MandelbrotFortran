program Mandelbrot
      use mpi
      use omp_lib
      implicit none

      ! x and y resolution
      integer(kind=4) :: nx = 1000
      integer(kind=4) :: ny = 1000
      ! Maximum number of iterations
      integer(kind=4) :: maxiter = 2000
      ! Use a two dimensional array to hold data for the global image
      integer(kind=4),dimension(:,:),allocatable :: MSet
      ! Use a one dimensional array to hold data for the global image
      integer(kind=4),dimension(:),allocatable :: myMSet
      ! Change maximum and minimum values of the window that the
      ! appear in the image. Increase the values and you zoom out
      ! decrease the values and you zoom in
      integer(kind=4) :: xmin = -2
      integer(kind=4) :: ymin = -2
      integer(kind=4) :: xmax = 2
      integer(kind=4) :: ymax = 2
      real(kind=8) :: threshold = 1
      real(kind=8) :: dist = 0
      integer(kind=4) :: ix, iy, ii, proc
      real(kind=8) :: cx, cy
      integer(kind=4) :: iter, i, ierror, allocatestatus
      real(kind=8) :: x,y,x2,y2
      real(kind=8) :: temp = 0.0
      real(kind=8) :: xder = 0.0
      real(kind=8) :: yder = 0.0
      real(kind=8),dimension(:),allocatable :: xorbit
      real(kind=8),dimension(:),allocatable :: yorbit
      integer(kind=4) :: hugenum = 100000
      logical :: flag = .false.
      integer(kind=4),parameter :: overflow = 2147483647
      real(kind=8) :: delta 
      ! chunk size for dynamic scheduling of shared memory loop
      integer(kind=4),parameter :: chunk = 5
      ! name of file to write image out to
      character ( len = 80 ) :: file_name = 'mandelbrot.ascii.pgm'
      delta = (threshold*(xmax-xmin))/real(nx-1)
      ! MPI variables
      integer(kind=4) :: ierr, myid, numprocs, mystart,myend,itemcount,procstart,procend
      integer(kind=4) stat(MPI_STATUS_SIZE) 
      ! initialisation of MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      if(myid.eq.0) then
        mystart=1+floor(myid*Nx*Ny/real(numprocs))
        myend=min(floor((myid+1.0d0)*Nx*Ny/real(numprocs)),Nx*Ny)
        allocate(MSet(1:nx,1:ny),myMSet(mystart:myend), xorbit(maxiter), yorbit(maxiter),stat=allocatestatus)
        if (allocatestatus .ne. 0) stop	
      else
        mystart=1+floor(myid*Nx*Ny/real(numprocs))
        myend=min(floor((myid+1.0d0)*Nx*Ny/real(numprocs)),Nx*Ny)
        allocate(myMSet(mystart:myend), xorbit(maxiter), yorbit(maxiter),stat=allocatestatus)
        if (allocatestatus .ne. 0) stop	
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      ! We use a nested loop here to effectively traverse over each part of the grid (pixel of the image)
      ! in sequence. First, the complex values of the points are determined and then used as the basis of 
      ! the computaion. Effectively, it will loop over each point (pixel) and according on how many iterations 
      ! it takes for the value that the mathematical function returns on each iteration it will determine 
      ! whether or not the point "escapes" to infinity (or an arbitrarily large number.) or not. If it takes 
      ! few iterations to escape then it will decide that this point is NOT part of the Mandelbrot set and 
      ! will put a 0 in that point's index in MSet. If it takes nearly all or all of the iterations to escape, 
      ! then it will decide that the point/pixel is part of the Mandelbrot set and instead put a 1 in its place 
      ! in myMSet.

      ! The use of the OpenMP pragma here will divide up the iterations between threads and execute them in parallel
      ! This region is VERY easily parallelized because there is NO data shared between the loop iterations.
      !$OMP PARALLEL DO PRIVATE(ii,ix,iy,cx,cy,iter,i,x,y,x2,y2,temp,xder,yder,dist,xorbit,yorbit,flag) &
      !$OMP SHARED(myMSet) SCHEDULE(DYNAMIC,Chunk)
      do ii=mystart,myend
        iy=1+floor(ii/real(Nx))
        ix=1+ii-(iy-1)*Nx
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
        ! This is the main loop that determins whether or not the point escapes or not. 
        ! It breaks out of the loop when it escapes
        do iter = 0,maxiter
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
        ! if the point escapes, find the distance from the set, just in case its close 
        ! to the set. if it is, it will make it part of the set.

        if (x2+y2 .ge. hugenum) then
          xder = 0
          yder = 0
          i = 0
          flag = .false.

          do i=0, iter
                
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
        ! Assign the appropriate values to MSet in the place relating to the point in question
        if (dist .lt. delta) then
          ! MSet(ix,iy) = 1
          myMSet(ix+Nx*iy) = 1
        else
          ! MSet(ix,iy) = 0
          myMSet(ix+Nx*iy) = 0
        endif
      enddo  
      !$OMP END PARALLEL DO

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ! Get Process 0 to collect data and write out a file. This is scalable to a moderate 
      ! number of processes, limited by memory on process 0.
      if (myid.eq.0) then	
        do ii=mystart,myend
          iy=1+floor(ii/real(Nx))
          ix=1+ii-(iy-1)*Nx
          MSet(ix,iy)=myMSet(ix+Nx*iy)
        end do

        if (numprocs.gt.1) then
          do proc=1,numprocs
            procstart=1+floor(proc*Nx*Ny/real(numprocs))
            procend=min(1+floor((proc+1.0d0)*Nx*Ny/real(numprocs)),Nx*Ny)
            itemcount=1+procend-procstart
            call MPI_RECEIVE(myMSet,itemcount,MPI_INTEGER,proc,proc,MPI_COMM_WORLD,stat,ierr)
            ! copy data into main array
            do ii=procstart,procend
              iy=1+floor(ii/real(Nx))
              ix=1+ii-(iy-1)*Nx
              MSet(ix,iy)=myMSet(mystart+ix+Nx*iy-procstart)
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
        ! Call Burkardt's Fortran code to write Mandelbrot picture to disk
        call pgma_write ( file_name, nx, ny, Mset, ierror )
        deallocate(MSet,myMSet,xorbit,yorbit,stat=allocatestatus)
        if (allocatestatus .ne. 0) stop	
      else
        deallocate(myMSet,xorbit,yorbit,stat=allocatestatus)
        if (allocatestatus .ne. 0) stop	
      end if
      call MPI_FINALIZE(ierr)	
end program Mandelbrot


subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
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
subroutine getint ( done, ierror, inunit, ival, string )

!*****************************************************************************80
!
!! GETINT reads an integer from a file.
!
!  Discussion:
!
!    The file, or at least the part read by GETINT, is assumed to
!    contain nothing but integers.  These integers may be separated
!    by spaces, or appear on separate lines.  Comments, which begin
!    with "#" and extend to the end of the line, may appear anywhere.
!
!    Each time GETINT is called, it tries to read the next integer
!    it can find.  It remembers where it was in the current line
!    of text.
!
!    The user should open a text file on FORTRAN unit INUNIT,
!    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
!    care of reading in a new STRING as necessary, and extracting
!    as many integers as possible from the line of text before 
!    reading in the next line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical DONE.
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!    On output, if there is no more data to be read from STRING,
!    then DONE is TRUE.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to read the integer.
!
!    Input, integer ( kind = 4 ) INUNIT, the FORTRAN unit from which to read.
!
!    Output, integer ( kind = 4 ) IVAL, the integer that was read.
!
!    Input/output, character ( len = * ) STRING, the text of the most recently 
!    read line of the file.
!
  implicit none

  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) string
  character ( len = 80 ) word

  do

    call word_next_rd ( string, word, done )

    if ( .not. done ) then
      exit
    end if

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    i = index ( string, '#' )
    if ( i /= 0 ) then
      string(i:) = ' '
    end if

  end do

  call s_to_i4 ( word, ival, ierror, last )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINT - Fatal error!'
    write ( *, '(a)' ) '  Error trying to convert string to integer.'
    stop
  end if

  return
end
subroutine pgma_check_data ( row_num, col_num, g_max, g, ierror )

!*****************************************************************************80
!
!! PGMA_CHECK_DATA checks ASCII PGM data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROW_NUM,  COL_NUM, the number of rows 
!    and columns of data.
!
!    Input, integer ( kind = 4 ) G_MAX, the maximum gray value.
!
!    Input, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), contains the gray data.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) g_max

  ierror = 0

  if ( minval ( g(1:row_num,1:col_num) ) < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
    write ( *, '(a)' ) '  At least one gray value is below 0.'
    ierror = 1
    stop
  end if

  if ( g_max < maxval ( g(1:row_num,1:col_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
    write ( *, '(a)' ) '  At least one gray value exceeds G_MAX.'
    write ( *, '(a,i12)' ) '  G_MAX = ', g_max
    ierror = 1
    stop
  end if

  return
end
subroutine pgma_example ( row_num, col_num, g )

!*****************************************************************************80
!
!! PGMA_EXAMPLE sets up sample ASCII PGM data.
!
!  Discussion:
!
!    The data is based on three periods of a sine curve.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.  A reasonable value is 200 for ROW_NUM and 600 
!    for COL_NUM.
!
!    Output, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), the gray scale data.
!
  implicit none

  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) col_num

  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: periods = 3
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do i = 1, row_num
    y = 2.0D+00 * real ( i - 1, kind = 8 ) &
      / real ( row_num - 1, kind = 8 ) - 1.0D+00
    do j = 1, col_num
      x = 2.0D+00 * pi * real ( periods * ( j - 1 ), kind = 8 ) &
        / real ( col_num - 1, kind = 8 )
      g(i,j) = int ( 20.0D+00 * ( sin ( x ) - y + 2.0D+00 ) )
    end do
  end do

  return
end
subroutine pgma_read_data ( file_in_unit, row_num, col_num, g )

!*****************************************************************************80
!
!! PGMA_READ_DATA reads the data in an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), the gray data.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 80 ) string

  ierror = 0
  done = .true.
  string = ' '

  do i = 1, row_num
    do j = 1, col_num

      call getint ( done, ierror, file_in_unit, g(i,j), string )

      if ( ierror /= 0 ) then
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G data.'
        stop
      end if

    end do
  end do

  return
end
subroutine pgma_read_header ( file_in_unit, row_num, col_num, g_max )

!*****************************************************************************80
!
!! PGMA_READ_HEADER reads the header of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer ( kind = 4 ) G_MAX, the maximum gray value.
!
  implicit none

  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 2 )  magic
  integer ( kind = 4 ) g_max
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  logical s_eqi
  character ( len = 80 ) string
!
!  Read the first line of data, which must begin with the magic number.
!
  read ( file_in_unit, '(a)', iostat = ios ) magic

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  End or error while reading file.'
    ierror = 2
    stop
  end if

  if ( .not. s_eqi ( magic, 'P2' ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error.'
    write ( *, '(a)' ) '  First two bytes are not magic number "P2".'
    write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
    stop
  end if
!
!  Now search for COL_NUM, ROW_NUM, and G_MAX.
!
  done = .true.
  string = ' '

  call getint ( done, ierror, file_in_unit, col_num, string )

  if ( ierror /= 0 ) then
    close ( unit = file_in_unit )
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading COL_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, row_num, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading ROW_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, g_max, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading G_MAX.'
    stop
  end if

  return
end
subroutine pgma_read_test ( file_in_name, ierror )

!*****************************************************************************80
!
!! PGMA_READ_TEST tests an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file 
!    containing the ASCII PGM data.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag which is nonzero if
!    there was an error.
!
  implicit none

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) g_max
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Read the header.
!
  call pgma_read_header ( file_in_unit, row_num, col_num, g_max )
!
!  Allocate the data.
!
  allocate ( g(row_num,col_num) )
!
!  Read the data.
!
  call pgma_read_data ( file_in_unit, row_num, col_num, g )

  close ( unit = file_in_unit )
!
!  Check the data.
!
  call pgma_check_data ( row_num, col_num, g_max, g, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST - Warning!'
    write ( *, '(a)' ) '  PGMA_CHECK_DATA did not approve the data.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_TEST:'
    write ( *, '(a)' ) '  PGMA_CHECK_DATA has approved the data from the file.'
  end if

  deallocate ( g )

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
!    01 March 2003
!
!  Author:
!
!    John Burkardt
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

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  logical, parameter :: debug = .false.
  character ( len = * )  file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
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
    write ( *, '(a)' ) '  The data was checked and written.'
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
!    28 February 2003
!
!  Author:
!
!    John Burkardt
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

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
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
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
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

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  character ( len = 2 ) :: magic = 'P2'
  integer ( kind = 4 ) g_max
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

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
subroutine pgma_write_test ( file_out_name )

!*****************************************************************************80
!
!! PGMA_WRITE_TEST tests the ASCII PGM write routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to contain the ASCII PGM data.
!
  implicit none

  character ( len = * ) file_out_name
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  row_num = 300
  col_num = 300
!
!  Allocate memory.
!
  allocate ( g(row_num,col_num) )
!
!  Set the data.
!
  call pgma_example ( row_num, col_num, g )
!
!  Write the data to the file.
!
  call pgma_write ( file_out_name, row_num, col_num, g, ierror )

  deallocate ( g );

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_WRITE_TEST - Fatal error!'
    write ( *, '(a)' ) '  PGMA_WRITE failed.'
    stop
  end if

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lens
  character ( len = * ) s

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( s )

  i = 0

  do

    i = i + 1

    c = s(i:i)

    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      exit
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      exit
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine word_next_rd ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_RD "reads" words from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    words separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!    If DONE is TRUE, then you need to provide a new LINE of data.
!
!  Local Parameters:
!
!    NEXT is the next location in LINE that should be searched.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lenl
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  character ( len = 1 ), parameter :: TAB = char(9)
  character ( len = * ) word

  lenl = len_trim ( line )

  if ( done ) then
    next = 1
    done = .false.
  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank.
!
  ilo = next

  do
!
!  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
!
    if ( lenl < ilo ) then
      word = ' '
      done = .true.
      next = lenl + 1
      return
    end if
!
!  ...If the current character is blank, skip to the next one.
!
    if ( line(ilo:ilo) /= ' ' .and. line(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  To get here, ILO must be the index of the nonblank starting
!  character of the next word.
!
!  Now search for the LAST nonblank character.
!
  next = ilo + 1

  do

    if ( lenl < next ) then
      word = line(ilo:next-1)
      return
    end if

    if ( line(next:next) == ' ' .or. line(next:next) == TAB ) then
      exit
    end if

    next = next + 1

  end do

  word = line(ilo:next-1)

  return
end
