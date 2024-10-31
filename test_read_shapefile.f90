program test_read_shapefile
use miscellaneous_functions, only: remove_extension
use shapefile_io_module, only: inquire_shapefile, read_shapefile
implicit none

    character(len=500):: fname
    integer:: shptype, nrec, npnttot, npntmax, nprttot, nprtmax
    double precision, dimension(4):: bbox
    integer, dimension(:), allocatable:: headlen, npoints, nparts, recpos, recpos_in_part, partpos, partlen
    double precision, dimension(:,:), allocatable:: points
    integer:: i, j, k, n


    ! Input shapefile (get command-line argument, or interactively ask)
    n = command_argument_count()
    if (n >= 1) then
        call get_command_argument(1, fname)
    else
        write(unit=*, fmt='(A)', advance='no') 'Enter input shapefile: '
        read(unit=*, fmt='(A)') fname
    end if
    call remove_extension(fname, list_ext=(/'.shp', '.shx', '.dbf', '.prj'/))
    print *, fname

    ! Initialization
    shptype = -1
    nrec = -1
    npnttot = -1
    nprttot = -1
    npntmax = -1
    nprtmax = -1
    bbox = (/-1d36, -1d36, -1d36, -1d36/)


    ! Inquire main attribute of shapefile
    call inquire_shapefile(trim(fname), nrec, headlen, npoints, nparts, &
                           shptype=shptype, npointstot=npnttot, npointsmax=npntmax, npartstot=nprttot, npartsmax=nprtmax, bbox=bbox)

    print *
    print *, 'ShapeType:        ', shptype
    print *, 'BBOX:             ', bbox
    print *, '# of records:     ', nrec
    print *, '# of parts (tot): ', nprttot
    print *, 'max # of parts:   ', nprtmax
    print *, '# of points (tot):', npnttot
    print *, 'max # of points:  ', npntmax
    print *
    print *, 'header length:', headlen
    print *
    print *, 'rec_nparts:           ', nparts
    print *, 'rec_npoints:          ', npoints
    print *

    ! Read shapefile data
    call read_shapefile(trim(fname), nrec, npoints, nparts, recpos, recpos_in_part, partpos, partlen, points, bbox)

    print *
    print *, '# of records:     ', nrec
    print *, '# of points (tot):', size(points,2)
    print *, '# of parts (tot): ', size(partpos)
    print *
    print *, 'records npoints:            ', npoints
    print *, 'records nparts:             ', nparts
    print *, 'records pos:                ', recpos
    print *, 'records pos in parts arrays:', recpos_in_part
    print *, 'parts pos:                  ', partpos
    print *, 'parts length:               ', partlen
    print *

    ! Write shapefile data in text file
    open(unit=1, file='toto.txt', status='replace', action='write')
    n = 0
    do i = 1,nrec
      do j = 1,nparts(i)
        do k = 1,partlen(recpos_in_part(i)+j-1)
          n = n + 1
          write(unit=1,fmt=*) points(1,n), points(2,n)
        end do
        write(unit=1,fmt=*) 'nan nan'
      end do
    end do
    print *, 'points => written into file "toto.txt"'
    print *, n, size(points,2)
    close(unit=1)

end program
