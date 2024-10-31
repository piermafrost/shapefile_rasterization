program shp2grid

use scan_arguments_mod, only: scan_arguments
use miscellaneous_functions, only: remove_extension
use shapefile_io_module, only: inquire_shapefile, read_shapefile
use netcdf_output, only: create_output_file, write_output, write_output_slice, close_file
use geography, only: get_earth_parameters, authalic_latitude, authalic_radius
use shp2grid_functions, only: iNA, NA, define_grid, get_grid_from_netcdf_file, get_pos, trapez_area, close_cell_ring, &
                              record_border_in, record_border_out, reinitialize_border, get_connected_points, fill_polygon_interior
implicit none


include 'path.inc' ! <-- get variable "root_path"

!############################# PARAMETERS: #############################!

double precision, parameter:: PI = acos(-1d0)
double precision, parameter:: DEG2RAD = PI/180d0
double precision, parameter:: DEFAULT_TOLERANCE = 10*epsilon(1d0)
character(len=*), parameter:: DEFAULT_EARTH_MODEL = 'WGS84'

! maximal length of character variables (optional arguments, file names, variable names)
integer, parameter:: ARGLEN=1000, FLEN=500, VLEN=80

!#######################################################################!

double precision, dimension(:,:),   allocatable:: grid_area, polygon_area, polygon_area_sum
double precision, dimension(:,:),   allocatable:: polygon, pol_part
double precision, dimension(:,:,:), allocatable:: last_border_pos
double precision, dimension(:),     allocatable:: lon, lat, lonbnd, latbnd, xbnd, ybnd
double precision, dimension(:,:),   allocatable:: xyadd
double precision, dimension(:),     allocatable:: xadd, yadd, xaddtiming, yaddtiming
logical, dimension(:,:,:),          allocatable:: is_last_border_out
logical, dimension(:,:),            allocatable:: visited, connected
integer, dimension(:,:),            allocatable:: list_visited, list_connected
integer, dimension(:,:),            allocatable:: ijadd
integer, dimension(:),              allocatable:: recpos, recpos_in_part, npoints, nparts, partpos, partlen
integer:: nvisited, nconnected

double precision:: y0, x1, y1, x2, y2, x2_add, y2_add, xin, yin, area, recarea, rastarea, totarea, totrastarea
double precision, dimension(4):: bbox, grid_bbox
integer:: i, j, k0, k, n, ishift, iadd, jadd, kadd, naddi, naddj, nadd, di, dj, dj0, irec, ipart, npts, nrec, listsize
integer:: ofid, ovarid
logical:: flip_x, flip_y, goon
character(len=ARGLEN):: dummychar

! Arguments and options
integer:: nx, ny
character(len=FLEN):: input_file, output_file, grid_file
character(len=VLEN):: xvarname, xbndvarname, yvarname, ybndvarname
character(len=VLEN):: earth_model
double precision:: dlon, dlat, radius, authradius, flat, ecc, tolerance
logical:: check_tot_area, check_rec_area, area_fract, global_grid, use_grid_file, use_shpf_bbox, split_records, split_computation, &
          split_connected_computation, check, split, print_progress

! Variables for subdivision in N task, print on screen when each task in completed
integer, parameter:: NTASK = 132
integer, dimension(NTASK+1):: list_task
integer:: itask
character(len=6):: barformat
character(len=NTASK) barprogress

! Variables concerning user options
integer, parameter:: NOPT = 24, NREQARG_MAX = 1
character(len=VLEN), dimension(NOPT):: options
integer, dimension(NOPT)::           nreqarg
character(len=ARGLEN), dimension(NOPT,NREQARG_MAX):: opt_arg
logical, dimension(NOPT):: got_opt
! Variables concerning mandatory arguments
integer, parameter:: NARG = 1
character(len=200), dimension(NARG):: list_arg




!##############################################################################################################################!
!=========================        list of programs options the user may pass to the executable:        ========================!
!
options(1)  = '-o' ! -o output_file_name
nreqarg(1)  = 1
!
options(2)  = '--size' ! --size nxy  OR  -s nx,ny
options(3)  = '-s'
nreqarg(2)  = 1
nreqarg(3)  = 1
!
options(4)  = '--res' ! --res dxy  OR  --res dx,dy
options(5)  = '-r'
nreqarg(4)  = 1
nreqarg(5)  = 1
!
options(6)  = '--match' ! --match file_name,xvar,yvar  OR  --match file_name,xvar,yvar,xbndvar,ybndvar
options(7)  = '-m'
nreqarg(6)  = 1
nreqarg(7)  = 1
!
options(8) = '--bbox' ! --bbox x0,y0,x1,y1
nreqarg(8) = 1
!
options(9)  = '--use-shp-bbox'
nreqarg(9)  = 0
!
options(10) = '--global-grid'
options(11) = '-g'
nreqarg(10) = 0
nreqarg(11) = 0
!
options(12) = '--fraction'
options(13) = '-f'
nreqarg(12) = 0
nreqarg(13) = 0
!
options(14) = '--earth-model' ! --earth-model model_name
nreqarg(14) = 1
!
options(15) = '--earth-flattening' ! --earth-flattening x
nreqarg(15) = 1
!
options(16) = '--earth-radius' ! --earth-radius x
nreqarg(16) = 1
!
options(17) = '--check-tot-area'
nreqarg(17) = 0
!
options(18) = '--check-rec-area'
nreqarg(18) = 0
!
options(19) = '--check'
nreqarg(19) = 0
!
options(20) = '--protect-adjacent-polygons'
nreqarg(20) = 0
!
options(21) = '--split-records'
nreqarg(21) = 0
!
options(22) = '--split-computation'
nreqarg(22) = 0
!
options(23) = '--tolerance' ! --tolerance x
nreqarg(23) = 1
!
options(24) = '--help'
nreqarg(24) = 0
!
!==============================================================================================================================!
!##############################################################################################################################!




! Display of progress bar:
write(barformat, fmt='(A2,I3.3,A1)') '(A',NTASK,')'
barprogress = ''
barprogress(1:4) = '| 0%'
barprogress(NTASK-5:NTASK) = '100% |'




!##############################################################################################################################!
!==============================================================================================================================!
!=========================               Get options passed to the program by the user                =========================!
!==============================================================================================================================!
!##############################################################################################################################!


print *

! Call subroutine retrieving regular arguments and options !
! ======================================================== !

call scan_arguments(list_arg, list_legal_opt=options, list_expected_n_optarg=nreqarg, list_optarg=opt_arg, got_opt=got_opt, &
                    inquire_missing_arg=.false.)


! Help message !
! ============ !

if (got_opt(24)) then
    open(unit=13, status='old', action='read', file=root_path//'/help_message.txt')
    i = 0
    dummychar = ''
    do while (i==0)
        print *, trim(dummychar)
        read(unit=13, fmt='(A)', iostat=i) dummychar
    end do
    close(unit=13)
    stop 0
end if


! Assign retrieved arguments to code variables !
! ============================================ !

! Input file
! ----------
input_file = list_arg(1)
if (input_file == '') then
    print *, 'ERROR: missing required argument "input shapefile"'
        print *, 'Program stopped'
    stop 1
end if
call remove_extension(input_file, list_ext=(/'.shp', '.shx', '.dbf', '.prj'/))

! Output file
! -----------
if (got_opt(1)) then
    output_file = opt_arg(1,1)(1:FLEN)
else
    output_file = trim(input_file)//'_rasterized.nc'
end if

! initialize with "undefined" values
use_grid_file = .false.
use_shpf_bbox = .false.
grid_file   = ''
xvarname    = ''
xbndvarname = ''
yvarname    = ''
ybndvarname = ''
nx = iNA
ny = iNA
dlon = NA
dlat = NA
grid_bbox = NA

! Match grid from another netCDF File
! -----------------------------------
if (got_opt(6) .or. got_opt(7)) then
    if (got_opt(6)) then
        i = 6
    else
        i = 7
    end if
    dummychar = opt_arg(i,1)
    ! count number of commas
    n = 0
    do k = 1,len_trim(dummychar)
        if (dummychar(k:k)==',') n = n+1
    end do
    if (n/=2 .and. n/=4) then
        print *, 'ERROR: unqualified statement:'
        print *, '        "'//trim(options(i))//' '//trim(dummychar)//'"'
        print *, 'Expect: "'//trim(options(i))//' file_name,x_var_name,y_var_name"'
        print *, '    or: "'//trim(options(i))//' file_name,x_var_name,y_var_name,xbnds_var_name,ybnds_var_name"'
        print *, 'Program stopped'
        stop 1
    end if
    k = 1
    k0 = k
    do while (dummychar(k:k)/=',')
        k = k+1
    end do
    grid_file = dummychar(k0:k-1)
    k = k+1
    k0 = k
    do while (dummychar(k:k)/=',')
        k = k+1
    end do
    xvarname = dummychar(k0:k-1)
    k = k+1
    k0 = k
    if (n==2) then
        yvarname = trim(dummychar(k0:))
    else ! --> n==4
        do while (dummychar(k:k)/=',')
            k = k+1
        end do
        yvarname = dummychar(k0:k-1)
        k = k+1
        k0 = k
        do while (dummychar(k:k)/=',')
            k = k+1
        end do
        xbndvarname = dummychar(k0:k-1)
        k = k+1
        k0 = k
        ybndvarname = trim(dummychar(k0:))
    end if
    use_grid_file = .true.

! Grid size
! ---------
elseif (got_opt(2) .or. got_opt(3)) then
    if (got_opt(2)) then
        i = 2
    else
        i = 3
    end if
    dummychar = opt_arg(i,1)
    ! count number of commas
    n = 0
    do k = 1,len_trim(dummychar)
        if (dummychar(k:k)==',') n = n+1
    end do
    select case (n)
        case (0)
            read(opt_arg(i,1), fmt=*) nx
            ny = nx
        case (1)
            read(opt_arg(i,1), fmt=*) nx,ny
        case default
            print *, 'ERROR: unqualified statement:'
            print *, '        "'//trim(options(i))//' '//trim(dummychar)//'"'
            print *, 'Expect: "'//trim(options(i))//' nxy" or "'//trim(options(i))//' nx,ny"'
            print *, 'Program stopped'
            stop 1
    end select

! Grid Resolution
! ---------------
elseif (got_opt(4) .or. got_opt(5)) then
    if (got_opt(4)) then
        i = 4
    else
        i = 5
    end if
    dummychar = opt_arg(i,1)
    ! count number of commas
    n = 0
    do k = 1,len_trim(dummychar)
        if (dummychar(k:k)==',') n = n+1
    end do
    select case (n)
        case (0)
            read(opt_arg(i,1), fmt=*) dlon
            dlat = dlon
        case (1)
            read(opt_arg(i,1), fmt=*) dlon,dlat
        case default
            print *, 'ERROR: unqualified statement:'
            print *, '        "'//trim(options(i))//' '//trim(dummychar)//'"'
            print *, 'Expect: "'//trim(options(i))//' dxy" or "'//trim(options(i))//' dx,dy"'
            print *, 'Program stopped'
            stop 1
    end select

else
    print *, 'ERROR: not enough information given about the grid.'
    print *, 'Use either "--res" (grid resolution), "--size" (grid size) or "--match" (to match'
    print *, 'the grid of a given file).'
    print *, 'Type `shp2grid --help` for more information.'
    print *, 'Program stopped'
    stop 1
end if

! Global (T) or local (F) grid
! ----------------------------
global_grid = (got_opt(10) .or. got_opt(11))

! Bounding-box
! ------------
if (got_opt(8)) then
    i = 8
    dummychar = opt_arg(i,1)
    ! count number of commas
    n = 0
    do k = 1,len_trim(dummychar)
        if (dummychar(k:k)==',') n = n+1
    end do
    if (n/=3) then
        print *, 'ERROR: unqualified statement:'
        print *, '        "'//trim(options(i))//' '//trim(dummychar)//'"'
        print *, 'Expect: "'//trim(options(i))//' x0,y0,x1,y1"'
        print *, 'Program stopped'
        stop 1
    end if
    k = 1
    k0 = k
    do n = 1,3
        do while (dummychar(k:k)/=',')
            k = k+1
        end do
        if (k0 /= k) read(dummychar(k0:k-1), fmt=*) grid_bbox(n)
        k = k+1
        k0 = k
    end do
    if (trim(dummychar(k0:)) /= '') read(dummychar(k0:), fmt=*) grid_bbox(4)

! Whether or not using input shapefile bounding-box
! -------------------------------------------------
elseif (got_opt(9)) then
    if (global_grid) then
        print *, 'Warning: option "--use-shp-bbox" is incompatible with option "--global-grid".'
        print *, 'Will be ignored' 
        print *
    else
        use_shpf_bbox = .true.
    end if
end if

if ((.not. global_grid) .and. (.not. use_grid_file) .and. (.not. use_shpf_bbox) .and. any(grid_bbox==NA)) then
    print *, 'ERROR: not enough information about the bounding box for the rasterization.'
    print *, 'Use either:'
    print *, '    "--bbox x0,y0,x1,y1" WITHOUT SKIPPED VALUES, to specify the bounding box'
    print *, '    "--use-shp-bbox" to keep the bounding box from the input shapefile'
    print *, '    "--match ..." to fully define the grid (bounding box is then optional)'
    print *, '    "--global-grid" to indicate that the grid must cover the entire Earth'
    print *, 'Type `shp2grid --help` for more information.'
    print *, 'Program stopped'
    stop 1
end if

! Write polygon area in grid cell or polygon fraction of grid cell
! ----------------------------------------------------------------
area_fract = (got_opt(12) .or. got_opt(13))

! Earth model
! -----------
if (got_opt(14)) then
    read(opt_arg(14,1), fmt=*) earth_model
else

    if (got_opt(15) .and. got_opt(16)) then

! Earth flattening
! ----------------
        read(opt_arg(15,1), fmt=*) flat

! Earth radius
! ------------
        read(opt_arg(16,1), fmt=*) radius

        earth_model = 'custom'

    else
        earth_model = DEFAULT_EARTH_MODEL
    end if
end if

! Whether or not checking the total area (sum of all records)
! -----------------------------------------------------------
check_tot_area = (got_opt(17) .or. got_opt(19))

! Whether or not checking Check the area of each record
! -----------------------------------------------------
check_rec_area = (got_opt(18) .or.got_opt(19))

! Whether or not splitting computation of connected points avoid errors with adjactent polygons
! ---------------------------------------------------------------------------------------------
split_connected_computation = got_opt(20)

! Whether or not writing each record separately
! ---------------------------------------------
split_records = got_opt(21)

! Whether or not processing each record independently, and sum all the areas at the end
! -------------------------------------------------------------------------------------
split_computation = got_opt(22)

! Tolerance (minimum distance for a point to be considered as not in a cell corner)
! ---------------------------------------------------------------------------------
if (got_opt(23)) then
    read(opt_arg(23,1), fmt=*) tolerance
else
    tolerance = DEFAULT_TOLERANCE
end if


! Actual code options:
check          = check_tot_area .or. check_rec_area
split          = split_records .or. split_computation
print_progress = (.not. split) .or. (.not. check_rec_area)

! consider global grid if user-defined bounding-box is consistent
if ((.not. global_grid) .and. all(grid_bbox/=NA)) then
    if (grid_bbox(2)==-90d0 .and. grid_bbox(4)==90d0 .and. grid_bbox(3)==grid_bbox(1)+360d0)  global_grid = .true.
end if




!##############################################################################################################################!
!==============================================================================================================================!
!=========================                   load input file / create output file                     =========================!
!==============================================================================================================================!
!##############################################################################################################################!


if (earth_model/='custom') then
    call get_earth_parameters(earth_model, radius=radius, eccentricity=ecc)
end if
authradius = authalic_radius(radius, ecc)


! Get shapefile bounding-box information
call inquire_shapefile(input_file, bbox=bbox)


! safety check
if (global_grid .and. bbox(3)>bbox(1)+360d0) then
    print *, 'ERROR:'
    print *, 'with "--global-grid" option, cannot handle shapefile data spending over more'
    print *, 'than 360 degrees of longitude'
    print *, 'Program stopped'
    stop 2
elseif (bbox(2)<-90d0 .or. bbox(4)>90d0) then
    print *, 'ERROR:'
    print *, 'cannot handle shapefile data exceeding 90°N or 90°S'
    print *, 'Program stopped'
    stop 2
end if


if (use_grid_file) then

    ! grid information are taken from an inpu netCDF file
    ! => fully define the grid and initialize polygon area before checking the shapefile bounding box 
    ! ***********************************************************************************************
    call get_grid_from_netcdf_file(grid_file, xvarname, yvarname, xbndvarname, ybndvarname, authradius, ecc, global_grid, bbox, &
                                   grid_bbox, nx, ny, lon, lat, lonbnd, latbnd, xbnd, ybnd, flip_x, flip_y, grid_area, polygon_area)

else

    if (global_grid) then

        ! longitude bounds
        if (grid_bbox(1)==NA .and. grid_bbox(3)==NA) then
            grid_bbox(1) = -180d0
            grid_bbox(3) = 180d0
        elseif (grid_bbox(3)==NA) then
            grid_bbox(3) = grid_bbox(1) + 360d0
        elseif (grid_bbox(1)==NA) then
            grid_bbox(1) = grid_bbox(3) - 360d0
        else
            grid_bbox(1) = (grid_bbox(1) + grid_bbox(3)-360d0)/2
            grid_bbox(3) = grid_bbox(1) + 360d0
        end if

        if (use_grid_file .or. (nx/=iNA .and. ny/=iNA)) then
            ! latituide bounds
            grid_bbox(2) = -90d0
            grid_bbox(4) = 90d0

        else!--> "dlon" and "dlat" were defined
            ! redefine longitude bounds:
            grid_bbox(1) = floor(grid_bbox(1)/dlon)*dlon
            grid_bbox(3) = grid_bbox(1) + ceiling(360d0/dlon)*dlon
            ! latitude bounds:
            grid_bbox(2) = floor(-90d0/dlat)*dlat
            grid_bbox(4) = ceiling(90d0/dlat)*dlat
            ! NOTE: lat_bnds array will be corrected to ensure that it doesn't exceed [-90,90]
        end if

    elseif (use_shpf_bbox) then

        if (use_grid_file .or. (nx/=iNA .and. ny/=iNA)) then
            grid_bbox = bbox
        else ! => "dlon" and "dlat" were defined
            grid_bbox = (/ floor(bbox(1)/dlon)*dlon, floor(bbox(2)/dlat)*dlat, &
                        ceiling(bbox(3)/dlon)*dlon, ceiling(bbox(4)/dlat)*dlat /)
        end if

    elseif ( (grid_bbox(1)/=NA .and. bbox(1)<grid_bbox(1)) .or. &
             (grid_bbox(2)/=NA .and. bbox(2)<grid_bbox(2)) .or. &
             (grid_bbox(3)/=NA .and. bbox(3)>grid_bbox(3)) .or. &
             (grid_bbox(4)/=NA .and. bbox(4)>grid_bbox(4)) ) then
        print *, 'ERROR: Shapefile bounding-box exceeds current grid'
        print *, 'Shapefile bbox:', bbox
        print *, 'Grid bbox:     ', grid_bbox
        print *, 'Program stopped'
        stop 2

    end if

    ! Grid definition and initialization of polygon area
    ! **************************************************
    call define_grid(grid_bbox, global_grid, nx, ny, dlon, dlat, authradius, ecc, lon, lat, lonbnd, latbnd, xbnd, ybnd, grid_area, &
                     polygon_area)
    flip_x = .false.
    flip_y = .false.

end if


! %%%%%%%%%%%%%%%%%%%% !
! Load input shapefile !
! %%%%%%%%%%%%%%%%%%%% !

call read_shapefile(input_file, nrec, npoints, nparts, recpos, recpos_in_part, partpos, partlen, polygon)


print *
print *, 'Pre-computation...'

! degree -> radian conversion:
polygon = DEG2RAD*polygon

! conversion in authalic latitude:
call authalic_latitude(polygon(2,:) , ecc)


! %%%%%%%%%% !
! Allocation !
! %%%%%%%%%% !

listsize = min( (nx+1)*ny , int(10*sqrt(real((nx+1)*ny))) )
allocate( pol_part(2,maxval(partlen)) )
allocate(       list_visited(2,(nx+1)*ny) )
allocate(     list_connected(2,(nx+1)*ny) )
allocate(            visited((nx+1),ny)   ) 
allocate(          connected((nx+1),ny)   ) 
allocate(    last_border_pos(4,(nx+1),ny) )
allocate( is_last_border_out(4,(nx+1),ny) )
allocate(      ijadd(2,2*listsize) )
allocate(      xyadd(2,2*listsize) )
allocate(       xadd(listsize)   )
allocate( xaddtiming(listsize)   )
allocate(       yadd(listsize)   )
allocate( yaddtiming(listsize)   )
if (split_computation) then
    allocate( polygon_area_sum(nx,ny) )
    polygon_area_sum = 0d0
end if



! Create output file
! ------------------

if (split_records) then
    call create_output_file(output_file, lon, lat, lonbnd, latbnd, flip_x, flip_y, ofid, ovarid, nrec=nrec, &
                            earth_model=earth_model, grid_area=grid_area(1:nx,:), area_fraction=area_fract, ifname=input_file)
else
    call create_output_file(output_file, lon, lat, lonbnd, latbnd, flip_x, flip_y, ofid, ovarid, &
                            earth_model=earth_model, grid_area=grid_area(1:nx,:), area_fraction=area_fract, ifname=input_file)
end if



! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
! If the grid is global: check longitudinal span of shapefile, compare to the defined grid. !
! If needed, shift the grid to ensure it covers the shapefile longitudinal span             !
! + + + + + + + + + + + + + + + + +  + + + + + ++ + + + + + + + + + + + + + + + + + + + + + !

ishift = 1

if (global_grid) then

    k = int(ceiling((bbox(1) - lonbnd(nx+1)) / 360d0))

    ! modulo-360 grid shifting
    ! modify "bbox" instead of "lonbnd" for the current computation,
    ! because the original "lonbnd" will written in output file, while "bbox" will not be used after
    bbox(1) = bbox(1) - k*360d0
    bbox(3) = bbox(3) - k*360d0
    xbnd = xbnd + k*2*PI

    do while (lonbnd(ishift+1)+360d0 < bbox(3))
        ishift = ishift + 1
    end do

    ! split the grid at ishift, and generate a shifted grid that covers the shapefile longitudinal span
    xbnd = (/ xbnd(ishift:nx), xbnd(1:ishift+1)+2*PI /)
    do j = 1,ny
        grid_area(:,j) = (/ grid_area(ishift:nx,j), grid_area(1:ishift,j) /)
    end do

end if

! inverse the shifting point, for reverse-shifting to write outputs
ishift = nx + 2 - ishift




!##############################################################################################################################!
!==============================================================================================================================!
!=========================                                 Computation                                =========================!
!==============================================================================================================================!
!##############################################################################################################################!



print *
print *, 'Computation...'


!===========================================================!
! - - - - - - - - - -  initialisation:  - - - - - - - - - - !

nvisited = 0
nconnected = 0

! "list_visited" contains {i,j} positions in of visited points (shape of "list_visited": (npoints,2))
! "nvisited" is the number of visited points

visited = .false.
!       "visited(i,j)" [logical] tells if the point {i,j} have been visited

connected = .false.

last_border_pos = -9d+99
!       contains the position on a cell edge of the right-most (looking from outside the cell) entry/exist point
!       ** A border (or edge) is considered empty is its right-most point position is -9d+99 **

! "is_last_border_out" tells if the the point of each edge of each cell is an exit point (= .true.) or an entry point
! (= .false.) 

!  "last_border_pos" and "is_last_border_out" have the shape (4,nx+1,ny)
!  1st dimension: cell edges: 1: left, 2: up, 3: right, 4: down

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!===========================================================!


! Subdivision in NTASK groups:
! ----------------------------
list_task = (/ ((k*nrec)/NTASK, k = 0,NTASK) /)
! make sure not to miss any record:
list_task(1) = 0
list_task(NTASK+1) = nrec


! Display head of progression message
! -----------------------------------
if (print_progress) then
  print *
  print *, 'Progress:'
  write(*, fmt=barformat) barprogress
end if


!-----------------------------------------!
! ====== Loop on shapefile records ====== !
!-----------------------------------------!
totarea = 0d0
totrastarea = 0d0
n = 0
do itask = 1,NTASK
    do irec = list_task(itask)+1,list_task(itask+1)


        !---------------------------------!
        ! ==== Loop on polygon parts ==== !
        !---------------------------------!
        recarea = 0d0
        do ipart = 1,nparts(irec)

            n = n + 1

            npts = partlen(n)
            pol_part(:,1:npts) = polygon(:, recpos(irec)+partpos(n):recpos(irec)+partpos(n)+npts-1)


            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
            ! - Preliminary loop on points of current polygon part to have a beginning point on a cell edge:  - !

            k0 = 1
            x2 = pol_part(1,k0)
            y2 = pol_part(2,k0)
            call get_pos(pol_part(:,k0), pol_part(:,k0+1), xbnd, ybnd, i, j)
            y0 = y2!ybnd(j)

            area = 0d0
            goon = .true.
            do while (goon)

                x1 = x2
                y1 = y2
                x2 = pol_part(1,k0+1)
                y2 = pol_part(2,k0+1)

                area = area + trapez_area(y0, x1, y1, x2, y2, authradius)

                if (y2 > ybnd(j+1)) then
                    dj = +1
                    naddj = 1
                    do while (y2 > ybnd(j+1+naddj))
                        naddj = naddj + 1
                    end do
                    j = j + naddj
                    yaddtiming(1) = (ybnd(j) - y1) / (y2 - y1)
                    xin = x1 + yaddtiming(1)*(x2 - x1)
                    goon = .false.
                elseif (y2 < ybnd(j)) then
                    dj = -1
                    naddj = 1
                    do while (y2 < ybnd(j-naddj))
                        naddj = naddj + 1
                    end do
                    j = j - naddj
                    yaddtiming(1) = (ybnd(j+1) - y1) / (y2 - y1)
                    xin = x1 + yaddtiming(1)*(x2 - x1)
                    goon = .false.
                else
                    dj = 0
                    yaddtiming(1) = -1d99
                end if

                if (x2 > xbnd(i+1)) then
                    di = +1
                    naddi = 1
                    do while (x2 > xbnd(i+1+naddi))
                        naddi = naddi + 1
                    end do
                    i = i + naddi
                    xaddtiming(1) = (xbnd(i) - x1) / (x2 - x1)
                    yin = y1 + xaddtiming(1)*(y2 - y1)
                    goon = .false.
                elseif (x2 < xbnd(i)) then
                    di = -1
                    naddi = 1
                    do while (x2 < xbnd(i-naddi))
                        naddi = naddi + 1
                    end do
                    i = i - naddi
                    xaddtiming(1) = (xbnd(i+1) - x1) / (x2 - x1)
                    yin = y1 + xaddtiming(1)*(y2 - y1)
                    goon = .false.
                else
                    di = 0
                    xaddtiming(1) = -1d99
                end if

                k0 = k0 + 1
                if (k0 == npts) then
                    goon = .false.
                end if

            end do
            ! - End of preliminary loop - !
            ! - - - - - - - - - - - - - - !


            if (.not. visited(i,j)) then
                visited(i,j) = .true.
                nvisited = nvisited + 1
                list_visited(:,nvisited) = (/i,j/)
            end if


            if (k0 == npts) then ! The entire polygon part is in the same cell

                polygon_area(i,j) = polygon_area(i,j) + area


            else


                if (check) then
                    ! Continute Preliminary loop to the end of the records to get its area
                    x2 = pol_part(1,k0)
                    y2 = pol_part(2,k0)
                    do k = k0,npts-1
                        x1 = x2
                        y1 = y2
                        x2 = pol_part(1,k+1)
                        y2 = pol_part(2,k+1)
                        area = area + trapez_area(y0, x1, y1, x2, y2, authradius)
                    end do
                end if
                

                ! new starting point:
                if (xaddtiming(1) > yaddtiming(1)+tolerance) then
                    ! consider 'x' shift if last jump was across one of cell 'x' boundaries
                    xin = xbnd(i + int((1-di)/2))
                    dj = 0
                elseif (yaddtiming(1) > xaddtiming(1)+tolerance) then
                    ! consider 'y' shift if last jump was across one of cell 'y' boundaries
                    yin = ybnd(j + int((1-dj)/2))
                    di = 0
                else
                    ! if last jump crossed one of cell corner (xaddtiming == yaddtiming +/- tolerance)
                    ! consider 'x' shift if negative diagonal motion (di = -dj)  or 'y' shift is positive diagonal motion
                    xin = xbnd(i + int((1-di)/2))
                    yin = ybnd(j + int((1-dj)/2))
                    if (di == -dj) then
                        iadd = iadd + 1
                    else
                        jadd = jadd + 1
                    end if
                end if


                ! Reorder polygon part with the new starting point
                pol_part(1, 1:npts) = (/ pol_part(1, k0:npts), pol_part(1, 2:k0) /)
                pol_part(2, 1:npts) = (/ pol_part(2, k0:npts), pol_part(2, 2:k0) /)

                ! Loop initial conditions:
                dj0 = dj
                x1 = xin
                y1 = yin
                x2 = pol_part(1,1)
                y2 = pol_part(2,1)
                y0 = ybnd(j)

                !------------------------------------!
                ! == Loop on polygon part points: == !
                !------------------------------------!
                do k = 2,npts

                    ! record area of last segment
                    polygon_area(i,j) = polygon_area(i,j)  + trapez_area(y0, x1, y1, x2, y2, authradius)

                    x1 = x2
                    y1 = y2
                    x2 = pol_part(1,k)
                    y2 = pol_part(2,k)

                    ! If (x,y)_k+1 is out of the current cell: CROSSING EVENT
                    ! --> define an additional point for every time the segment [(x,y)_k, (x,y)_k+1]
                    ! crosses a cell edge (both "x-direction" and "y-direction" edges).
                    ! The "timing" of the crossing is the position of the crossing point
                    ! on the segment [(x,y)_k, (x,y)_k+1] (i.e., 0 at (x,y)_k, 1 at (x,y)_k+1)
                    if (y2 > ybnd(j+1)) then
                        dj = +1
                        naddj = 1
                        yadd(naddj) = ybnd(j+naddj)
                        yaddtiming(naddj) = (yadd(naddj) - y1) / (y2 - y1)
                        do while (y2 > ybnd(j+1+naddj))
                            naddj = naddj + 1
                            yadd(naddj) = ybnd(j+naddj)
                            yaddtiming(naddj) = (yadd(naddj) - y1) / (y2 - y1)
                        end do
                    elseif (y2 < ybnd(j)) then
                        dj = -1
                        naddj = 1
                        yadd(naddj) = ybnd(j-naddj+1)
                        yaddtiming(naddj) = (yadd(naddj) - y1) / (y2 - y1)
                        do while (y2 < ybnd(j-naddj))
                            naddj = naddj + 1
                            yadd(naddj) = ybnd(j-naddj+1)
                            yaddtiming(naddj) = (yadd(naddj) - y1) / (y2 - y1)
                        end do
                    else
                        if (y1==ybnd(j+1) .and. y2==ybnd(j+1) .and. x2<xbnd(i)) then
                            dj = +1
                            naddj = 1
                            yadd(1) = ybnd(j+1)
                            yaddtiming(1) = 0d0
                        elseif (y1==ybnd(j) .and. y2==ybnd(j) .and. x2>xbnd(i+1)) then
                            dj = -1
                            naddj = 1
                            yadd(1) = ybnd(j)
                            yaddtiming(1) = 0d0
                        else
                            dj = 0
                            naddj = 0
                        end if
                    end if

                    if (x2 > xbnd(i+1)) then
                        di = +1
                        naddi = 1
                        xadd(naddi) = xbnd(i+naddi)
                        xaddtiming(naddi) = (xadd(naddi) - x1) / (x2 - x1)
                        do while (x2 > xbnd(i+1+naddi))
                            naddi = naddi + 1
                            xadd(naddi) = xbnd(i+naddi)
                            xaddtiming(naddi) = (xadd(naddi) - x1) / (x2 - x1)
                        end do
                    elseif (x2 < xbnd(i)) then
                        di = -1
                        naddi = 1
                        xadd(naddi) = xbnd(i-naddi+1)
                        xaddtiming(naddi) = (xadd(naddi) - x1) / (x2 - x1)
                        do while (x2 < xbnd(i-naddi))
                            naddi = naddi + 1
                            xadd(naddi) = xbnd(i-naddi+1)
                            xaddtiming(naddi) = (xadd(naddi) - x1) / (x2 - x1)
                        end do
                    else
                        if (x1==xbnd(i+1) .and. x2==xbnd(i+1) .and. y2>ybnd(j+1)) then
                            di = +1
                            naddi = 1
                            xadd(1) = xbnd(i+1)
                            xaddtiming(1) = 0d0
                        elseif (x1==xbnd(i) .and. x2==xbnd(i) .and. y2<ybnd(j)) then
                            di = -1
                            naddi = 1
                            xadd(1) = xbnd(i)
                            xaddtiming(1) = 0d0
                        else
                            di = 0
                            naddi = 0
                        end if
                    end if

                    ! In case next point (x,y)_k+1 jumped across one or several cell edges (i.e., naddi>0 or naddj>0),
                    ! determine the order of the crossings of the "x" and "y" cell edges between (x,k)_k and (x,y)_k+1
                    ! by sorting the "timings" (= position of the crossing points on the segment [(x,k)_k, (x,y)_k+1])
                    xaddtiming(naddi+1) = 1d99
                    yaddtiming(naddj+1) = 1d99
                    iadd = 1
                    jadd = 1
                    nadd = 0
                    do while (iadd <= naddi .or. jadd <= naddj)
                        nadd = nadd + 1
                        if (xaddtiming(iadd) < yaddtiming(jadd)-tolerance) then
                            ! consider 'x' shift if current jump crossed one of cell 'x' boundaries first
                            xyadd(:,nadd) = (/ xadd(iadd), y1 + xaddtiming(iadd)*(y2-y1) /)
                            iadd = iadd + 1
                        elseif (yaddtiming(jadd) < xaddtiming(iadd)-tolerance) then
                            ! consider 'y' shift if current jump crossed one of cell 'y' boundaries first
                            xyadd(:,nadd) = (/ x1 + yaddtiming(jadd)*(x2-x1), yadd(jadd) /)
                            jadd = jadd + 1
                        else
                            ! if jump crossed one of cell corners (xaddtiming == yaddtiming +/- tolerance)
                            ! consider 'x' shift if positive diagonal motion (di==dj) or 'y' shift if neg. diagonal motion (di==-dj)
                            xyadd(:,nadd) = (/ xadd(iadd), yadd(jadd) /)
                            if (di == dj) then
                                iadd = iadd + 1
                            else
                                jadd = jadd + 1
                            end if
                        end if
                        ijadd(:,nadd) = (/ i+di*(iadd-1), j+dj*(jadd-1) /)
                    end do


                    ! record the trapezium areas of the additional points (cell edge crossings),
                    ! and record in/out points on borders of the cells visited:
                    do kadd = 1,nadd

                        x2_add = xyadd(1,kadd)
                        y2_add = xyadd(2,kadd)
                        di = ijadd(1,kadd) - i
                        dj = ijadd(2,kadd) - j

                        polygon_area(i,j) = polygon_area(i,j)  + trapez_area(y0, x1, y1, x2_add, y2_add, authradius)

                        call close_cell_ring(polygon_area(i,j), xbnd(i:i+1), ybnd(j:j+1), xin, x2_add, y0, dj0, dj, authradius)

                        call record_border_out(di, dj, x2_add, y2_add, last_border_pos(:,i,j), is_last_border_out(:,i,j))

                        dj0 = dj
                        i = ijadd(1,kadd)
                        j = ijadd(2,kadd)
                        x1 = xyadd(1,kadd)
                        xin = x1
                        y1 = xyadd(2,kadd)
                        yin = y1
                        y0 = ybnd(j)

                        if (.not. visited(i,j)) then
                            visited(i,j) = .true.
                            nvisited = nvisited + 1
                            list_visited(:,nvisited) = (/i,j/)
                        end if

                        call record_border_in(di, dj, x2_add, y2_add, last_border_pos(:,i,j), is_last_border_out(:,i,j))

                    end do

                end do
                !-----------------------------!
                ! == end of loop on points == !
                !-----------------------------!


            end if
            recarea = recarea + area


        end do
        !--------------------------------!
        ! ==== end of loop on parts ==== !
        !--------------------------------!



        if (split) then


            ! areas modulo:
            do k = 1,nvisited
                i = list_visited(1,k)
                j = list_visited(2,k)
                do while (polygon_area(i,j)  <  -tolerance * grid_area(i,j))
                    polygon_area(i,j) = polygon_area(i,j) + grid_area(i,j)
                end do
                do while (polygon_area(i,j)  >  (1d0+tolerance) * grid_area(i,j))
                    polygon_area(i,j) = polygon_area(i,j) - grid_area(i,j)
                end do
                if (polygon_area(i,j) < 0d0)             polygon_area(i,j) = 0d0
                if (polygon_area(i,j) > grid_area(i,j))  polygon_area(i,j) = grid_area(i,j)
            end do


            ! fill polygon interior:
            call get_connected_points(list_visited, nvisited, last_border_pos, is_last_border_out, &
                                      connected, list_connected, nconnected)
            call fill_polygon_interior(polygon_area, grid_area, visited, connected, list_connected, nconnected)

            ! combine duplicated cells
            polygon_area(1,:) = polygon_area(1,:) + polygon_area(nx+1,:)


            if (split_records) then
                ! Write output
                if (area_fract) then
                    call write_output_slice(ofid, ovarid, nx, ny, flip_x, flip_y, ishift, polygon_area(1:nx,:)/grid_area(1:nx,:), &
                                            irec)
                else
                    call write_output_slice(ofid, ovarid, nx, ny, flip_x, flip_y, ishift, polygon_area(1:nx,:), irec)
                end if
            else
                ! Add current record area to total
                polygon_area_sum = polygon_area_sum + polygon_area(1:nx,:)
            end if


            if (check) then
                if (check_rec_area) then
                    rastarea = sum(polygon_area(1:nx,:))
                    totrastarea = totrastarea + rastarea
                    ! Display summary
                    print *
                    print *, '---------------------------------------------------'
                    print *, 'Record #', irec, '/', nrec
                    print *, 'Record area (m2):          ', recarea
                    print *, 'Record area in raster (m2):', rastarea
                    print *, 'Minimum area ratio:        ', minval(polygon_area(1:nx,:)/grid_area(1:nx,:))
                    print *, 'Maximum area ratio:        ', maxval(polygon_area(1:nx,:)/grid_area(1:nx,:))
                    print *, '---------------------------------------------------'
                end if
                totarea = totarea + recarea
            end if


            ! reinitialization
            call reinitialize_border(list_visited, nvisited, last_border_pos)
            polygon_area = 0d0
            nvisited = 0
            nconnected = 0
            visited = .false.
            connected = .false.



        else

            if (split_connected_computation) then
                ! Record connected points:
                call get_connected_points(list_visited, nvisited, last_border_pos, is_last_border_out, &
                                          connected, list_connected, nconnected)
                call reinitialize_border(list_visited, nvisited, last_border_pos)
            end if

            if (check)  totarea = totarea + recarea

        end if



    end do
    if (print_progress) write(*, fmt='(A1)', advance='no') '|'
end do
if (print_progress) print *
!------------------------------------------------!
! ====== End of loop on shapefile records ====== !
!------------------------------------------------!



if (.not. split_records) then


    if (split_computation) then
        polygon_area(1:nx,:) = polygon_area_sum
    else

        ! areas modulo:
        do k = 1,nvisited
            i = list_visited(1,k)
            j = list_visited(2,k)
            do while (polygon_area(i,j)  <  -tolerance * grid_area(i,j))
                polygon_area(i,j) = polygon_area(i,j) + grid_area(i,j)
            end do
            do while (polygon_area(i,j)  >  (1d0+tolerance) * grid_area(i,j))
                polygon_area(i,j) = polygon_area(i,j) - grid_area(i,j)
            end do
            if (polygon_area(i,j)  <  tolerance * grid_area(i,j))        polygon_area(i,j) = 0d0
            if (polygon_area(i,j)  >  (1d0-tolerance) * grid_area(i,j))  polygon_area(i,j) = grid_area(i,j)
        end do

        ! fill polygon interior:
        if (.not. split_connected_computation)  call get_connected_points(list_visited, nvisited, last_border_pos, &
                                                                          is_last_border_out, connected, list_connected, nconnected)
        call fill_polygon_interior(polygon_area, grid_area, visited, connected, list_connected, nconnected)

        ! combine duplicated cells
        polygon_area(1,:) = polygon_area(1,:) + polygon_area(nx+1,:)

    end if


    ! Write output
    if (area_fract) then
        call write_output(ofid, ovarid, nx, ny, flip_x, flip_y, ishift, polygon_area(1:nx,:)/grid_area(1:nx,:))
    else
        call write_output(ofid, ovarid, nx, ny, flip_x, flip_y, ishift, polygon_area(1:nx,:))
    end if


    totrastarea = sum(polygon_area(1:nx,:))


end if



if (check) then
    ! Display summary
    print *
    print *, '---------------------------------------------------------'
    print *, 'Earth area (m2):                 ', sum(grid_area(1:nx,:))
    print *, 'Shapefile Total area (m2):       ', totarea
    print *, 'Total area of raster output (m2):', totrastarea
    print *, 'Minimum area ratio:              ', minval(polygon_area(1:nx,:)/grid_area(1:nx,:))
    print *, 'Maximum area ratio:              ', maxval(polygon_area(1:nx,:)/grid_area(1:nx,:))
    print *, '---------------------------------------------------------'
end if




!===================!
! Close output file !
!===================!

call close_file(ofid)

print *
print *, 'Created '//trim(output_file)//' successfully.'




end program
