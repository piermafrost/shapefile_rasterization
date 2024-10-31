module shp2grid_functions
implicit none

double precision, parameter :: DEG2RAD = acos(-1d0)/180d0

integer, parameter          :: iNA = -1
double precision, parameter :: NA  = -6d66

contains


    !==================================================!


    function get_pos_regular(x, x0, dx)
        integer:: get_pos_regular
        double precision, intent(in):: x, x0, dx
        if (x==x0) then
            get_pos_regular = 1
        else
            get_pos_regular = ceiling((x-x0)/dx)
        end if
    end function

    !----------!

    function get_pos_dichotomy(x, vec)
        integer:: get_pos_dichotomy
        double precision, intent(in):: x, vec(:)
        integer:: i1, i2, i

        i2 = size(vec)
        i1 = 1
        do while (i2-i1 > 1)
            i = (i1+i2) / 2 
            !
            if (x >= vec(i)) then
                i1 = i
            else
                i2 = i
            end if
        end do

        get_pos_dichotomy = i1

    end function

    !----------!

    subroutine get_pos(xy, xy_next, xbnd, ybnd, i, j)
        double precision, dimension(2), intent(in):: xy, xy_next
        double precision, dimension(:), intent(in):: xbnd, ybnd
        integer, intent(out):: i, j

        i = get_pos_dichotomy(xy(1), xbnd)
        j = get_pos_dichotomy(xy(2), ybnd)

        if (xy(1)==xbnd(i) .and. xy_next(1)==xbnd(i) .and. xy_next(2)<xy(2))   i = i - 1
        if (xy(2)==ybnd(j) .and. xy_next(2)==ybnd(j) .and. xy_next(1)>xy(1))   j = j - 1

    end subroutine


    !==================================================!


    function trapez_area(y0, x1, y1, x2, y2, RE)
        double precision:: trapez_area
        double precision, intent(in):: y0, x1, y1, x2, y2, RE
        if (abs(y1-y2)<1d-14) then
            trapez_area = (RE**2) * (x2-x1) * (sin((y1+y2)/2) - sin(y0))
        else
            trapez_area = (RE**2) * (x2-x1) * ((cos(y1)-cos(y2))/(y2-y1) - sin(y0))
        end if
        !trapez_area = (x2-x1) * ((y1+y2)/2 - y0)
    end function


    !==================================================!

    subroutine close_cell_ring(pol_cell_area, cell_xbnd, cell_ybnd, xin, x2, y0, dj0, dj, RE)
        ! correct polygon area in 1 grid cell by connecting the entry and exit point of the polygon (ie: close the ring)
        double precision, intent(inout):: pol_cell_area
        double precision, dimension(2), intent(in):: cell_xbnd, cell_ybnd
        double precision, intent(in):: xin, x2, y0, RE
        integer, intent(in):: dj0, dj
        !
        if (dj==+1 .or. dj0==-1) then
            if (dj==-1) then
                pol_cell_area = pol_cell_area  +  trapez_area(y0, cell_xbnd(1), cell_ybnd(2), xin, cell_ybnd(2), RE)
            elseif (dj0==+1) then
                pol_cell_area = pol_cell_area  +  trapez_area(y0, x2, cell_ybnd(2), cell_xbnd(2), cell_ybnd(2), RE)
            else
                pol_cell_area = pol_cell_area  +  trapez_area(y0, x2, cell_ybnd(2), xin, cell_ybnd(2), RE)
            end if
        end if
        !
    end subroutine

    !==================================================!


!!!! TO: OPTIMIZE THIS FUNCTION
    subroutine get_ibord_in(di, dj, x, y, ibord, pos)
        integer, intent(in):: di, dj
        double precision, intent(in):: x, y
        integer, intent(out):: ibord
        double precision, intent(out):: pos
        double precision, dimension(4):: posvec
        !
        ! di==+1 & dj==0  =>  ibord = 1 & pos = y
        ! di==0 & dj==-1  =>  ibord = 2 & pos = x
        ! di==-1 & dj==0  =>  ibord = 3 & pos = -y
        ! di==1 & dj==+1  =>  ibord = 4 & pos = -x
        posvec = (/y, x, -y, -x/)
        ibord = di * (2*(1 - di) + 2)  +  dj * (2*(1-dj) + 2)
        if (di==+1) then
            ibord = 1
            pos = y
        elseif (dj==-1) then
            ibord = 2
            pos = x
        elseif (di==-1) then
            ibord = 3
            pos = -y
        else !dj==+1
            ibord = 4
            pos = -x
        end if
        !
    end subroutine

    !----------!

    subroutine get_ibord_out(di, dj, x, y, ibord, pos)
        integer, intent(in):: di, dj
        double precision, intent(in):: x, y
        integer, intent(out):: ibord
        double precision, intent(out):: pos
        !
        if (di==-1) then
            ibord = 1
            pos = y
        elseif (dj==+1) then
            ibord = 2
            pos = x
        elseif (di==+1) then
            ibord = 3
            pos = -y
        else !dj==-1
            ibord = 4
            pos = -x
        end if
        !
    end subroutine

    !----------!

    subroutine record_border_in(di, dj, x, y, last_border_pos, is_last_border_out)
        integer, intent(in):: di, dj
        double precision, intent(in):: x, y
        double precision, dimension(4), intent(out):: last_border_pos
        logical, dimension(4), intent(out):: is_last_border_out
        integer:: ibord
        double precision:: pos

        call get_ibord_in(di, dj, x, y, ibord, pos)

        if (pos > last_border_pos(ibord)) then
            last_border_pos(ibord) = pos
            is_last_border_out(ibord) = .false.
        end if

    end subroutine

    !----------!

    subroutine record_border_out(di, dj, x, y, last_border_pos, is_last_border_out)
        integer, intent(in):: di, dj
        double precision, intent(in):: x, y
        double precision, dimension(4), intent(out):: last_border_pos
        logical, dimension(4), intent(out):: is_last_border_out
        integer:: ibord
        double precision:: pos

        call get_ibord_out(di, dj, x, y, ibord, pos)

        if (pos > last_border_pos(ibord)) then
            last_border_pos(ibord) = pos
            is_last_border_out(ibord) = .true.
        end if

    end subroutine


    !==================================================!


    subroutine reinitialize_border(list_visited, nvisited, last_border_pos)
        integer, dimension(:,:), intent(in):: list_visited
        integer, intent(in):: nvisited
        double precision, dimension(:,:,:), intent(inout):: last_border_pos
        integer:: i, j, k
        do k = 1,nvisited
            i = list_visited(1,k)
            j = list_visited(2,k)
            last_border_pos(:,i,j) = -9d99
        end do
    end subroutine


    !==================================================!


    subroutine connect_point(i, j, connected, list_connected_points, nconnected)
        integer, intent(in):: i, j
        logical, dimension(:,:), intent(inout):: connected
        integer, dimension(:,:), intent(inout):: list_connected_points
        integer, intent(inout):: nconnected
        integer:: nx, ny
        !
        nx = size(connected,1)
        ny = size(connected,2)
        !
        if (i>=1 .and. i<=nx .and. j>=1 .and. j<=ny) then
            if (.not. connected(i,j)) then
                connected(i,j) = .true.
                nconnected = nconnected + 1
                list_connected_points(:,nconnected) = (/i,j/)
            end if
        end if
        !
    end subroutine

    !----------!

    subroutine get_connected_points(list_visited, nvisited, last_border_pos, is_last_border_out, &
                                    connected, list_connected_points, nconnected)
        integer, dimension(:,:), intent(in):: list_visited
        integer, intent(in):: nvisited
        double precision, dimension(:,:,:), intent(in):: last_border_pos
        logical, dimension(:,:,:), intent(in):: is_last_border_out
        logical, dimension(:,:), intent(inout):: connected
        integer, dimension(:,:), intent(inout):: list_connected_points
        integer, intent(inout):: nconnected
        !
        integer, dimension(7), parameter:: modulolist = (/1,2,3,4,1,2,3/)
        integer, dimension(4), parameter:: prev_bord = (/4,1,2,3/)
        integer, dimension(2,4), parameter:: dirlist = reshape( (/  -1,0,  0,+1,  +1,0,  0,-1 /),   shape=(/2,4/) )
        !                                                           left,   up,   right, down
        integer:: ibord, ibord0, idummy, i0, j0, i, j, k
        logical:: enterloop
        logical, dimension(4):: connected_bord

        !nconnected = 0

        do k = 1,nvisited

            i0 = list_visited(1,k)
            j0 = list_visited(2,k)

            enterloop = .false.
            do ibord = 1,4
                if (last_border_pos(ibord,i0,j0) > -9d+99) then
                    ibord0 = ibord
                    enterloop = .true.
                end if
            end do


            if (enterloop) then

                ! sort list of entry/exit points:
                connected_bord = (/ .false., .false., .false., .false. /)

                i = i0 + dirlist(1,ibord0)
                j = j0 + dirlist(2,ibord0)
                call connect_point(i, j, connected, list_connected_points, nconnected)
                connected_bord(ibord0) = .true.

                do idummy = ibord0+1,ibord0+3

                    ibord = modulolist(idummy) ! ibord = ibord + 1 modulo 4, as modulolist = (/1,2,3,4,1,2,3/)

                    ! border is connected if there is any entry/exit point on it:
                    if (last_border_pos(ibord,i0,j0) > -9d99) then
                        i = i0 + dirlist(1,ibord)
                        j = j0 + dirlist(2,ibord)
                        call connect_point(i, j, connected, list_connected_points, nconnected)
                        connected_bord(ibord) = .true.

                    ! border connected if empty and...
                    elseif (last_border_pos(prev_bord(ibord),i0,j0) > -9d+99) then
                        ! ... last point of previous border is an exit point:
                        if (is_last_border_out(prev_bord(ibord),i0,j0)) then
                            i = i0 + dirlist(1,ibord)
                            j = j0 + dirlist(2,ibord)
                            call connect_point(i, j, connected, list_connected_points, nconnected)
                            connected_bord(ibord) = .true.
                        end if

                    else
                        !  ... or previous border is empty and connected:
                        if (connected_bord(prev_bord(ibord))) then
                            i = i0 + dirlist(1,ibord)
                            j = j0 + dirlist(2,ibord)
                            call connect_point(i, j, connected, list_connected_points, nconnected)
                            connected_bord(ibord) = .true.
                        end if

                    end if

                end do
            end if

        end do

    end subroutine

    !----------!

    subroutine fill_polygon_interior(polygon_area, grid_area, visited, connected, list_active, nactive)
        double precision, dimension(:,:), intent(inout):: polygon_area
        double precision, dimension(:,:), intent(in):: grid_area
        logical, dimension(:,:), intent(in):: visited
        logical, dimension(:,:), intent(inout):: connected
        integer, dimension(:,:), intent(inout):: list_active
        integer, dimension(:,:), allocatable:: new_list_active
        integer, intent(inout):: nactive
        integer, dimension(2,4), parameter:: dirlist = reshape( (/  -1,0,  0,+1,  +1,0,  0,-1 /),   shape=(/2,4/) )
        !                                                           left,   up,   right, down
        integer:: i, j, i0, j0, k, l, n, nx, ny

        nx = size(polygon_area,1)
        ny = size(polygon_area,2)

        allocate( new_list_active(2, size(list_active,2)) )

        do while (nactive > 0)

            n = 0

            do k = 1,nactive

                i0 = list_active(1,k)
                j0 = list_active(2,k)

                if (polygon_area(i0,j0)==0d0) then

                    polygon_area(i0,j0) = grid_area(i0,j0)

                    if (.not. visited(i0,j0)) then
                        do l = 1,4
                            i = i0 + dirlist(1,l)
                            j = j0 + dirlist(2,l)
                            call connect_point(i, j, connected, new_list_active, n)
                        end do
                    end if

                end if

            end do

            nactive = n
            list_active(:,1:nactive) = new_list_active(:,1:nactive)

        end do

    end subroutine


    !==================================================!


    subroutine define_grid(bbox, global_grid, nx, ny, dlon, dlat, RE, eccentricity, &
                           lon, lat, lonbounds, latbounds, xbounds, ybounds, grid_area, init_polygon_area)
        use geography, only: authalic_latitude
        double precision, dimension(4), intent(in):: bbox
        logical, intent(in):: global_grid
        integer, intent(inout):: nx, ny
        double precision, intent(inout):: dlon, dlat
        double precision, intent(in):: RE, eccentricity
        double precision, dimension(:), allocatable, intent(out):: lon, lat, lonbounds, latbounds, xbounds, ybounds
        double precision, dimension(:,:), allocatable, intent(out):: grid_area, init_polygon_area
        integer:: i, j

        if (nx/=iNA .and. ny/=iNA) then
            dlon = (bbox(3)-bbox(1))/nx
            dlat = (bbox(4)-bbox(2))/ny
        elseif (dlon/=NA .and. dlat/=NA) then
            nx = ceiling((bbox(3)-bbox(1)) / dlon)
            ny = ceiling((bbox(4)-bbox(2)) / dlat)
        else
            print *, 'ERROR: Cannot define grid. Need grid size or resolution information'
            stop
        end if

        allocate(lon(nx))
        allocate(lat(ny))
        allocate(lonbounds(nx+2)) ! --> extra "lon" cell to duplicate the 1st cell
        allocate(latbounds(ny+1))
        allocate(xbounds(nx+2))
        allocate(ybounds(ny+1))
        allocate(init_polygon_area(nx+1,ny))
        allocate(grid_area(nx+1,ny))

        lonbounds(1:nx+1) = (/ (bbox(1) + i*dlon,  i=0,nx) /)
        latbounds(1:ny+1) = (/ (bbox(2) + j*dlat,  j=0,ny) /)
        ! safeguard correction:
        latbounds(1)    = max(-90d0, latbounds(1))
        latbounds(ny+1) = min(90d0,  latbounds(ny+1))
        ! extra cell:
        if (global_grid) then
            lonbounds(nx+2) = lonbounds(2) + 360d0
        else
            lonbounds(nx+2) = lonbounds(nx+1) ! empty cell
        end if

        lon = (lonbounds(1:nx) + lonbounds(2:nx+1)) / 2
        lat = (latbounds(1:ny) + latbounds(2:ny+1)) / 2

        ! degree -> radian conversion:
        xbounds = DEG2RAD * lonbounds
        ybounds = DEG2RAD * latbounds

        ! conversion in authalic latitude:
        call authalic_latitude(ybounds , eccentricity)
        !xbounds = lonbounds
        !ybounds = latbounds

        ! Total grid cell area:
        do j = 1,ny
            do i = 1,nx+1
                grid_area(i,j)  =  (RE**2) * (xbounds(i+1)-xbounds(i)) * (sin(ybounds(j+1))-sin(ybounds(j)))
                !grid_area(i,j)  =  (xbounds(i+1)-xbounds(i)) * (ybounds(j+1)-ybounds(j))
            end do
        end do

        ! Initialise polygon area (zero)
        init_polygon_area = 0d0

    end subroutine


    !==================================================!


    subroutine get_grid_from_netcdf_file(grid_file, xvar_name, yvar_name, xbndvar_name, ybndvar_name, RE,eccentricity, global_grid,&
                   shpf_bbox, grid_bbox, nx, ny, lon, lat, lonbnd, latbnd, xbnd, ybnd, flip_x, flip_y, grid_area, init_polygon_area)
        use netcdf
        use netcdf_output, only: nf90check
        use geography, only: authalic_latitude

        character(len=*), intent(in):: grid_file, xvar_name, yvar_name, xbndvar_name, ybndvar_name
        double precision, intent(in):: RE, eccentricity
        logical, intent(inout):: global_grid
        double precision, dimension(4), intent(in):: shpf_bbox, grid_bbox
        integer, intent(out):: nx, ny
        double precision, dimension(:), allocatable, intent(out):: lon, lat, lonbnd, xbnd, latbnd, ybnd
        double precision, dimension(:,:), allocatable, intent(out):: grid_area, init_polygon_area
        logical, intent(out):: flip_x, flip_y
        
        double precision, dimension(:,:), allocatable :: loc_bnd
        integer :: fid, varid, ndims, xdimid(1), ydimid(1), dimids(2), ierr
        integer :: k, i, j, nb(2)
        logical :: swap_bnd


        ! Open input netcdf 'grid" file
        ! -----------------------------
        ierr = nf90_open(grid_file, NF90_NOWRITE, fid)
        call nf90check(ierr, 'Error while openning input file '//trim(grid_file))


        ! Get "lon" variable
        ! ------------------
        ierr = nf90_inq_varid(fid, trim(xvar_name), varid)
        call nf90check(ierr, 'Error while inquiring ID of variable "'//trim(xvar_name)//'" in input file '//trim(grid_file))
        ierr = nf90_inquire_variable(fid, varid, ndims=ndims)
        call nf90check(ierr, 'Error while inquiring number of dimensions of variable "'//trim(xvar_name)//'" in input file ' &
                             //trim(grid_file))
        if (ndims /= 1) then
            print *, 'ERROR: x axis variable "'//trim(xvar_name)//'" in input file '//trim(grid_file)//' must have 1 dimension.'
            print *, 'Got: ', ndims, 'dimensions.'
            stop 2
        end if
        ierr = nf90_inquire_variable(fid, varid, dimids=xdimid(1:1))
        call nf90check(ierr, 'Error while IDs of dimensions of variable "'//trim(xvar_name)//'" in input file '//trim(grid_file))
        ierr = nf90_inquire_dimension(fid, xdimid(1), len=nx)
        call nf90check(ierr, 'Error while inquiring length of dimension in input file '//trim(grid_file))
        allocate( lon(nx) )
        ierr = nf90_get_var(fid, varid, lon)
        call nf90check(ierr, 'Error while getting variable "'//trim(xvar_name)//'" of file '//trim(grid_file))
        !
        if (lon(nx) < lon(1)) then
            flip_x = .true.
            lon = lon(nx:1:-1)
        else
            flip_x = .false.
        end if

        ! Get "lat" variable
        ! ------------------
        ierr = nf90_inq_varid(fid, trim(yvar_name), varid)
        call nf90check(ierr, 'Error while inquiring ID of variable "'//trim(yvar_name)//'" in input file '//trim(grid_file))
        ierr = nf90_inquire_variable(fid, varid, ndims=ndims)
        call nf90check(ierr, 'Error while inquiring number of dimensions of variable "'//trim(yvar_name)//'" in input file ' &
                             //trim(grid_file))
        if (ndims /= 1) then
            print *, 'ERROR: y axis variable "'//trim(yvar_name)//'" in input file '//trim(grid_file)//' must have 1 dimension.'
            print *, 'Got: ', ndims, 'dimensions.'
            stop 2
        end if
        ierr = nf90_inquire_variable(fid, varid, dimids=ydimid(1:1))
        call nf90check(ierr, 'Error while IDs of dimensions of variable "'//trim(yvar_name)//'" in input file '//trim(grid_file))
        ierr = nf90_inquire_dimension(fid, ydimid(1), len=ny)
        call nf90check(ierr, 'Error while inquiring length of dimension in input file '//trim(grid_file))
        allocate( lat(ny) )
        ierr = nf90_get_var(fid, varid, lat)
        call nf90check(ierr, 'Error while getting variable "'//trim(yvar_name)//'" of file '//trim(grid_file))
        !
        if (lat(ny) < lat(1)) then
            flip_y = .true.
            lat = lat(ny:1:-1)
        else
            flip_y = .false.
        end if


        allocate( lonbnd(nx+2) ) ! --> extra cell to duplicate 1st cell
        allocate(   xbnd(nx+1) )
        allocate( latbnd(ny+1) )
        allocate(   ybnd(ny+1) )


        ! Get "lon bounds" variable
        ! -------------------------
        if (trim(xbndvar_name) /= '') then ! -- > get variable from file
            !
            ierr = nf90_inq_varid(fid, trim(xbndvar_name), varid)
            call nf90check(ierr, 'Error while inquiring ID of variable "'//trim(xbndvar_name)//'" in input file '//trim(grid_file))
            ierr = nf90_inquire_variable(fid, varid, ndims=ndims)
            call nf90check(ierr, 'Error while inquiring number of dimensions of variable "'//trim(xbndvar_name) &
                                 //'" in input file '//trim(grid_file))
            if (ndims /= 2) then
                print *, 'ERROR: x bounds variable "'//trim(xbndvar_name)//'" in input file '//trim(grid_file) &
                         //'must have 2 dimensions.'
                print *, 'Got: ', ndims, 'dimensions.'
                stop 2
            end if
            ierr = nf90_inquire_variable(fid, varid, dimids=dimids(1:2))
            call nf90check(ierr, 'Error while IDs of dimensions of variable "'//trim(xbndvar_name)//'" in input file ' &
                                  //trim(grid_file))
            do k = 1,2
                ierr = nf90_inquire_dimension(fid, dimids(k), len=nb(k))
                call nf90check(ierr, 'Error while inquiring length of dimension in input file '//trim(grid_file))
            end do
            if (dimids(1)==xdimid(1) .and. nb(2)==2) then
                swap_bnd = .false.
                allocate( loc_bnd(nx,2) )
            elseif (dimids(2)==xdimid(1) .and. nb(1)==2) then
                swap_bnd = .true.
                allocate( loc_bnd(2,nx) )
            else
                print *, 'ERROR: inconsistent "x bound" variable.'
                print *, 'Must be 2-dimensional, 1 dimension being the same as "x" variable, the second being length-2.'
                stop 2
            end if
            ierr = nf90_get_var(fid, varid, loc_bnd)
            call nf90check(ierr, 'Error while getting variable "'//trim(xbndvar_name)//'" of file '//trim(grid_file))
            !
            if (swap_bnd) then
                lonbnd(2:nx) = (loc_bnd(1,2:nx) + loc_bnd(2,1:nx-1)) / 2
                lonbnd(1)    = loc_bnd(1,1)
                lonbnd(nx+1) = loc_bnd(2,nx)
            else
                lonbnd(2:nx) = (loc_bnd(2:nx,1) + loc_bnd(1:nx-1,2)) / 2
                lonbnd(1)    = loc_bnd(1,1)
                lonbnd(nx+1) = loc_bnd(nx,2)
            end if
            deallocate( loc_bnd )
            !
            if (flip_x) lonbnd(1:nx+1) = lonbnd(nx+1:1:-1)
            !
        else ! --> default grid boundaries: middle of grid points
            lonbnd(2:nx) = (lon(1:nx-1) + lon(2:nx)) / 2
            lonbnd(1)    = lon(1)  - (lonbnd(2)-lon(1))
            lonbnd(nx+1) = lon(nx) + (lon(nx)-lonbnd(nx-1))
        end if


        ! Get "lat bounds" variable
        ! -------------------------
        if (trim(ybndvar_name) /= '') then ! -- > get variable from file
            !
            ierr = nf90_inq_varid(fid, trim(ybndvar_name), varid)
            call nf90check(ierr, 'Error while inquiring ID of variable "'//trim(ybndvar_name)//'" in input file '//trim(grid_file))
            ierr = nf90_inquire_variable(fid, varid, ndims=ndims)
            call nf90check(ierr, 'Error while inquiring number of dimensions of variable "'//trim(ybndvar_name) &
                                 //'" in input file '//trim(grid_file))
            if (ndims /= 2) then
                print *, 'ERROR: y bounds variable "'//trim(ybndvar_name)//'" in input file '//trim(grid_file) &
                         //'must have 2 dimensions.'
                print *, 'Got: ', ndims, 'dimensions.'
                stop 2
            end if
            ierr = nf90_inquire_variable(fid, varid, dimids=dimids(1:2))
            call nf90check(ierr, 'Error while IDs of dimensions of variable "'//trim(ybndvar_name)//'" in input file ' &
                                  //trim(grid_file))
            do k = 1,2
                ierr = nf90_inquire_dimension(fid, dimids(k), len=nb(k))
                call nf90check(ierr, 'Error while inquiring length of dimension in input file '//trim(grid_file))
            end do
            if (dimids(1)==ydimid(1) .and. nb(2)==2) then
                swap_bnd = .false.
                allocate( loc_bnd(ny,2) )
            elseif (dimids(2)==ydimid(1) .and. nb(1)==2) then
                swap_bnd = .true.
                allocate( loc_bnd(2,ny) )
            else
                print *, 'ERROR: inconsistent "y bound" variable.'
                print *, 'Must be 2-dimensional, 1 dimension being the same as "y" variable, the second being length-2.'
                stop 2
            end if
            ierr = nf90_get_var(fid, varid, loc_bnd)
            call nf90check(ierr, 'Error while getting variable "'//trim(ybndvar_name)//'" of file '//trim(grid_file))
            !
            if (swap_bnd) then
                latbnd(2:ny) = (loc_bnd(1,2:ny) + loc_bnd(2,1:ny-1)) / 2
                latbnd(1)    = loc_bnd(1,1)
                latbnd(ny+1) = loc_bnd(2,ny)
            else
                latbnd(2:ny) = (loc_bnd(2:ny,1) + loc_bnd(1:ny-1,2)) / 2
                latbnd(1)    = loc_bnd(1,1)
                latbnd(ny+1) = loc_bnd(ny,2)
            end if
            deallocate( loc_bnd )
            !
            if (flip_y) latbnd = latbnd(ny+1:1:-1)
            !
        else ! --> default grid boundaries: middle of grid points
            latbnd(2:ny) = (lat(1:ny-1) + lat(2:ny)) / 2
            latbnd(1)    = max(-90d0,  lat(1)  - (latbnd(2)-lat(1)))
            latbnd(ny+1) = min(90d0,   lat(ny) + (lat(ny)-latbnd(ny-1)))
        end if


        ! Close input netcdf 'grid" file
        ! ------------------------------
        ierr = nf90_close(fid)
        call nf90check(ierr, 'Error while closing input file '//trim(grid_file))


        if (global_grid) then

            latbnd(1) = -90d0
            latbnd(ny+1) = 90d0
            lonbnd(1) = (lon(1) + lon(nx)-360d0)/2
            lonbnd(nx+1) = lonbnd(1) + 360d0
            lonbnd(nx+2) = lonbnd(2) + 360d0 ! extra cell duplicating the 1st cell

            ! do not check the shapefile bounding-box, for the grid will be aligned to it by the main program

        elseif (latbnd(1)==-90d0 .and. latbnd(ny+1)==90d0 .and. lonbnd(nx+1)==lonbnd(1)+360d0) then
            ! --> consider "global grid" behaviour if the grid loaded from the netCDF file matches the conditions

            global_grid = .true.
            lonbnd(nx+1) = lonbnd(1) + 360d0
            lonbnd(nx+2) = lonbnd(2) + 360d0 ! extra cell duplicating the 1st cell

        else

            ! add optional bounding-box information
            if ( (grid_bbox(1)/=NA .and. any(lon<grid_bbox(1))) .or. &
                (grid_bbox(2)/=NA .and. any(lat<grid_bbox(2))) .or. &
                (grid_bbox(3)/=NA .and. any(lon>grid_bbox(3))) .or. &
                (grid_bbox(4)/=NA .and. any(lat>grid_bbox(4))) ) then
                print *, 'ERROR: '//trim(xvar_name)//' and '//trim(yvar_name)//' axes loaded from'
                print *, 'netCDF file '//trim(grid_file)
                print *, 'exceed specified bounding box:', grid_bbox
                print *, 'Program stopped'
                stop 2
            end if
            if (grid_bbox(1)/=NA) lonbnd(1)    = grid_bbox(1)
            if (grid_bbox(2)/=NA) latbnd(1)    = grid_bbox(2)
            if (grid_bbox(3)/=NA) lonbnd(nx+1) = grid_bbox(3)
            if (grid_bbox(4)/=NA) latbnd(ny+1) = grid_bbox(4)
            lonbnd(nx+2) = lonbnd(nx+1) ! empty extra cell

            if ( shpf_bbox(1)<lonbnd(1) .or. shpf_bbox(2)<latbnd(1) .or. &
                 shpf_bbox(3)>lonbnd(nx+1) .or. shpf_bbox(4)>latbnd(ny+1) ) then
                print *, 'ERROR: Shapefile bounding-box exceeds the grid loaded from netCDF file'
                print *, trim(grid_file)
                print *, 'Shapefile bbox:        ', shpf_bbox
                print *, 'netCDF grid file bbox: ', lonbnd(1), latbnd(1), lonbnd(nx+1), latbnd(ny+1)
                print *, 'Program stopped'
                stop 2
            end if

         end if

        ! Conversion in radians
        xbnd = DEG2RAD*lonbnd
        ybnd = DEG2RAD*latbnd

        ! conversion in authalic latitude:
        call authalic_latitude(ybnd , eccentricity)

        allocate(init_polygon_area(nx+1,ny))
        allocate(grid_area(nx+1,ny))

        ! Total grid cell area:
        do j = 1,ny
            do i = 1,nx+1
                grid_area(i,j)  =  (RE**2) * (xbnd(i+1)-xbnd(i)) * (sin(ybnd(j+1))-sin(ybnd(j)))
            end do
        end do

        ! Initialise polygon area (zero)
        init_polygon_area = 0d0


    end subroutine


    !==================================================!



end module

