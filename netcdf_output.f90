module netcdf_output
use netcdf
implicit none

contains


    subroutine create_output_file(fname, x, y, xbounds, ybounds, flip_x, flip_y, fid, varid, nrec, earth_model, radius, flattening,&
                                  grid_area, area_fraction, ifname)

        use miscellaneous_functions, only: utc_time_string
        use geography, only: get_earth_parameters

        character(len=*), intent(in):: fname
        double precision, dimension(:), intent(in):: x, y, xbounds, ybounds
        logical, intent(in):: flip_x, flip_y
        integer, intent(out):: fid, varid
        integer, intent(in), optional:: nrec
        character(len=*), intent(in), optional:: earth_model
        double precision, intent(in), optional:: radius, flattening
        double precision, dimension(:,:), intent(in), optional:: grid_area
        logical, intent(in), optional:: area_fraction
        character(len=*), optional:: ifname

        double precision:: loc_radius, loc_flattening, eccentricity
        logical:: fract, write_earth_param
        character(len=1000):: command_line
        integer:: dimid(4), loc_varid(7), nx, ny, ierr, i



        if (present(area_fraction)) then
            fract = area_fraction
        else
            fract = .false.
        end if


        ! Data dimension:
        ! ---------------

        nx = size(x)
        ny = size(y)


        ! Create output file
        ! ------------------

        ierr = nf90_create(fname, NF90_CLOBBER, fid)
        call nf90check(ierr, 'Error while creating netCDF output file '//trim(fname))


        ! Global attributes
        ! ----------------------

        ierr = nf90_put_att(fid, NF90_GLOBAL, 'title', 'Area fraction of polygons on a longitude-latitude grid')
        call nf90check(ierr, 'Error while putting attribute "title" in '//trim(fname))

        if (present(earth_model)) then

            call get_earth_parameters(earth_model, loc_radius, loc_flattening, eccentricity)

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_model', earth_model)
            call nf90check(ierr, 'Error while putting attribute "earth_model" in '//trim(fname))

            write_earth_param = .true.

        elseif (present(radius) .and. present(flattening)) then

            loc_radius = radius
            loc_flattening = flattening
            eccentricity = sqrt(1d0 - (1d0-loc_flattening)**2)

            write_earth_param = .true.

        else

            write_earth_param = .false.

        end if

        if (write_earth_param) then

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_equatorial_radius', loc_radius)
            call nf90check(ierr, 'Error while putting attribute "earth_equatorial_radius" in '//trim(fname))

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_equatorial_radius_units', 'm')
            call nf90check(ierr, 'Error while putting attribute "earth_equatorial_radius_units" in '//trim(fname))

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_polar_radius', loc_radius*(1-loc_flattening))
            call nf90check(ierr, 'Error while putting attribute "earth_equatorial_radius" in '//trim(fname))

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_polar_radius_units', 'm')
            call nf90check(ierr, 'Error while putting attribute "earth_polar_radius_units" in '//trim(fname))

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_flattening', loc_flattening)
            call nf90check(ierr, 'Error while putting attribute "earth_flattening" in '//trim(fname))

            ierr = nf90_put_att(fid, NF90_GLOBAL, 'earth_first_eccentricity', eccentricity)
            call nf90check(ierr, 'Error while putting attribute "earth_eccentricity" in '//trim(fname))

        end if

        if (present(ifname)) then
            ierr = nf90_put_att(fid, NF90_GLOBAL, 'source_shapefile', ifname)
            call nf90check(ierr, 'Error while putting attribute "source_shapefile" in '//trim(fname))
        end if

        call get_command(command_line)
        ierr = nf90_put_att(fid, NF90_GLOBAL, 'history', utc_time_string()//': '//trim(command_line))
        call nf90check(ierr, 'Error while putting attribute "history" in '//trim(fname))


        ! Define dimensions and variables
        ! -------------------------------

        ! Dimensions
        ierr = nf90_def_dim(fid, name='longitude', len=nx, dimid=dimid(1))
        call nf90check(ierr, 'Error while defining dimension "longitude" in '//trim(fname))
        ierr = nf90_def_dim(fid, name='latitude', len=ny, dimid=dimid(2))
        call nf90check(ierr, 'Error while defining dimension "latitude" in '//trim(fname))
        ierr = nf90_def_dim(fid, name='bounds', len=2, dimid=dimid(3))
        call nf90check(ierr, 'Error while defining dimension "latitude" in '//trim(fname))
        if (present(nrec)) then
            ierr = nf90_def_dim(fid, name='records', len=nrec, dimid=dimid(4))
            call nf90check(ierr, 'Error while defining dimension "records" in '//trim(fname))
        end if

        ! Define lon variable
        ierr = nf90_def_var(fid, 'longitude', NF90_REAL, dimid(1:1), loc_varid(1))
        call nf90check(ierr, 'Error while defining variable "longitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(1), 'axis', 'X')
        call nf90check(ierr, 'Error while putting attribute "axis" in variable "longitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(1), 'long_name', 'longitude')
        call nf90check(ierr, 'Error while putting attribute "long_name" in variable "longitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(1), 'standard_name', 'longitude')
        call nf90check(ierr, 'Error while putting attribute "standard_name" in variable "longitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(1), 'units', 'degrees_east')
        call nf90check(ierr, 'Error while putting attribute "units" in variable "longitude" in '//trim(fname))

        ! Define lat variable
        ierr = nf90_def_var(fid, 'latitude', NF90_REAL, dimid(2:2), loc_varid(2))
        call nf90check(ierr, 'Error while defining variable "latitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(2), 'axis', 'Y')
        call nf90check(ierr, 'Error while putting attribute "axis" in variable "latitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(2), 'long_name', 'latitude')
        call nf90check(ierr, 'Error while putting attribute "long_name" in variable "latitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(2), 'standard_name', 'latitude')
        call nf90check(ierr, 'Error while putting attribute "standard_name" in variable "latitude" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(2), 'units', 'degrees_north')
        call nf90check(ierr, 'Error while putting attribute "units" in variable "latitude" in '//trim(fname))

        ! Define bounds variable
        ierr = nf90_def_var(fid, 'bounds', NF90_INT, dimid(3:3), loc_varid(3))
        call nf90check(ierr, 'Error while defining variable "bounds" in '//trim(fname))

        ! Define record variable
        if (present(nrec)) then
            ierr = nf90_def_var(fid, 'records', NF90_INT, dimid(4:4), loc_varid(4))
            call nf90check(ierr, 'Error while defining variable "records" in '//trim(fname))
            ierr = nf90_put_att(fid, loc_varid(4), 'long_name', 'record number in source shapefile')
            call nf90check(ierr, 'Error while putting attribute "long_name" in variable "records" in '//trim(fname))
        end if

        ! Define lon bounds variable
        ierr = nf90_def_var(fid, 'longitude_bounds', NF90_REAL, (/dimid(1),dimid(3)/), loc_varid(5))
        call nf90check(ierr, 'Error while defining variable "longitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(5), 'long_name', 'longitude boundaries')
        call nf90check(ierr, 'Error while putting attribute "long_name" in variable "longitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(5), 'standard_name', 'lon_bounds')
        call nf90check(ierr, 'Error while putting attribute "standard_name" in variable "longitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(5), 'units', 'degrees_east')
        call nf90check(ierr, 'Error while putting attribute "units" in variable "longitude_bounds" in '//trim(fname))

        ! Define lat bounds variable
        ierr = nf90_def_var(fid, 'latitude_bounds', NF90_REAL, dimid(2:3), loc_varid(6))
        call nf90check(ierr, 'Error while defining variable "latitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(6), 'long_name', 'latitude boundaries')
        call nf90check(ierr, 'Error while putting attribute "long_name" in variable "latitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(6), 'standard_name', 'lat_bounds')
        call nf90check(ierr, 'Error while putting attribute "standard_name" in variable "latitude_bounds" in '//trim(fname))
        ierr = nf90_put_att(fid, loc_varid(6), 'units', 'degrees_north')
        call nf90check(ierr, 'Error while putting attribute "units" in variable "latitude_bounds" in ' &
                            //trim(fname))

        ! Define area variable
        if (present(grid_area)) then
            ierr = nf90_def_var(fid, 'grid_area', NF90_REAL, dimid(1:2), loc_varid(7))
            call nf90check(ierr, 'Error while defining variable "grid_area" in '//trim(fname))
            ierr = nf90_put_att(fid, loc_varid(7), 'long_name', 'grid cell area')
            call nf90check(ierr, 'Error while putting attribute "long_name" in variable "area" in '//trim(fname))
            ierr = nf90_put_att(fid, loc_varid(7), 'standard_name', 'area')
            call nf90check(ierr, 'Error while putting attribute "standard_name" in variable "area" in '//trim(fname))
            ierr = nf90_put_att(fid, loc_varid(7), 'units', 'm2')
            call nf90check(ierr, 'Error while putting attribute "units" in variable "area" in '//trim(fname))
        end if

        ! Define polygon fraction variable
        if (fract) then
            if (present(nrec)) then
                ierr = nf90_def_var(fid, 'polygon_fraction', NF90_REAL, (/dimid(1:2),dimid(4)/), varid)
            else
                ierr = nf90_def_var(fid, 'polygon_fraction', NF90_REAL, dimid(1:2), varid)
            end if
            call nf90check(ierr, 'Error while defining variable "polygon_fraction" in '//trim(fname))
            ierr = nf90_put_att(fid, varid, 'long_name', 'fraction of grid cell covered by polygon')
            call nf90check(ierr, 'Error while putting attribute "long_name" in variable "polygon_fraction" in '//trim(fname))
            ierr = nf90_put_att(fid, varid, 'units', '-')
            call nf90check(ierr, 'Error while putting attribute "units" in variable "polygon_fraction" in '//trim(fname))
        else
            if (present(nrec)) then
                ierr = nf90_def_var(fid, 'polygon_area', NF90_REAL, (/dimid(1:2),dimid(4)/), varid)
            else
                ierr = nf90_def_var(fid, 'polygon_area', NF90_REAL, dimid(1:2), varid)
            end if
            call nf90check(ierr, 'Error while defining variable "polygon_area" in '//trim(fname))
            ierr = nf90_put_att(fid, varid, 'long_name', 'area grid cell covered by polygon')
            call nf90check(ierr, 'Error while putting attribute "long_name" in variable "polygon_area" in '//trim(fname))
            ierr = nf90_put_att(fid, varid, 'units', 'm2')
            call nf90check(ierr, 'Error while putting attribute "units" in variable "polygon_area" in '//trim(fname))
        end if

        ! End definition mode:
        ierr = nf90_enddef(fid)
        call nf90check(ierr, 'error while end of definition of file '//trim(fname))


        ! Put variable
        ! ------------

        ! Longitude
        if (flip_x) then
            ierr = nf90_put_var(fid, loc_varid(1), x(nx:1:-1))
        else
            ierr = nf90_put_var(fid, loc_varid(1), x)
        end if
        call nf90check(ierr, 'Error while putting variable "longitude" in '//trim(fname))

        ! Latitude
        if (flip_y) then
            ierr = nf90_put_var(fid, loc_varid(2), y(ny:1:-1))
        else
            ierr = nf90_put_var(fid, loc_varid(2), y)
        end if
        call nf90check(ierr, 'Error while putting variable "latitude" in '//trim(fname))

        ! Bounds
        ierr = nf90_put_var(fid, loc_varid(3), (/1,2/))
        call nf90check(ierr, 'Error while putting variable "bounds" in '//trim(fname))

        ! Records
        if (present(nrec)) then
            ierr = nf90_put_var(fid, loc_varid(4), (/(i,i=1,nrec)/))
            call nf90check(ierr, 'Error while putting variable "record" in '//trim(fname))
        end if

        ! Longitude bounds
        if (flip_x) then
            ierr = nf90_put_var(fid, loc_varid(5), xbounds(nx+1:2:-1),   start=(/1,1/), count=(/nx,1/))
            ierr = nf90_put_var(fid, loc_varid(5), xbounds(nx:1:-1),     start=(/1,2/), count=(/nx,1/))
        else
            ierr = nf90_put_var(fid, loc_varid(5), xbounds(1:nx),   start=(/1,1/), count=(/nx,1/))
            ierr = nf90_put_var(fid, loc_varid(5), xbounds(2:nx+1), start=(/1,2/), count=(/nx,1/))
        end if
        call nf90check(ierr, 'Error while putting variable "longitude_bounds" in '//trim(fname))

        ! Latitude bounds
        if (flip_y) then
            ierr = nf90_put_var(fid, loc_varid(6), ybounds(ny+1:2:-1),   start=(/1,1/), count=(/ny,1/))
            ierr = nf90_put_var(fid, loc_varid(6), ybounds(ny:1:-1),     start=(/1,2/), count=(/ny,1/))
        else
            ierr = nf90_put_var(fid, loc_varid(6), ybounds(1:ny),   start=(/1,1/), count=(/ny,1/))
            ierr = nf90_put_var(fid, loc_varid(6), ybounds(2:ny+1), start=(/1,2/), count=(/ny,1/))
        end if
        call nf90check(ierr, 'Error while putting variable "latitude_bounds" in '//trim(fname))

        ! Area
        if (present(grid_area)) then
            ierr = nf90_put_var(fid, loc_varid(7), grid_area)
            call nf90check(ierr, 'Error while putting variable "area" in '//trim(fname))
        end if



    end subroutine


    !======================================================================!


    subroutine write_output(fid, varid, nx, ny, flip_x, flip_y, i0, outvar)
        integer, intent(in):: fid, varid, nx, ny, i0
        logical, intent(in):: flip_x, flip_y
        double precision, dimension(nx,ny), intent(in):: outvar
        integer:: ierr
        if (flip_x .and. flip_y) then
            ierr = nf90_put_var(fid, varid, outvar(i0-1:1:-1, ny:1:-1), start=(/1,1/),   count=(/i0-1,ny/))
            ierr = nf90_put_var(fid, varid, outvar(nx:i0:-1,  ny:1:-1), start=(/i0,1/),  count=(/nx+1-i0,ny/))
        elseif (flip_x) then
            ierr = nf90_put_var(fid, varid, outvar(i0-1:1:-1, :), start=(/1,1/),   count=(/i0-1,ny/))
            ierr = nf90_put_var(fid, varid, outvar(nx:i0:-1,  :), start=(/i0,1/),  count=(/nx+1-i0,ny/))
        elseif (flip_y) then
            ierr = nf90_put_var(fid, varid, outvar(i0:nx,  ny:1:-1), start=(/1,1/),       count=(/nx+1-i0,ny/))
            ierr = nf90_put_var(fid, varid, outvar(1:i0-1, ny:1:-1), start=(/nx+2-i0,1/), count=(/i0-1,ny/))
        else
            ierr = nf90_put_var(fid, varid, outvar(i0:nx,  :), start=(/1,1/),       count=(/nx+1-i0,ny/))
            ierr = nf90_put_var(fid, varid, outvar(1:i0-1, :), start=(/nx+2-i0,1/), count=(/i0-1,ny/))
        end if
        call nf90check(ierr, 'Error while putting variable')
    end subroutine

    !======================================================================!


    subroutine write_output_slice(fid, varid, nx, ny, flip_x, flip_y, i0, outvar, recnum)
        integer, intent(in):: fid, varid, recnum, nx, ny, i0
        logical, intent(in):: flip_x, flip_y
        double precision, dimension(nx,ny), intent(in):: outvar
        integer:: ierr
        if (flip_x .and. flip_y) then
            ierr = nf90_put_var(fid, varid, outvar(i0-1:1:-1, ny:1:-1), start=(/1,1,recnum/),   count=(/i0-1,ny,1/))
            ierr = nf90_put_var(fid, varid, outvar(nx:i0:-1,  ny:1:-1), start=(/i0,1,recnum/),  count=(/nx+1-i0,ny,1/))
        elseif (flip_x) then
            ierr = nf90_put_var(fid, varid, outvar(i0-1:1:-1, :), start=(/1,1,recnum/),   count=(/i0-1,ny,1/))
            ierr = nf90_put_var(fid, varid, outvar(nx:i0:-1,  :), start=(/i0,1,recnum/),  count=(/nx+1-i0,ny,1/))
        elseif (flip_y) then
            ierr = nf90_put_var(fid, varid, outvar(i0:nx,  ny:1:-1), start=(/1,1,recnum/),       count=(/nx+1-i0,ny,1/))
            ierr = nf90_put_var(fid, varid, outvar(1:i0-1, ny:1:-1), start=(/nx+2-i0,1,recnum/), count=(/i0-1,ny,1/))
        else
            ierr = nf90_put_var(fid, varid, outvar(i0:nx,  :), start=(/1,1,recnum/),       count=(/nx+1-i0,ny,1/))
            ierr = nf90_put_var(fid, varid, outvar(1:i0-1, :), start=(/nx+2-i0,1,recnum/), count=(/i0-1,ny,1/))
        end if
        call nf90check(ierr, 'Error while putting variable slice')
    end subroutine


    !======================================================================!


    subroutine close_file(fid)
        integer, intent(in):: fid
        integer:: ierr
        ierr = nf90_close(fid)
        call nf90check(ierr, 'Error while closing file')
    end subroutine


    !======================================================================!


    subroutine nf90check(ierr, message, nokill)
        integer, intent(in):: ierr
        character(len=*), intent(in):: message
        logical, intent(in), optional:: nokill
        if (ierr/=NF90_NOERR) then
            print *, message
            print *, nf90_strerror(ierr)
            if (present(nokill)) then
                if (.not. nokill) then
                    stop 2
                end if
            else
                stop 2
            end if
        end if
    end subroutine


end module
