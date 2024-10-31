module shapefile_io_module

!##########################################################################################################################!
!#                                                                                                                        #!
!#                                          #=================================#                                           #!
!#                                          # NOTE: ORGANIZATION OF shp FILE: #                                           #!
!#                                          #=================================#                                           #!
!#                                                                                                                        #!
!#  Default record header (common to all shape types):                                                                    #!
!#      Byte 0: record number                                         [endian=-1, int32bit]                               #!
!#      Byte 4: record length (in 2-bytes words)                      [endian=-1, int32bit]                               #!
!#      Byte 8: shape type (same as in header) or `0' if null record  [endian=+1, int32bit]                               #!
!#                                                                                                                        #!
!#  Polygon (shape type = 5) specific record header                                                                       #!
!#      Byte 8+4:  xmin             [endian=+1, float64bit]                                                               #!
!#      Byte 8+12: ymin             [endian=+1, float64bit]                                                               #!
!#      Byte 8+20: xmax             [endian=+1, float64bit]                                                               #!
!#      Byte 8+28: ymax             [endian=+1, float64bit]                                                               #!
!#      Byte 8+36: number of parts  [endian=+1, int32bit]                                                                 #!
!#      Byte 8+40: number of points [endian=+1, int32bit]                                                                 #!
!#  Polygon record:                                                                                                       #!
!#      Byte 8+44,8+48,...8+4*(nparts-1): [endian=+1, float64bit]                                                         #!
!#                        Position of first points of each part (starting from 0 for the first point of the first part)   #!
!#                                                                                                                        #!
!#      Byte 8+44+4*npart,8+44+4*nparts+16,...,8+44+4*nparts+16*(npoints-1): {x,y} values (2 x float64bit) [endian=+1]    #!
!#                                                                                                                        #!
!##########################################################################################################################!

implicit none

contains


    subroutine check_shapefile_reading(filecode, fileversion, filename)
        integer, intent(in):: filecode, fileversion
        character(len=*), intent(in), optional:: filename

        if (filecode/=9994 .or. fileversion/=1000) then

            if (filecode==170328064 .and. fileversion==-402456576) then
                if (present(filename)) then
                    print *, "Error: reversed endianness in binary file "//trim(filename)
                else
                    print *, "Error: reversed endianness"
                end if
            else
                if (present(filename)) then
                    print *, "Error: can't read properly binary file "//trim(filename)
                else
                    print *, "Error: can't read properly the binary file"
                end if
            end if

            print *, 'Error report:'
            write(*, fmt='(A30,I12)') 'File code read (EXPECT 9994):   ', filecode
            write(*, fmt='(A30,I12)') 'File version read (EXPECT 1000):', fileversion

            stop

        end if

    end subroutine


    !================================================!


    subroutine inquire_shapefile(filename, nrec, rec_headlen, rec_npoints, rec_nparts, &
                                 shptype, npointstot, npointsmax, npartstot, npartsmax, bbox, Zbnds, Mbnds)

        character(len=*), intent(in):: filename
        integer, intent(out), optional:: nrec
        integer, dimension(:), allocatable, intent(out), optional:: rec_headlen, rec_npoints, rec_nparts
        integer, intent(out), optional:: shptype,     npointstot,     npointsmax,     npartstot,     npartsmax
        integer::                    loc_shptype, loc_npointstot, loc_npointsmax, loc_npartstot, loc_npartsmax
        double precision, dimension(4), intent(out), optional:: bbox
        double precision, dimension(4)::                    loc_bbox
        double precision, dimension(2), intent(out), optional:: Zbnds,     Mbnds
        double precision, dimension(2)::                    loc_Zbnds, loc_Mbnds
        logical, dimension(:), allocatable:: is_nullshape
        integer, dimension(:,:), allocatable:: shxdata
        integer, dimension(:), allocatable:: recpos, reclen
        integer:: filecode, fileversion, shpflen, shxflen, npointspos, npartspos, n
        integer, dimension(25):: dummy32bit
        character(len=15), dimension(34):: shape_types
        integer,           dimension(34):: shape_headlen, shape_npointspos, shape_npartspos

        integer, parameter:: NTASK = 132 ! subdivision in N task, print on screen when each task in completed
        integer, dimension(NTASK+1):: list_task
        integer:: itask
        character(len=6):: barformat
        character(len=NTASK) barprogress


        ! Display of progress bar:
        write(barformat, fmt='(A2,I3.3,A1)') '(A',NTASK,')'
        barprogress = ''
        barprogress(1:4) = '| 0%'
        barprogress(NTASK-5:NTASK) = '100% |'


        ! Shape types:
        shape_types(:)  = 'Unknown'
        shape_types(1)  = 'Null Shape'
        shape_types(2)  = 'Point'
        shape_types(4)  = 'PolyLine'
        shape_types(6)  = 'Polygon'
        shape_types(9)  = 'MultiPoint'
        shape_types(12) = 'PointZ'
        shape_types(14) = 'PolyLineZ'
        shape_types(16) = 'PolygonZ'
        shape_types(19) = 'MultiPointZ'
        shape_types(22) = 'PointM'
        shape_types(24) = 'PolyLineM'
        shape_types(26) = 'PolygonM'
        shape_types(29) = 'MultiPointM'
        shape_types(32) = 'MultiPatch'

        ! Header length for in each shape type (in number of 4-bytes (32-bit) words)
        shape_headlen(:)  = 0  !Unknown
        shape_headlen(1)  = 3  !Null Shape
        shape_headlen(2)  = 3  !Point
        shape_headlen(4)  = 13 !PolyLine
        shape_headlen(6)  = 13 !Polygon
        shape_headlen(9)  = 12 !MultiPoint
        shape_headlen(12) = 3  !PointZ
        shape_headlen(14) = 13 !PolyLineZ
        shape_headlen(16) = 13 !PolygonZ
        shape_headlen(19) = 12 !MultiPointZ
        shape_headlen(22) = 3  !PointM
        shape_headlen(24) = 13 !PolyLineM
        shape_headlen(26) = 13 !PolygonM
        shape_headlen(29) = 12 !MultiPointM
        shape_headlen(32) = 13 !MultiPatch

        ! NOTE: the "header" of a record, as considered here, is the length-invariant part encountered at the beginning of each
        !       record of the corresponding shape type. In other words, the number of 4-bytes words before the first 'point' or
        !       'part' information.
        !       It is NOT the shp header stricto sensus, it includes:
        !       the shp header (8bytes) + the shapetype (4bytes) + potential bbox (32bytes), npoints (4bytes) and nparts (4bytes)

        ! Position of 'number of points' information in record "header: for each shape type (in number of 4-bytes (32-bit) words)
        shape_npointspos(:)  = 0  !Unknown shape types or shape types that don't have 'number of points'
        shape_npointspos(4)  = 13 !PolyLine
        shape_npointspos(6)  = 13 !Polygon
        shape_npointspos(9)  = 12 !MultiPoint
        shape_npointspos(14) = 13 !PolyLineZ
        shape_npointspos(16) = 13 !PolygonZ
        shape_npointspos(19) = 12 !MultiPointZ
        shape_npointspos(24) = 13 !PolyLineM
        shape_npointspos(26) = 13 !PolygonM
        shape_npointspos(29) = 12 !MultiPointM
        shape_npointspos(32) = 13 !MultiPatch

        ! Position of 'number of parts' information in record "header: for each shape type (in number of 4-bytes (32-bit) words)
        shape_npartspos(:)  = 0  !Unknown shape types or shape types that don't have 'number of parts'
        shape_npartspos(4)  = 12 !PolyLine
        shape_npartspos(6)  = 12 !Polygon
        shape_npartspos(14) = 12 !PolyLineZ
        shape_npartspos(16) = 12 !PolygonZ
        shape_npartspos(24) = 12 !PolyLineM
        shape_npartspos(26) = 12 !PolygonM
        shape_npartspos(32) = 12 !MultiPatch


        print *
        print *, '====================='
        print *, 'Shapefile pre-reading'
        print *, '====================='
        print *


        ! Read shp and shp files header
        ! -----------------------------

        print *, 'Reading files header...'
        print *

        ! shx file code and length (Big Endian)
        open(unit=10, file=trim(filename)//'.shx', &
             form='unformatted', access='direct', recl=4, status='old', action='read', convert='big_endian')
        read(unit=10, rec=1) filecode
        read(unit=10, rec=7) shxflen
        close(unit=10)

        ! shx file version (Little Endian)
        open(unit=10, file=trim(filename)//'.shx', &
             form='unformatted', access='direct', recl=4, status='old', action='read', convert='little_endian')
        read(unit=10, rec=8) fileversion
        close(unit=10)

        call check_shapefile_reading(filecode, fileversion, trim(filename)//'.shx')

        ! shp file code and length (Big Endian)
        open(unit=11, file=trim(filename)//'.shp', &
             form='unformatted', access='direct', recl=4, status='old', action='read', convert='big_endian')
        read(unit=11, rec=1) filecode
        read(unit=11, rec=7) shpflen
        close(unit=11)

        ! shx file version (Little Endian)
        open(unit=11, file=trim(filename)//'.shx', &
             form='unformatted', access='direct', recl=4, status='old', action='read', convert='little_endian')
        read(unit=11, rec=8) fileversion
        close(unit=11)

        call check_shapefile_reading(filecode, fileversion, trim(filename)//'.shp')

        ! read rest of shx header (stream access, little-endian mode):
        open(unit=10, file=trim(filename)//'.shx', &
             form='unformatted', access='stream', status='old', action='read', convert='little_endian')
        read(unit=10) dummy32bit(1:8), loc_shptype, loc_bbox, loc_Zbnds, loc_Mbnds
        close(unit=10)


        if (present(nrec)) then


            ! convert file length from number of '16-bit words' to number of records:
            !     1) substract header's 50 '16-bit words'
            !     2) divide by 4, as each shx record span on 8 bytes (ie: 4 '16-bit words')
            !+++++++++++++++++++!
            nrec = (shxflen-50)/4
            !+++++++++++++++++++!

            if  (present(rec_headlen) .and. present(rec_npoints) .and. present(rec_nparts)) then


                ! Shape type dependent values
                npointspos = shape_npointspos(loc_shptype+1)
                npartspos  = shape_npartspos(loc_shptype+1)


                ! Allocate record-dependent variables
                ! -----------------------------------

                print *, '# of record:', nrec
                print *

                allocate(    shxdata(2,nrec) )
                allocate(       recpos(nrec) )
                allocate(       reclen(nrec) )
                allocate( is_nullshape(nrec) )
                allocate(  rec_headlen(nrec) )
                allocate(  rec_npoints(nrec) )
                allocate(   rec_nparts(nrec) )


                ! Read shx information
                ! --------------------

                print *, 'Reading shx file...'
                print *

                ! open shx file in stream access, big-endian mode
                open(unit=10, file=trim(filename)//'.shx', &
                    form='unformatted', access='stream', status='old', action='read', convert='big_endian')

                ! read shx header (100 bytes -> 25 '32-bit words')
                read(unit=10) dummy32bit(1:25)

                ! read shx records
                ! - - - - - - - - - !
                read(unit=10) shxdata
                ! - - - - - - - - - !

                ! close file
                close(unit=10)

                ! convert records position and length from '16-bit words' to '4-bytes words' (ie: 32-bit)
                recpos = shxdata(1,:) / 2
                reclen = shxdata(2,:) / 2

                is_nullshape = (reclen == 1) ! Null shape records have a total content length of 1 because only contains ShapeType info


                ! Read shp information (number of points and parts of each record)
                ! ----------------------------------------------------------------

                loc_npointstot = 0
                loc_npointsmax = 0
                loc_npartstot = 0
                loc_npartsmax = 0

                if (npointspos>0) then ! If current shape type have several points

                    ! Subdivision in NTASK groups:
                    list_task = (/ ((n*nrec)/NTASK, n = 0,NTASK) /)
                    ! make sure not to miss any record:
                    list_task(1) = 0
                    list_task(NTASK+1) = nrec

                    print *, 'Reading shp records header'
                    print *, 'Progress:'
                    write(*, fmt=barformat) barprogress

                    open(unit=11, file=trim(filename)//'.shp', &
                         form='unformatted', access='direct', recl=4, status='old', action='read', convert='little_endian')

                    if (npartspos>0) then ! If current shape type also have several parts

                        do itask = 1,NTASK
                            do n = list_task(itask)+1,list_task(itask+1)

                                if (is_nullshape(n)) then
                                    rec_headlen(n) = shape_headlen(1)
                                    rec_npoints(n) = 0
                                    rec_nparts(n) = 0
                                else
                                    rec_headlen(n) = shape_headlen(loc_shptype+1)
                                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
                                    read(unit=11, rec=recpos(n)+npointspos) rec_npoints(n)
                                    read(unit=11, rec=recpos(n)+npartspos) rec_nparts(n)
                                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
                                    loc_npointstot = loc_npointstot + rec_npoints(n)
                                    loc_npartstot  = loc_npartstot  + rec_nparts(n)
                                    if (rec_npoints(n) > loc_npointsmax)  loc_npointsmax = rec_npoints(n)
                                    if (rec_nparts(n)  > loc_npartsmax)   loc_npartsmax  = rec_nparts(n)
                                end if

                            end do
                            write(*, fmt='(A1)', advance='no') '|'
                        end do
                        print *

                    else ! If current shape type doesn't have parts

                        do itask = 1,NTASK
                            do n = list_task(itask)+1,list_task(itask+1)

                                rec_nparts(n) = 0

                                if (is_nullshape(n)) then
                                    rec_headlen(n) = shape_headlen(1)
                                    rec_npoints(n) = 0
                                else
                                    rec_headlen(n) = shape_headlen(loc_shptype+1)
                                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
                                    read(unit=11, rec=recpos(n)+npointspos) rec_npoints(n)
                                    ! - - - - - - - - - - - - - - - - - - - - - - - - - - !
                                    loc_npointstot = loc_npointstot + rec_npoints(n)
                                    if (rec_npoints(n) > loc_npointsmax)  loc_npointsmax = rec_npoints(n)
                                end if

                            end do
                            write(*, fmt='(A1)', advance='no') '|'
                        end do
                        print *

                    end if

                    close(unit=11)

                else ! If current shape type have only 1 point (and no parts):

                    do n = 1,nrec
                        if (is_nullshape(n)) then
                            rec_headlen(n) = shape_headlen(1)
                            rec_npoints(n) = 0
                        else
                            rec_headlen(n) = shape_headlen(loc_shptype+1)
                            rec_npoints(n) = 1
                            loc_npointstot = loc_npointstot + 1
                        end if
                    end do
                    if (loc_npointstot > 0)  loc_npointsmax = 1

                    rec_nparts = 0

                end if

                ! Save other optional variables dependent on "rec_headlen", "rec_npoints" and "rec_nparts":
                ! -----------------------------------------------------------------------------------------

                if (present(npointstot))  npointstot = loc_npointstot
                if (present(npointsmax))  npointsmax = loc_npointsmax
                if (present(npartstot))   npartstot  = loc_npartstot
                if (present(npartsmax))   npartsmax  = loc_npartsmax


            else

                if (present(rec_headlen) .or. present(rec_npoints) .or. present(rec_nparts)) then
                    print *, 'WARNING: improper call signature of subroutine "inquire_shapefile":'
                    print *, 'Cannot assign values to output optional arguments "rec_headlen", "rec_npoints"'
                    print *, 'and "rec_nparts" if all 3 of them are not in the the subroutine call.'
                end if

                if (present(npointstot) .or. present(npointsmax) .or. present(npartstot) .or. present(npartsmax)) then
                    print *, 'WARNING: improper call signature of subroutine "inquire_shapefile":'
                    print *, 'Cannot assign values to output optional arguments "npointstot", "npointsmax",'
                    print *, '"npartstot" or "npartsmax" if the 3 optional arguments "rec_headlen",'
                    print *, '"rec_npoints" and "rec_nparts" are not in the the subroutine call.'
                end if

            end if

        else

            if (present(rec_headlen) .or. present(rec_npoints) .or. present(rec_nparts) .or. present(npointstot) .or. &
                present(npointsmax) .or. present(npartstot) .or. present(npartsmax)) then
                print *, 'WARNING: improper call signature of subroutine "inquire_shapefile":'
                print *, 'Cannot assign values to output optional arguments "npointstot", "npointsmax",'
                print *, '"npartstot", "npartsmax", "rec_headlen", "rec_npoints" or "rec_nparts" if'
                print *, 'optional argument "nrec" is not in the the subroutine call.'
            end if

        end if


        ! Save other independant optional variables:
        ! ------------------------------------------

        if (present(shptype))     shptype    = loc_shptype
        if (present(bbox))        bbox       = loc_bbox
        if (present(Zbnds))       Zbnds      = loc_Zbnds
        if (present(Mbnds))       Mbnds      = loc_Mbnds



    end subroutine



    !==================================================!



    subroutine read_shapefile(filename, nrec, rec_npoints, rec_nparts, rec_pos, rec_pos_in_partlist, part_pos, part_len, points, &
                              bbox)!, rec_bbox)
        ! Load all the records of a shapefile type 5 (polygon) into the following variables:
        ! nrec                      [INTEGER]: number of records
        ! rec_npoints(nrec)         [INTEGER]: number of points for each records
        ! rec_nparts(nrec)          [INTEGER]: number of parts for each records
        ! rec_pos(nrec)             [INTEGER]: positions of the record in the array of points
        ! rec_pos_in_partlist(nrec) [INTEGER]: positions of the record in the arrays of parts ('part_pos' and 'part_len')
        ! part_pos(npartstot)       [INTEGER]: positions of each part of each record (all records in a row) in the array of points
        ! part_len(npartstot)       [INTEGER]: number of points of each part of each record (all records in a row)
        ! points(2,npointstot)      [DOUBLE]:  array of points: {x,y} value of each points (all parts and records in a row)
        ! bbox(4)              [OPT, DOUBLE]:  shapefile bounding-box
        ! rec_bbox(4,nrec)     [OPT, DOUBLE]:  bounding box of each record

        ! Parameters:
        integer, parameter:: EXPECTED_SHAPETYPE = 5 ! POLYGON

        ! I/O variables:
        character(len=*), intent(in):: filename
        integer, intent(out):: nrec
        integer, dimension(:), allocatable, intent(out):: rec_npoints, rec_nparts, &
                                                          rec_pos, rec_pos_in_partlist, part_pos, part_len
        double precision, dimension(:,:), allocatable, intent(out):: points
        ! Optinal I/O variables:
        double precision, dimension(4), intent(out), optional:: bbox
        !double precision, dimension(:,:), allocatable, intent(out), optional:: rec_bbox

        ! Other variables:
        integer, dimension(:), allocatable:: rec_headlen
        integer:: k, kpart, kpoint, n, npartstot, npointstot, shptype
        double precision, dimension(4):: loc_bbox
        !double precision, dimension(:,:), allocatable:: loc_rec_bbox
        character(len=15), dimension(34):: shape_types
        integer, dimension(25):: dummy32bit

        ! * * * * * * * * * * * * * * * *
        integer, parameter:: NTASK = 132 ! subdivision in N task, print on screen when each task in completed
        ! * * * * * * * * * * * * * * * *
        integer, dimension(NTASK+1):: list_task
        integer:: itask
        character(len=6):: barformat
        character(len=NTASK) barprogress


        ! Display of progress bar:
        write(barformat, fmt='(A2,I3.3,A1)') '(A',NTASK,')'
        barprogress = ''
        barprogress(1:4) = '| 0%'
        barprogress(NTASK-5:NTASK) = '100% |'

        ! Shape types:
        shape_types(:)  = 'Unknown'
        shape_types(1)  = 'Null Shape'
        shape_types(2)  = 'Point'
        shape_types(4)  = 'PolyLine'
        shape_types(6)  = 'Polygon'
        shape_types(9)  = 'MultiPoint'
        shape_types(12) = 'PointZ'
        shape_types(14) = 'PolyLineZ'
        shape_types(16) = 'PolygonZ'
        shape_types(19) = 'MultiPointZ'
        shape_types(22) = 'PointM'
        shape_types(24) = 'PolyLineM'
        shape_types(26) = 'PolygonM'
        shape_types(29) = 'MultiPointM'
        shape_types(32) = 'MultiPatch'


        ! Shapefile pre-reading:
        ! ----------------------

        call inquire_shapefile(filename, nrec, rec_headlen, rec_npoints, rec_nparts, &
                               shptype=shptype, npointstot=npointstot, npartstot=npartstot, bbox=loc_bbox)

        if (present(bbox))  bbox = loc_bbox


        ! Check shape type:
        if (shptype/=EXPECTED_SHAPETYPE) then
            print *, 'ERROR: subroutine only handles shape of type', EXPECTED_SHAPETYPE, &
                     '('//trim(shape_types(EXPECTED_SHAPETYPE+1))//')'
            print *, 'Got shape type', shptype, '('//trim(shape_types(shptype+1))//')'
            stop
        end if


        ! Allocate variables:
        ! -------------------

        allocate(             rec_pos(nrec)         )
        allocate( rec_pos_in_partlist(nrec)         )
        allocate(         part_pos(npartstot)    )
        allocate(         part_len(npartstot)    )
        allocate(              points(2,npointstot) )
        !allocate(        loc_rec_bbox(nrec)         )

        ! position of records in array of points
        rec_pos(1) = 1
        rec_pos_in_partlist(1) = 1
        do n = 1,nrec-1
            rec_pos(n+1) = rec_pos(n) + rec_npoints(n)
            rec_pos_in_partlist(n+1) = rec_pos_in_partlist(n) + rec_nparts(n)
        end do


        ! Open shapefile (stream access, little-endian mode):
        ! ---------------------------------------------------

        open(unit=9, file=trim(filename)//'.shp', &
             form='unformatted', access='stream', status='old', action='read', convert='little_endian')


        ! +++++++++++++++++++++++
        ! Read shapefile records:
        ! +++++++++++++++++++++++

        print *
        print *
        print *, '================='
        print *, 'Shapefile reading'
        print *, '================='
        print *


        ! shp header (25 blocs of 4 bytes):
        read(unit=9) dummy32bit(1:25)


        ! Subdivision in NTASK groups:
        list_task = (/ ((k*nrec)/NTASK, k = 0,NTASK) /)
        ! make sure not to miss any record:
        list_task(1) = 0
        list_task(NTASK+1) = nrec

        print *, 'Progress:'
        write(*, fmt=barformat) barprogress

        rec_pos(1) = 1
        do itask = 1,NTASK
            do n = list_task(itask)+1,list_task(itask+1)

                kpart = rec_pos_in_partlist(n)
                kpoint = rec_pos(n)

                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - !
                read(unit=9) dummy32bit(1:rec_headlen(n)), &
                             part_pos(kpart:kpart+rec_nparts(n)-1), &
                             points(:, kpoint:kpoint+rec_npoints(n)-1)
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - !

                part_len(kpart:kpart+rec_nparts(n)-2) =   part_pos(kpart+1:kpart+rec_nparts(n)-1) &
                                                           - part_pos(kpart:kpart+rec_nparts(n)-2)

            end do
            write(*, fmt='(A1)', advance='no') '|'
        end do
        print *


        ! Put length of last part of each records (npoints - pos_of_last_part) if non null shape (npoints > 0):
        do n = 1,nrec
            kpart = rec_pos_in_partlist(n)
            if (rec_npoints(n) > 0)  part_len(kpart+rec_nparts(n)-1) = rec_npoints(n) - part_pos(kpart+rec_nparts(n)-1)
        end do


        ! Close shp file
        ! --------------

        close(unit=9)



    end subroutine



end module

