module grd_output
implicit none

contains

    subroutine write_grd_output()

        if (output_files(ifile) is '.record()') then
            outname = output_dir + def_name + '_' + key + '.grd'
        else
            outname = output_dir + output_files(ifile) + '_' + key + '.grd'
        end if

        outf = open( outname , mode='w' )

        ! write header:
        outf.write('ncols         '+str(nx)+'\n')
        outf.write('nrows         '+str(ny)+'\n')
        outf.write('xllcorner     '+str(xbounds(0))+'\n')
        outf.write('yllcorner     '+str(ybounds(0))+'\n')
        if (dlon==dlat) then
            outf.write('cellsize      '+str(dlon)+'\n')
        else
            outf.write('cellsize      '+str(dlon)+','+str(dlat)+'\n')
        end if
        outf.write('NODATA_value  -1.\n')

        ! write data:
        do j = ny,1,-1
            do i = 1,nx
                outf.write(str(polygon_area(key)(i,j)))
                outf.write(' ')
            end do
            outf.write('\n')
        end do

        outf.close()

    end subroutine

end module
