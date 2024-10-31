module miscellaneous_functions
implicit none
contains


    subroutine switch_dble(x,y)
        double precision, intent(inout):: x, y
        double precision:: dummy
        dummy = x
        x = y
        y = dummy
    end subroutine


    !==============================================================================================================================!


    function utc_time_string()
        character(len=23):: utc_time_string
        integer, dimension(8):: time
        integer, dimension(12):: month_length
        integer:: curr_month, prev_month

        month_length = (/31,00,31,30,31,30,31,31,30,31,30,31/)

        call date_and_time(values=time)

        ! time(1): The year
        ! time(2): The month
        ! time(3): The day of the month
        ! time(4): Time difference with UTC in minutes
        ! time(5): The hour of the day
        ! time(6): The minutes of the hour
        ! time(7): The seconds of the minute
        ! time(8): The milliseconds of the second

        if (leap_year(time(1))) then
            month_length(2) = 29
        else
            month_length(2) = 28
        end if

        time(6) = time(6)-time(4) ! UTC time, in minutes (as time(4) = difference with UTC (min) and time(6) = minutes)

        do while (time(6) < 0) ! minute/hour
            time(5) = time(5) - 1
            time(6) = time(6) + 60
        end do
        do while (time(6) >= 60 )
            time(5) = time(5) + 1
            time(6) = time(6) - 60
        end do
        !
        do while (time(5) < 0) ! hour/day
            time(3) = time(3) - 1
            time(5) = time(5) + 24
        end do
        do while (time(5) >= 24 )
            time(3) = time(3) + 1
            time(5) = time(5) - 24
        end do
        !
        do while (time(3) <= 0) ! day/month
            prev_month = modulo(time(2)-2, 12) + 1
            time(2) = time(2) - 1
            time(3) = time(3) + month_length(prev_month)
        end do
        curr_month = modulo(time(2)-1, 12) + 1
        do while (time(3) > month_length(curr_month) )
            time(2) = time(2) + 1
            time(3) = time(3) - month_length(curr_month)
            curr_month = modulo(time(2)-1, 12) + 1
        end do
        !
        do while (time(2) <= 0) ! month/year
            time(1) = time(1) - 1
            time(2) = time(2) + 12
        end do
        do while (time(2) > 12 )
            time(1) = time(1) + 1
            time(2) = time(2) - 12
        end do

        utc_time_string = '1992-03-22 20:30:00 GMT'
        write(utc_time_string(1:4),   fmt='(I4.4)') time(1) ! year
        write(utc_time_string(6:7),   fmt='(I2.2)') time(2) ! month
        write(utc_time_string(9:10),  fmt='(I2.2)') time(3) ! day
        write(utc_time_string(12:13), fmt='(I2.2)') time(5) ! hour
        write(utc_time_string(15:16), fmt='(I2.2)') time(6) ! minute
        write(utc_time_string(18:19), fmt='(I2.2)') time(7) ! second

    end function


    !==============================================================================================================================!


    function leap_year(year)
        logical:: leap_year
        integer, intent(in):: year
        leap_year = (  ( (modulo(year,4)==0) .and. (modulo(year,100)/=0) )  .or.  (modulo(year,400)==0)  )
    end function


    !==============================================================================================================================!


    subroutine remove_extension(filename, ext, list_ext)
        character(len=*), intent(inout):: filename
        character(len=*), intent(in), optional:: ext
        character(len=*), dimension(:), optional:: list_ext
        character(len=100):: loc_ext
        integer:: flen, flen2, extlen, next, n
        logical:: go_on

        flen   = len(filename)
        flen2  = len_trim(filename)

        if (present(ext)) then

            extlen = len(ext)

            if (flen>extlen) then
                if (filename(flen2-extlen+1:flen2) == ext) then
                    filename(flen2-extlen+1:flen2) = ' '
                end if
            end if

        elseif (present(list_ext)) then

            next = size(list_ext)
            do n = 1,next
                loc_ext = list_ext(n)
                extlen = len_trim(loc_ext)
                if (flen>extlen) then
                    if (filename(flen2-extlen+1:flen2) == loc_ext(1:extlen)) then
                        filename(flen2-extlen+1:flen2) = ' '
                    end if
                end if
            end do

        else ! if nothing specified, remove any posible extension (first '.' found):

          n = 1
          go_on = .true.
          do while (go_on)
              if (filename(n:n) == '.') then
                  go_on = .false.
              else
                  n = n + 1
              end if
              if (n==flen2+1) go_on = .false.
          end do
          filename(n:flen2) = ' '


        end if

    end subroutine


    !==============================================================================================================================!


    subroutine remove_path(filename)
        character(len=*), intent(inout):: filename
        integer:: flen, k, k0

        flen = len_trim(filename)

        k0 = 1
        do k = 1,flen-1
            if (filename(k:k) == '/') k0 = k+1
        end do

        filename = filename(k0:flen)

    end subroutine


end module
