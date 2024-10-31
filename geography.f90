module geography
implicit none

contains

    subroutine get_earth_parameters(earth_model, radius, flattening, eccentricity, inverse_flattening)

        character(len=*), intent(in):: earth_model
        double precision, intent(out), optional:: radius, flattening, eccentricity, inverse_flattening

        !*******************!
        !*   PARAMETERS:   *!
        !*******************!
        ! Earth ellipsoidal semi-major axis
        double precision, parameter:: WGS84_RADIUS = 6378137.0d0 !m
        ! Earth ellipsoidal excentricity
        double precision, parameter:: WGS84_FM1 = 298.257223563d0 ! inverse flattening
        ! Earth Authalic Radius:
        double precision, parameter:: AUTH_RADIUS = 6371007.2d0 !m

        double precision:: a, flat


        if (earth_model=='WGS84' .or. earth_model=='WGS1984' .or. earth_model=='WGS-1984') then

            flat = 1d0 / WGS84_FM1
            a    = WGS84_RADIUS

        elseif (earth_model=='spherical') then

            flat = 0d0
            a    = AUTH_RADIUS

        else

            print *, 'ERROR: unkown Earth model '//trim(earth_model)
            print *, 'The only implemented Earth model area "WGS84" and "spherical"'
            stop

        end if


        if (present(radius))              radius = a
        if (present(flattening))          flattening = flat
        if (present(eccentricity))        eccentricity = sqrt(1d0 - (1d0-flat)**2)
        if (present(inverse_flattening))  inverse_flattening = 1d0 / flat


    end subroutine



    !=======================================================================!



    function authalic_radius(a, ecc)
        double precision:: authalic_radius
        double precision, intent(in):: a, ecc
        if (ecc>0) then
            authalic_radius    =    a * sqrt( 0.5*(    1    +    ((1-ecc**2) / ecc)  *  log( (1+ecc) / sqrt(1-ecc**2) )   ))
        else
            authalic_radius = a
        end if
    end function


    !----------------------------------------------------------------------!


    function authlat(phi, ecc)
        double precision:: authlat
        double precision, intent(in):: phi, ecc
        double precision:: qp, q
        qp   =   1   -   ( (1 - ecc**2) /(2*ecc) ) * log( (1 - ecc) / (1 + ecc) )
        q    =     (1 - ecc**2)*sin(phi) / (1 - (ecc*sin(phi))**2)   &
               -   ( (1 - ecc**2) /(2*ecc) ) * log( (1 - ecc*sin(phi)) / (1 + ecc*sin(phi)) )
        authlat = asin( q / qp )
    end function


    !----------------------------------------------------------------------!


    subroutine authalic_latitude(phi, ecc)
        ! Compute the authalic latitude (in radians) of the ellipsoidal latitude 'phi' (in RADIANS) for a given Earth eccentricity
        ! 'ecc'
        double precision, dimension(:), intent(inout)::  phi
        double precision, intent(in):: ecc
        double precision:: PI_OVER_2 = acos(-1d0)/2d0
        integer:: i
        if (ecc > 0d0) then
            do i = 1,size(phi,1)
                if (abs(phi(i)) <= (1-1d-15)*PI_OVER_2) then ! if not North or South Pole, at approx. 1d-15
                    phi(i) = authlat(phi(i), ecc)
                end if
            end do
        end if
    end subroutine


end module
