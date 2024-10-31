module scan_arguments_mod
implicit none

contains

    subroutine scan_arguments(list_arg, list_arg_name, list_legal_opt, list_expected_n_optarg, list_optarg, got_opt, &
                              inquire_missing_arg)
        !
        ! This subroutine retrieves the arguments passed to the main program and
        ! stores them in the variable [intent(out)] `list_arg'
        ! 
        ! The number of expected arguments is assumed to be the length of the
        ! character array `list_arg'. It is the ONLY mandatory variable. All
        ! others are optional.
        ! 
        ! By default, the subroutine will ask the program user if not enough
        ! arguments were given. To disable this behaviour, specify
        ! `inquire_missing_arg=.false.`
        ! If the variable `list_arg_name' is present, the subroutine will name
        ! the missing expected argument accordingly when asking for them.
        ! Extra arguments are ignored by the subroutine, and it is mentionned in
        ! a warning message.
        ! 
        ! This subroutine also allows to ask for optional arguments, or options.
        ! To do so, enter the arrays `list_legal_opt' [intent(in)]
        ! It contains the list of possible options (eg: /('-g','-f','--all'/)).
        ! So if a command argument passed to the main program matches one of the
        ! possible options, it is considered as an option, and removed from the
        ! list of arguments.
        ! 
        ! Options may require specific arguments. For instance, you may want to
        ! allow the user to entre `-o name_of_file`, or nothing. So the option
        ! `-o' require 1 additional argument immediately following.
        ! This is the role of the variable [intent(in)] `list_expected_n_optarg'
        ! It specified the number of arguments expected for each options.
        ! If this variable is not present, the subroutine assumes 0 needed
        ! argument for each option.
        ! If not, arguments needed to the options will be store in the variable
        ! [intent(out), dimension(:,:)] `list_optarg'. The first dimension of
        ! `list_optarg' correspond to the options (and MUST match the size of
        ! `list_legal_opt' and `list_expected_n_optarg'). The second dimension
        ! is for the needed arguments of each option, so it should be consistent
        ! with max(list_expected_n_optarg).
        ! 
        ! Finally, if the logical variable [intent(out)] `got_opt' is present,
        ! it stores for each option .true./.false. if the option given by the
        ! user to the main program.
        ! 
        ! EXAMPLE OF USE:
        ! Assuming that the variables:
        !     list_legal_opt         = (/ '-o', '--all', '-i', '-f', '--size' /)
        !     list_expected_n_optarg = (/ 1,    0,       0,    0,    3        /)
        ! are passed to the subroutine "scan_arguments" and that 2 arguments are
        ! expected [size(list_arg) == 2],
        ! The shell command:
        !     `./fortran_executable '/usr/local/bin' -o a.out input_file.txt --size 10 12 9 -f`
        ! Will yield the following [intent(out)] variables after calling the
        ! subroutine "scan_arguments":
        !    list_arg    = (/ '/usr/local/bin', 'input_file.txt' /)
        !    got_opt     = (/ .true., .false., .false., .true., .true. /) [if present]
        !    list_optarg = (/  (/ 'a.out', '', '' /),
        !                      (/ '', '', '' /)     ,
        !                      (/ '', '', '' /)     ,
        !                      (/ '', '', '' /)     ,
        !                      (/ '10', '12', '9' /)   /)
        !

        character(len=*), dimension(:),   intent(out)           :: list_arg
        character(len=*), dimension(:),   intent(in),  optional :: list_arg_name
        character(len=*), dimension(:),   intent(in),  optional :: list_legal_opt
        integer,          dimension(:),   intent(in),  optional :: list_expected_n_optarg
        character(len=*), dimension(:,:), intent(out), optional :: list_optarg
        logical,          dimension(:),   intent(out), optional :: got_opt
        logical, intent(in), optional:: inquire_missing_arg

        integer, dimension(:), allocatable:: loc_list_exp_n_optarg
        logical, dimension(:), allocatable:: loc_got_opt
        logical:: loc_inq_miss_arg
        character(len=500):: arg
        integer:: k, i, karg, kopt
        integer:: ntotarg, nopt, nexparg, noptarg


        ntotarg = command_argument_count()
        nexparg = size(list_arg,1)

        ! initialisation
        list_arg(:) = ''
        if (present(list_optarg)) list_optarg(:,:) = ''



        !-------------------------------------!
        ! If user expects optional arguments: !
        !-------------------------------------!
        if (present(list_legal_opt)) then


            nopt = size(list_legal_opt,1)

            allocate(loc_list_exp_n_optarg(nopt))
            if (present(list_expected_n_optarg)) then
                loc_list_exp_n_optarg = list_expected_n_optarg
            else
                loc_list_exp_n_optarg = 0
            end if

            allocate(loc_got_opt(nopt))
            loc_got_opt = .false.
            

            karg = 0
            k = 0

            ! loop on all the arguments passed to the main program (options and
            ! regular arguments)
            do while (k < ntotarg)
 

                k = k + 1


                ! test whether the current argument it is an option:
                call get_command_argument(k, arg)
                kopt = 0
                do i = 1,nopt
                    if ( trim(arg) == trim(list_legal_opt(i)) ) then
                        kopt = i
                    end if
                end do

                ! if it matches one of the possible options:
                if (kopt>0) then
                    if (loc_got_opt(kopt)) then ! if this option was already passed
                        print *, 'Warning: option '//trim(arg)//' was already given.'
                        print *, 'Will be ignored as well as the',loc_list_exp_n_optarg(kopt),&
                                 'expected following optional argument(s)'
                    else
                        loc_got_opt(kopt) = .true.
                        if (k+loc_list_exp_n_optarg(kopt) > ntotarg) then
                            print *, 'ERROR: not enough arguments for the option: "'//trim(list_legal_opt(kopt))//'"'
                            noptarg = ntotarg - k
                        else
                            noptarg = loc_list_exp_n_optarg(kopt)
                        end if
                        if (present(list_optarg)) then
                            do i = 1,noptarg
                                call get_command_argument(k+i, list_optarg(kopt,i))
                            end do
                        end if
                    end if
                    k = k + loc_list_exp_n_optarg(kopt)

                ! if no match found: considered as a regular argument:
                else
                        
                    karg = karg + 1
                    if (karg <= nexparg) then
                        call get_command_argument(k, list_arg(karg))
                    else
                        print *, 'Warning, unexpected argument '//trim(arg)
                        print *, 'Will be ignored'
                    end if

                end if



            end do

            if (present(got_opt)) got_opt = loc_got_opt



        !---------------------------------------------!
        ! If user does not expect optional arguments: !
        !---------------------------------------------!
        else

            ! loop on all arguments passed to the main program
            do karg = 1,ntotarg
                if (karg <= nexparg) then
                    call get_command_argument(karg, list_arg(karg))
                else
                    call get_command_argument(karg, arg)
                    print *, 'Warning, unexpected argument '//trim(arg)
                    print *, 'Will be ignored'
                end if
            end do

        end if



        !====================================!
        ! In cas not enough input arguments: !
        !====================================!

        if (present(inquire_missing_arg)) then
            loc_inq_miss_arg = inquire_missing_arg
        else
            loc_inq_miss_arg = .true.
        end if

        if (loc_inq_miss_arg) then
            if (present(list_arg_name)) then
                do k = karg+1,nexparg
                    print *, 'Enter argument "'//trim(list_arg_name(k))//'":'
                    read *, list_arg(k)
                end do
            else
                do k = karg+1,nexparg
                    print *, 'Enter argument #',k
                    read *, list_arg(k)
                end do
            end if
        end if



    end subroutine

end module
