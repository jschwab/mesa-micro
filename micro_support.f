module micro_support

  implicit none

contains

  subroutine setup_eos(handle)
    ! allocate and load the eos tables
    use eos_def
    use eos_lib
    integer, intent(out) :: handle

    character (len=256) :: eos_file_prefix
    integer :: ierr
    logical, parameter :: use_cache = .true.

    character (len=8) :: eos_names(num_eos_basic_results)
    
    eos_file_prefix = 'mesa'

    call eos_init(eos_file_prefix, '', '', '', use_cache, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in eos_init'
       stop 1
    end if
    call eos_result_names(eos_names)

    handle = alloc_eos_handle(ierr)
    if (ierr /= 0) then
       write(*,*) 'failed trying to allocate eos handle'
       stop 1
    end if

  end subroutine setup_eos


  subroutine shutdown_eos(handle)
    use eos_def
    use eos_lib
    integer, intent(in) :: handle
    call free_eos_handle(handle)
    call eos_shutdown
  end subroutine shutdown_eos


  subroutine setup_kap(handle)

    use kap_def
    use kap_lib

    integer, intent(out) :: handle
    integer :: ierr

    character(len=256) :: kappa_file_prefix, kappa_CO_prefix, &
         kappa_lowT_prefix
    logical :: use_cache, cubic_X_interpolation, &
         cubic_Z_interpolation, include_electron_conduction, &
         use_Zbase_for_Type1, use_Type2_opacities

    real(dp) :: Mstar, Xc, Xn, Xo, Xne, xheavy, &
         xc_base, xn_base, xo_base, xne_base, &
         frac_Type2, lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT
    real(dp) :: Z_init, dXC, dXO, type2_logT_lower_bdy, Zbase
    real(dp) :: kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy
    real(dp) :: kap_Type2_full_off_X = 0d0, kap_Type2_full_on_X = 0d0
    real(dp) :: kap_Type2_full_off_dZ = 0d0, kap_Type2_full_on_dZ = 0d0

    kappa_file_prefix = 'gn93'
    kappa_lowT_prefix = 'lowT_fa05_gn93'

    kappa_CO_prefix = '' ! use default
    kappa_lowT_prefix = '' ! use default
    kappa_blend_logT_upper_bdy = 0 ! use default
    kappa_blend_logT_lower_bdy = 0 ! use default
    use_cache = .false.
    ierr = 0

    call kap_init( &
         kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
         kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
         use_cache, '', '', ierr) 
    if(ierr/=0) stop 'problem in kap_init'

    handle = alloc_kap_handle(ierr) 
    if (ierr /= 0) stop 'problem in alloc_kap_handle'

    !the final step is to set the type of X and Z interpolation, either linear or cubic spline
    !this is achieved by two logical variables, in this example both X and Z use cubic splines
    cubic_X_interpolation = .true.
    cubic_Z_interpolation = .true.
    include_electron_conduction = .true.
    use_Zbase_for_Type1 = .false.
    use_Type2_opacities = .false.
    call kap_set_choices( &
         handle, cubic_X_interpolation, cubic_Z_interpolation, &
         include_electron_conduction, &
         use_Zbase_for_Type1, use_Type2_opacities, &
         kap_Type2_full_off_X, kap_Type2_full_on_X, &
         kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
         ierr)

    ! using type 2 is like
    ! use_Type2_opacities = .true.
    ! call kap_set_choices( &
    !      handle2, cubic_X_interpolation, cubic_Z_interpolation, &
    !      include_electron_conduction, &
    !      use_Zbase_for_Type1, use_Type2_opacities, &
    !      kap_Type2_full_off_X, kap_Type2_full_on_X, &
    !      kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
    !      ierr)

    if (ierr/=0) stop 'problem in kap_set_interpolation_choices'
    !end of initialization and setup

  end subroutine setup_kap


  subroutine shutdown_kap(handle)
    use kap_def
    use kap_lib
    integer, intent(in) :: handle
    call free_kap_handle(handle)
    call kap_shutdown
  end subroutine shutdown_kap


  subroutine setup_net( &
       net_file, handle, which_rates, which_rates_choice, &
       species, chem_id, net_iso, num_reactions, lwork, ierr)
    use net_lib
    use rates_def, only: rates_reaction_id_max
    use rates_lib, only: set_which_rate_1212

    character (len=*), intent(in) :: net_file
    integer, intent(in) :: which_rates_choice
    integer, pointer :: which_rates(:) ! will be allocated
    integer, pointer :: chem_id(:), net_iso(:) ! set, but not allocated
    integer, intent(out) :: handle, species, num_reactions, lwork, ierr

    ierr = 0
    handle = alloc_net_handle(ierr)
    if (ierr /= 0) then
       write(*,*) 'alloc_net_handle failed'
       return
    end if

    call net_start_def(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_start_def failed'
       return
    end if

    write(*,*) 'load ' // trim(net_file)
    call read_net_file(net_file, handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'read_net_file failed ', trim(net_file)
       return
    end if

    call net_finish_def(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_finish_def failed'
       return
    end if

    allocate(which_rates(rates_reaction_id_max))
    which_rates(:) = which_rates_choice

    call set_which_rate_1212(which_rates, 1)

    call net_set_which_rates(handle, which_rates, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_set_which_rate_f17pg failed'
       return
    end if

    call net_setup_tables(handle, '', ierr)
    if (ierr /= 0) then
       write(*,*) 'net_setup_tables failed'
       return
    end if

    species = net_num_isos(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in net_num_isos'
       return
    end if

    call get_chem_id_table_ptr(handle, chem_id, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in get_chem_id_table_ptr'
       return
    end if

    call get_net_iso_table_ptr(handle, net_iso, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in get_net_iso_table_ptr'
       return
    end if

    num_reactions = net_num_reactions(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in net_num_reactions'
       return
    end if

    lwork = net_work_size(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'failed in net_work_size'
       return
    end if

  end subroutine setup_net
      

  subroutine shutdown_net(handle)
    use net_def
    use net_lib
    integer, intent(in) :: handle
    call free_net_handle(handle)
  end subroutine shutdown_net


end module
