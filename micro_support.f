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

    eos_file_prefix = 'mesa'

    call eos_init(eos_file_prefix, '', '', use_cache, ierr)
    if (ierr /= 0) then
       write(*,*) 'eos_init failed in Setup_eos'
       stop 1
    end if

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
         cubic_Z_interpolation, include_electron_conduction

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
    type2_logT_lower_bdy = 0 ! use default
    use_cache = .false.
    ierr = 0

    call kap_init( &
         kappa_file_prefix, kappa_CO_prefix, kappa_lowT_prefix, &
         kappa_blend_logT_upper_bdy, kappa_blend_logT_lower_bdy, &
         type2_logT_lower_bdy, use_cache, '', ierr) 
    if(ierr/=0) stop 'problem in kap_init'


    handle = alloc_kap_handle(ierr) 
    if (ierr /= 0) stop 'problem in alloc_kap_handle'

    !the final step is to set the type of X and Z interpolation, either linear or cubic spline
    !this is achieved by two logical variables, in this example both X and Z use cubic splines
    cubic_X_interpolation = .true.
    cubic_Z_interpolation = .true.
    include_electron_conduction = .true.
    call kap_set_choices( &
         handle, cubic_X_interpolation, cubic_X_interpolation, &
         include_electron_conduction, &
         kap_Type2_full_off_X, kap_Type2_full_on_X, &
         kap_Type2_full_off_dZ, kap_Type2_full_on_dZ, &
         ierr)
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



  subroutine setup_net(net_file_in, handle, g, chem_id, net_iso)

    use const_def, only : dp

    use net_def
    use net_lib
    use chem_def
!    use chem_lib
    use rates_def

    character (len=*), intent(in) :: net_file_in
    integer, intent(out) :: handle
    type (Net_General_Info), pointer, intent(out)  :: g
    integer, dimension(:), pointer, intent(out) :: net_iso, chem_id

    integer :: species, num_reactions    
    integer :: info, i, ierr

    character (len=64) :: net_file
    character (len=256) :: cache_suffix


    integer :: which_rates_choice

    integer, pointer ::  &
         reaction_id(:), reaction_table(:), which_rates(:)

    net_file = net_file_in

    call net_init(ierr)
    if (ierr /= 0) stop 1

    handle = alloc_net_handle(ierr)
    if (ierr /= 0) then
       write(*,*) 'alloc_net_handle failed'
       stop 2
    end if
         
    call net_start_def(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_start_def failed'
       stop 2
    end if
         
    call read_net_file(net_file, handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'read_net_file failed ', trim(net_file)
       stop 2
    end if
         
    call net_finish_def(handle, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_finish_def failed'
       stop 2
    end if
   	
    allocate(reaction_id(rates_reaction_id_max))
    allocate(reaction_table(rates_reaction_id_max))
    allocate(which_rates(rates_reaction_id_max))

    which_rates_choice = rates_NACRE_if_available      
    which_rates(:) = which_rates_choice

    call net_set_which_rates(handle, which_rates, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_set_which_rate_f17pg failed'
       stop 2
    end if

    cache_suffix = ''
    call net_setup_tables(handle, 'rate_tables', cache_suffix, info)
    if (ierr /= 0) then
       write(*,*) 'net_setup_tables failed'
       stop 2
    end if

    call net_ptr(handle, g, ierr)
    if (ierr /= 0) then
       write(*,*) 'net_ptr failed'
       stop 2
    end if
         
    species = g % num_isos
    num_reactions = g% num_reactions

    allocate(net_iso(num_chem_isos))
    allocate(chem_id(num_chem_isos))
         
    call get_chem_id_table(handle, species, chem_id, ierr)
    if (ierr /= 0) then
       write(*,*) 'get_chem_id_table failed'
       stop 2
    end if

    call get_net_iso_table(handle, net_iso, ierr)
    if (ierr /= 0) then
       write(*,*) 'get_net_iso_table failed'
       stop 2
    end if

    call get_reaction_id_table(handle, num_reactions, reaction_id, ierr)
    if (ierr /= 0) then
       write(*,*) 'get_reaction_id_table failed'
       stop 2
    end if

    call get_net_reaction_table(handle, reaction_table, ierr)
    if (ierr /= 0) then
       write(*,*) 'get_net_reaction_table failed'
       stop 2
    end if
         
  end subroutine setup_net
      

  subroutine shutdown_net(handle)
    use net_def
    use net_lib
    integer, intent(in) :: handle
    call free_net_handle(handle)
  end subroutine shutdown_net


end module
