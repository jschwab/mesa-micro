program mesa_micro

  ! provide convenient access to MESA microphysics

  use crlibm_lib

  use const_def
  use const_lib

  use chem_def
  use chem_lib

  use eos_def
  use eos_lib

  use kap_def
  use kap_lib

  use net_def
  use net_lib

  use neu_def
  use neu_lib

  use rates_def
  use rates_lib

  use interp_1d_def
  use interp_1d_lib

  use micro_support

  implicit none

  ! variables needed for set up
  integer :: eos_handle, kap_handle, net_handle
  type (Net_General_Info), pointer  :: g
  integer, pointer, dimension(:) :: net_iso, chem_id

  ! variables needed for composition
  real(dp) :: X, Z, Y, abar, zbar, z2bar, ye
  real(dp) :: frac, sumx, xh, xhe, mass_correction
  real(dp), dimension(:), pointer :: xa(:), dabar_dx(:), dzbar_dx(:), dmc_dx(:)

  ! variables needed for eos
  real(dp) :: Rho, T, log10Rho, log10T
  real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar

  ! variables needed for neu
  real(dp) :: log10_Tlim = 7.5
  logical :: flags(num_neu_types)  = .true.
  real(dp) :: loss(num_neu_rvs)
  real(dp) :: sources(num_neu_types, num_neu_rvs)

  ! variables needed for net
  !character(len=256) :: net_name = "approx21.net"
  character(len=256) :: net_name = "ecapture.net"

  type (Net_Info), target :: net_info_target
  type (Net_Info), pointer :: net_info_pointer

  integer :: lwork
  real(dp), dimension(:,:), pointer :: eps_nuc_categories
  real(dp), pointer :: work(:), rate_factors(:)
                  
  real(dp) :: weak_rate_factor, eta, d_eta_dlnT, d_eta_dlnRho, eps_neu_total
  real(dp), dimension(:), pointer ::  &
       d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT
  real(dp), pointer :: d_dxdt_dx(:, :)

  real(dp) :: eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT

  real(dp), parameter :: theta_e_for_graboske_et_al = 1 ! for nondegenerate
  real(dp) :: theta_e = theta_e_for_graboske_et_al

  integer :: screening_mode = extended_screening

  real(dp), dimension(:,:), pointer :: rate_screened, rate_raw
  real(dp) :: max_old_rate_div_new_rate
  integer :: num_rates_reduced

  ! for error handling
  integer :: ierr

  ierr = 0

  ! set up clrlibm
  call crlibm_init

  ! set up const
  call const_init('',ierr)
  if (ierr /= 0) then
     write(*,*) 'const_init failed'
     stop 1
  end if

  ! set up chem
  call chem_init('isotopes.data', ierr)
  if (ierr /= 0) then
     write(*,*) 'chem_init failed'
     stop 1
  end if

  ! set up rates
  call rates_init('reactions.list', '', '', ierr)
  if (ierr /= 0) then
     write(*,*) 'rates_init failed'
     stop 1
  end if

  ! set up net
  call net_init(ierr)
  if (ierr /= 0) then
     write(*,*) 'net_init failed'
     stop 1
  end if

  ! set up eos
  call setup_eos(eos_handle)

  ! set up kap
  call setup_kap(kap_handle)

  ! set up net
  call setup_net(net_name, net_handle, g, chem_id, net_iso)

  ! now that we know the net, allocate for other things
  allocate(xa(g% num_isos), dabar_dx(g% num_isos), dzbar_dx(g% num_isos), dmc_dx(g% num_isos))

  lwork = net_work_size(net_handle, ierr)
  net_info_pointer => net_info_target

  allocate(work(lwork),  &
       rate_factors(g% num_reactions), &
       eps_nuc_categories(num_rvs, num_categories),  &
       rate_screened(num_rvs, g% num_reactions), & 
       rate_raw(num_rvs, g% num_reactions), &
       stat=ierr)
  if (ierr /= 0) stop 1

  rate_factors(:) = 1
  weak_rate_factor = 1

  allocate( &
       d_eps_nuc_dx(g% num_isos),  &
       dxdt(g% num_isos), d_dxdt_dRho(g% num_isos), d_dxdt_dT(g% num_isos), &
       d_dxdt_dx(g% num_isos, g% num_isos))
  if (ierr /= 0) stop 1

  ! call your special function
  call do_micro

  ! tear down rates
  call rates_shutdown

  ! tear down eos
  call shutdown_eos(eos_handle)

  ! tear down kap
  call shutdown_kap(kap_handle)

  ! tear down kap
  call shutdown_net(net_handle)

contains

  !  include 'sample.f'
  !  include 'degen-CO.f'
  ! include 'mesa3-weak.f'
  ! include 'capture-Ne.f'
  ! include 'mesa3-weak-benchmark.f'
  ! include 'eps-weak.f'
  ! include 'ignition-r6596.f'
  ! include 'runaway.f'
  ! include 'list-gaps.f'
  incude 'lifetime.f'
  
end program mesa_micro
