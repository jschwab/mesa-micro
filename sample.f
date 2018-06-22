subroutine do_micro

  real(dp) :: beta
  
  ! define the composition

  xa = 0
  xa(net_iso(io16)) = 0.5
  xa(net_iso(ine20)) = 0.5

  call composition_info( &
       species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, &
       ye, mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
  ! call composition_info( &
  !      g% num_isos, chem_id, xa, xh, xhe, abar, zbar, z2bar, ye,  &
  !      mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

  ! define the temperature and density

  log10Rho = 9.6
  Rho = exp10_cr(log10Rho)
  log10T = 8.5
  T = exp10_cr(log10T)
    
  ! get the eos results

  call eosDT_get( &
       eos_handle, Z, X, abar, zbar,  &
       species, chem_id, net_iso, xa, &
       Rho, log10_cr(Rho), T, log10T,  &
       res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

  eta = res(i_eta)
  d_eta_dlnRho = d_dlnd(i_eta)
  d_eta_dlnT = d_dlnT(i_eta)

  beta = (0.511) / (8.617e-11 * T)
  
  write(*,*) res(i_chiT), eta * (8.617e-11) * T, (4 / zbar /(eta) )

  ! get some neu results
  
  call neu_get(T, log10T, Rho, log10Rho, abar, zbar, z2bar, log10_Tlim, flags, &
       loss, sources, ierr)

  ! get some net results

  
  ! call net_get(net_handle, .false., net_info_pointer, g% num_isos, g% num_reactions, &
  !      xa, T, log10T, Rho, log10Rho, &
  !      abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
  !      rate_factors, weak_rate_factor, &
  !      std_reaction_Qs, std_reaction_neuQs, .false., .false., &
  !      eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
  !      dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
  !      screening_mode, theta_e, &
  !      eps_nuc_categories, eps_neu_total, &
  !      lwork, work, ierr)

  call net_get( &
       net_handle, .false., net_info_pointer, species, num_reactions, &
       xa, T, log10T, Rho, log10Rho, & 
       abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
       rate_factors, weak_rate_factor, & 
       std_reaction_Qs, std_reaction_neuQs, .false., .false., &
       eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
       dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
       screening_mode, theta_e_for_graboske_et_al, &     
       eps_nuc_categories, eps_neu_total, & 
       lwork, work, ierr)


    ! write(*,*) eps_nuc, loss(ineu)

end subroutine do_micro
