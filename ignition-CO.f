subroutine do_micro

  implicit none
  
  real(dp) :: logTm, logTp, logT0, lhs, rhs, dlhs_dT, drhs_dT, logTnew, T9
  integer :: i, j

  integer, parameter :: max_iter = 32
  
  character(len=256) :: outfile
  integer :: funit = 13

  ! define the composition

  xa = 0
  xa(net_iso(io16)) = 0.5
  xa(net_iso(ic12)) = 0.5

  call composition_info( &
       g% num_isos, chem_id, xa, xh, xhe, abar, zbar, z2bar, ye,  &
       mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
  
  outfile = "tburn-CO-3e7s.dat"
  open(unit=funit, file=trim(outfile))

  do i = 0, 200

     log10Rho = 0 + 0.05 * i
     Rho = exp10_cr(log10Rho)

     logTm = 7.5
     logTp = 9.5
     logT0 = sqrt(logTm * logTp)

     do j = 1, max_iter

        log10T = sqrt(logTm * logTp)
        T = exp10_cr(log10T)

        call eosDT_get( &
             kap_handle, Z, X, abar, zbar,  &
             g% num_isos, chem_id, net_iso, xa, &
             Rho, log10_cr(Rho), T, log10T,  &
             res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)

        eta = res(i_eta)
        d_eta_dlnRho = d_dlnd(i_eta)
        d_eta_dlnT = d_dlnT(i_eta)

        call net_get(net_handle, .false., net_info_pointer, g% num_isos, g% num_reactions, &
             xa, T, log10T, Rho, log10Rho, &
             abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
             rate_factors, weak_rate_factor, &
             std_reaction_Qs, std_reaction_neuQs, .false., .false., &
             eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
             dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
             screening_mode, theta_e, &
             eps_nuc_categories, eps_neu_total, &
             lwork, work, ierr)


        call neu_get(T, log10T, Rho, log10Rho, abar, zbar, z2bar, log10_Tlim, flags, &
             loss, sources, ierr)

        lhs = eps_nuc
        lhs = 3e7
        dlhs_dT = d_eps_nuc_dT

        rhs = loss(ineu)
        rhs = res(i_Cp) * T / eps_nuc 
        drhs_dT = loss(idneu_dT)

        ! just do bisection
        if (lhs .lt. rhs) then 
           logTm = log10T
        else
           logTp = log10T
        endif

     enddo

     write(funit,'(5(ES12.4))') log10Rho, log10T, lhs, rhs
  enddo

  close(funit)

end subroutine do_micro
