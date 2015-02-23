subroutine do_micro

  implicit none

  real(dp) :: logT, logTm, Tp, Tm, logTp, dT, dlogT
  real(dp) :: f, df, lhs, rhs, dlhs, drhs, rerr
  integer :: i, j

  character(len=256) :: outfile
  integer :: funit = 13

  integer, parameter :: maxiter = 100
  real(dp), parameter :: rtol = 1d-6

  ! define the composition

  xa(net_iso(ic12)) = 1.0d0
  !xa(net_iso(io16)) = 1.0d0

  call composition_info( &
       g% num_isos, chem_id, xa, xh, xhe, abar, zbar, z2bar, ye,  &
       mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

  outfile = "carbon_burn.dat"
!  outfile = "oxygen_burn.dat"
  open(unit=funit, file=trim(outfile))

  do i = 0, 20 * 10

     log10Rho = 1 + i*0.05
     Rho = exp10_cr(log10Rho)

     Tm = 1d7
     Tp = 1d10
     
     logTm = log_cr(Tm)
     logTp = log_cr(Tp)


     logT = 0.5 * (logTm + logTp)
     T = exp_cr(logT)
     
     do j = 1, maxiter

!        T = exp_cr(logT)
        log10T = log10_cr(T)

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

        rhs = eps_nuc
        drhs = d_eps_nuc_dT

        lhs = loss(ineu)
        dlhs = loss(idneu_dT)

        f = (rhs - lhs)
        df = (drhs - dlhs)
        
        rerr = abs(lhs - rhs) / min(abs(rhs), abs(lhs))

        ! if ((T+dT).gt.Tp) dT = Tp - T
        ! if ((T-dT).lt.Tm) dT = Tm - T

        ! T = T - dT
        
        if (f .lt. 0) then
           Tm = T
        else
           Tp = T
        endif

        T = sqrt(Tm * Tp)

        if (rerr .lt. rtol) exit
        ! write(*,'(ES12.4, I4, 6ES12.4)') log10Rho, j, f/df, Tm, T, Tp, rerr
        
     enddo

     if (j .lt. maxiter) then 
        write(funit,'(5(ES12.4))') log10Rho, log10T, rhs
     endif

  enddo

  close(funit)

end subroutine do_micro
