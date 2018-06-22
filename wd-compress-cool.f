subroutine do_micro

  implicit none

  real(dp) :: logT, logTm, Tp, Tm, logTp, dT, dlogT
  real(dp) :: f, df, lhs, rhs, dlhs, drhs, rerr
  integer :: i, j

  character(len=256) :: outfile
  integer :: funit = 13

  integer, parameter :: nlines = 805
  real(dp), dimension(nlines)  :: zc, zeta1, phip1, &
       rhoc, mwd, rwd, dlnmdlnd, dlnrdlnd, dlnmdlnr

  real(dp) :: mwd_interp, dlnmdlnd_interp
  real(dp), dimension(:, :), pointer :: f1, f2
  real(dp), dimension(:), pointer :: f1_1, f2_1
  real(dp), dimension(:), pointer :: scratch1a, scratch1b
  real(dp), dimension(:), pointer :: scratch2a, scratch2b
  
  integer :: ierr
  character(len=256) :: dbg
  character(len=256) :: header

  real(dp) :: mdot, tcompress
  
  integer, parameter :: maxiter = 100
  real(dp), parameter :: rtol = 1d-6

  ! read in zero temperature white dwarf data
  open(unit=funit, file="/tmp/ZeroTemperatureWhiteDwarf.dat")
  do i = 1, 23
     read(funit,*) header
  enddo
  do i = 1, nlines
     read(funit, *) zc(i), zeta1(i), phip1(i), &
          rhoc(i), mwd(i), rwd(i), &
          dlnmdlnd(i), dlnrdlnd(i), dlnmdlnr(i)
  enddo
  close(funit)

  ! set up interpolation
  ! want to know Mwd(rhoc) and dlnmdlnr(rhoc)
  allocate(scratch1a(4*nlines), scratch2a(4*nlines))
  allocate(scratch1b(nlines*mp_work_size), scratch2b(nlines*mp_work_size))

  f1_1 => scratch1a
  f1(1:4,1:nlines) => f1_1(1:4*nlines)
  f1(1,1:nlines) = mwd

  f2_1 => scratch2a
  f2(1:4,1:nlines) => f2_1(1:4*nlines)
  f2(1,1:nlines) = dlnmdlnd

  
  call interp_m3a(log10(rhoc), nlines, f1_1, mp_work_size, scratch1b, dbg, ierr)
  call interp_m3a(log10(rhoc), nlines, f2_1, mp_work_size, scratch2b, dbg, ierr)
  

  ! define the composition
  ! xa(net_iso(io16)) = 0.50d0
  ! xa(net_iso(ine20)) = 0.45d0
  ! xa(net_iso(img24)) = 0.05d0

  xa(net_iso(ic12)) = 0.5d0
  xa(net_iso(io16)) = 0.5d0

  call composition_info( &
       g% num_isos, chem_id, xa, xh, xhe, abar, zbar, z2bar, ye,  &
       mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)

  mdot = 3e-8
  outfile = "CO-3em8.dat"
  open(unit=funit, file=trim(outfile))

  do i = 0, 20 * 5

     log10Rho = 6 + i*0.05
     Rho = exp10_cr(log10Rho)

     ! get the value of tcompress
     call interp_value(log10(rhoc), nlines, f1_1, log10Rho, mwd_interp, ierr)
     call interp_value(log10(rhoc), nlines, f2_1, log10Rho, dlnmdlnd_interp, ierr)

     tcompress = (mwd_interp/mdot) * dlnmdlnd_interp * secyer

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

        rhs = tcompress
        drhs = 0

        lhs = res(i_Cp) * T / loss(ineu)
        dlhs = 0

        f = (rhs - lhs)
        df = (drhs - dlhs)
        
        rerr = abs(lhs - rhs) / min(abs(rhs), abs(lhs))

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
        write(funit,'(7(ES11.3))') log10Rho, log10T, &
             res(i_Cp) * T, loss(ineu), &
             mwd_interp, dlnmdlnd_interp, tcompress
     endif

  enddo

  close(funit)

end subroutine do_micro
