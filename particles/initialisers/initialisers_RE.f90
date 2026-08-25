!> module containing initialization and pdf (probability distribution function) generation
!> subroutines specifically useful for runaway electrons
module initialisers_RE
  use mod_particle_types
  use mod_particle_sim
  use mod_rng
  use initialisers_base
  use constants, only: EL_CHG, ATOMIC_MASS_UNIT, SPEED_OF_LIGHT, MASS_ELECTRON, TWOPI
  use mod_pusher_tools, only: get_orthonormals
  use mod_coordinate_transforms, only: vector_cylindrical_to_cartesian
  use mod_model_settings
  use phys_module, only: CENTRAL_DENSITY, CENTRAL_MASS, ATOMIC_MASS_UNIT, MU_ZERO
  use equil_info
  implicit none

  !> ---------------------------------------------------------------------
  !> State of the kinetic RE drift-surface equilibrium initialization
  !> ('equilibrium' init_function): per-class data and the common profile
  !> function Nprof, read from the file written by the equilibrium solver
  !> (mod_re_kinetic_equilibrium / re_eq_write_output, format version 1).
  !> The marker spatial density of class s is the stationary guiding-centre
  !> density n_s(R,Z) ~ w_s * Nprof(Ahat_s)/R with the per-class invariant
  !> label Ahat_s = (alpha_s R - psi - A_axis)/(A_edge - A_axis); no further
  !> relaxation correction is needed (Bandaru & Hoelzl, PoP 30, 092508
  !> (2023); Bergstroem et al., PPCF (2025) for the full-f PIC model).
  !> ---------------------------------------------------------------------
  character(len=*), parameter :: req_file_name = 're_equilibrium.dat'
  integer                     :: req_n_class = 0    !< number of RE classes
  integer                     :: req_n_l     = 0    !< Nprof table size
  real*8, allocatable         :: req_cl_ekin(:)     !< class kinetic energy [eV]
  real*8, allocatable         :: req_cl_xi(:)       !< class pitch p_par/p
  real*8, allocatable         :: req_cl_w(:)        !< class weights (sum = 1)
  real*8, allocatable         :: req_cl_alpha(:)    !< gamma m_e v_par/e [Wb/m]
  real*8, allocatable         :: req_cl_A_axis(:)   !< label normalization [Wb]
  real*8, allocatable         :: req_cl_A_edge(:)   !< label normalization [Wb]
  real*8, allocatable         :: req_nprof_l(:)     !< Nprof l grid
  real*8, allocatable         :: req_nprof(:)       !< Nprof values [m^-2]
  real*8                      :: req_edge_taper = 0.2d0 !< edge truncation width (read from the file;
                                                        !< MUST match the equilibrium solver)
  real*8                      :: req_I_RE = 0.d0    !< RE current of the equilibrium [A]; the marker
                                                    !< weights are normalized to carry exactly this
                                                    !< current (num_re is ignored for this init)
  real*8                      :: req_l_beam = 1.d0  !< beam-edge label (from the file; MUST match
                                                    !< the equilibrium solver)
  real*8                      :: req_l_beam_width = 0.1d0 !< beam-edge roll-off width (from the file)
  real*8                      :: req_sup_pdf = 1.d0    !< sup of Nprof/R for rejection normalization
  !> Label parameters of the population currently being position-sampled.
  !> re_eq_marker_pdf is `pure` and runs under OpenMP, so it cannot carry
  !> per-marker state; these are set ONCE before each initialise_particles
  !> call and read-only inside it. The class path sets them from
  !> req_cl_*(s); the continuous path sets them per alpha stratum.
  real*8                      :: req_act_alpha  = 0.d0
  real*8                      :: req_act_A_axis = 0.d0
  real*8                      :: req_act_A_edge = 1.d0

  ! --- continuous-distribution descriptor (format version 2). When
  !     req_dist_fmt is not 'ekin_xi_w' the classes in the file are merely
  !     the QUADRATURE the equilibrium used; the markers are drawn from the
  !     underlying continuous f(p,xi) instead.
  character(len=32)           :: req_dist_fmt = 'ekin_xi_w'
  real*8                      :: req_av_zeff  = 1.d0
  real*8                      :: req_av_ztot  = 1.d0
  real*8                      :: req_av_eec   = 10.d0
  real*8                      :: req_av_lnl   = 15.d0
  real*8                      :: req_psi_axis = 0.d0
  integer                     :: req_n_alpha  = 0
  !> A_axis(alpha), A_edge(alpha) on a uniform alpha grid spanning both
  !> signs, written by the equilibrium. Needed because a marker's alpha is
  !> xi p/e with xi from the FULL pitch distribution, which reaches alpha -> 0
  !> (and, at broad pitch, the opposite sign) -- a range the log-spaced class
  !> table cannot cover.
  integer                     :: req_n_atab = 0
  real*8, allocatable         :: req_atab_alpha(:), req_atab_Aax(:), req_atab_Aed(:)
  !> Tabulated CDF of the momentum spectrum, for inverse-CDF sampling
  !> (see req_build_momentum_cdf for why not rejection).
  integer, parameter          :: REQ_NCDF = 4001
  real*8                      :: req_cdf(REQ_NCDF) = 0.d0
  real*8                      :: req_u_lo = 1.d0, req_u_hi = 1.d0, req_dlnu = 0.d0
  logical                     :: req_cdf_ready = .false.

  contains

! Quick and rough function to sample markers based on RZ-coordinates
pure function RZ_pdf(var) result(p)
  real*8, intent(in)  :: var(2) ! var(1)=j
  real*8              :: p
  real*8              :: minor_r
  real*8              :: R_ax, Z_ax

  R_ax = ES%R_axis
  Z_ax = ES%Z_axis

  minor_r = sqrt((var(1)-Z_ax)**2 + (var(2)-R_ax)**2)

  p = 1 / (1 + minor_r)**2
end function RZ_pdf

pure function analytical_pdf(var) result(p)
  real*8, intent(in)  :: var(2) ! var(1)=j
  real*8              :: p
  real*8              :: minor_r
  real*8              :: nu
  real*8              :: R_ax, Z_ax
  real*8              :: LCFS_a

  R_ax   = ES%R_axis
  Z_ax   = ES%Z_axis
  LCFS_a = ES%LCFS_a
  nu     = 2.d0

  minor_r = sqrt((var(1)-Z_ax)**2 + (var(2)-R_ax)**2)

  p = (1.d0 - (minor_r/LCFS_a)**2)**nu

end function analytical_pdf

! Quick and rough function to sample markers proportionally to toroidal current density
pure function current_pdf(var) result(p)
  real*8, intent(in)  :: var(2) ! var(2)=j, var(1) = R
  !real*8, intent(in)  :: var(1) ! var(1)=j
  real*8              :: p
  real*8              :: jzmin, jzmax 

  !> temporarily hard coded, but should be able to be obtained from fluid restart file
  jzmax = 3.0 / 10.0 !1.173 / 10
  jzmin = 0.0001239 / 11.0 !0.0003166 / 11

  p = (var(2)/var(1)-jzmin)/(jzmax-jzmin)

end function current_pdf
    
subroutine basic_initialization(sim, group_num, rng, init_pdf, energy, pitch, std_energy)
  use phys_module, only: tstep_particles
  use mod_kinetic_relativistic
  use mod_sampling, only: boxmueller_transform

  type(particle_sim),                   intent(inout) :: sim
  integer,                              intent(in)    :: group_num
  class(type_rng),                      intent(in)    :: rng
  character(len=50),                    intent(in)    :: init_pdf
  real*8,                               intent(in)    :: energy, pitch ! Kinetic energy in units of eV and pitch
  real*8,               optional,       intent(in)    :: std_energy
  real*8,               allocatable                   :: p_tot(:), p_par(:), p_perp(:)
  integer                                             :: j
  real(kind=8)                                        :: psi, U, gyro_angle
  real(kind=8),         dimension(3)                  :: E, B, B_cart, B_norm
  real*8                                              :: e1(3), e2(3) 
  integer                                             :: num_part
  real*8,               allocatable                   :: ran_uniform(:), ran_gaussian(:)

  select case (trim(init_pdf))
    case ("RZ")
      call initialise_particles(sim%groups(group_num)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-2,-1], transform=RZ_pdf)
    case ("current")
      call initialise_particles(sim%groups(group_num)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-1,var_zj], transform=current_pdf) 
    case ("analytical")
      call initialise_particles(sim%groups(group_num)%particles, sim%fields%node_list, sim%fields%element_list, rng, variables=[-2,-1], transform=analytical_pdf)
    case default
      if (sim%my_id == 0) then
        write(*,*) "ERROR: ", trim(init_pdf), " is not a valid pdf/transform function for "
        write(*,*) "  for group '", sim%groups(group_num)%id, "' when using the 'basic_initialization' function" 
        endif
      stop 1
  end select

  num_part = size(sim%groups(group_num)%particles,1) 

  allocate(p_tot(num_part))
  allocate(p_par(num_part))
  allocate(p_perp(num_part))

  ! Generate gassian distributed energy
  if (present(std_energy)) then

    allocate(ran_uniform(num_part + mod(num_part,2)))
    allocate(ran_gaussian(num_part + mod(num_part,2)))

    call random_number(ran_uniform)
    ran_gaussian = boxmueller_transform(ran_uniform)

    p_tot               = sqrt(((energy + std_energy*ran_gaussian(:num_part))*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

    deallocate(ran_uniform)
    deallocate(ran_gaussian)

  else

    p_tot               = sqrt((energy*EL_CHG/SPEED_OF_LIGHT + MASS_ELECTRON*SPEED_OF_LIGHT)**2 - (MASS_ELECTRON*SPEED_OF_LIGHT)**2)/ATOMIC_MASS_UNIT ! [AMU*m/s]

  end if

  p_par               = pitch * p_tot
  p_perp              = sqrt(p_tot**2 - p_par**2)

  ! Set particle momentum
  select type (particles => sim%groups(group_num)%particles)
  type is (particle_kinetic_relativistic)
    !$omp parallel do default(none) &
    !$omp private(E, B, psi, U, B_cart, B_norm, e1, e2, gyro_angle, j) &
    !$omp shared (sim, tstep_particles, p_par, p_perp, num_part)
    do j=1,num_part

      ! Extract magnetic field and convert to cartesian coordinates
      call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)
      B_cart = vector_cylindrical_to_cartesian(particles(j)%x(3),B)

      B_norm = B_cart/norm2(B_cart)

      ! Generate perpendicular component based on sampled gyro angle
      call get_orthonormals(B_norm, e1, e2)
      call random_number(gyro_angle)
      gyro_angle = gyro_angle * TWOPI

      particles(j)%p = p_par(j) * B_norm + p_perp(j)*(e1*cos(gyro_angle) + e2*sin(gyro_angle))

    end do
    !$omp end parallel do 
  end select

  deallocate(p_tot)
  deallocate(p_par)
  deallocate(p_perp)

end subroutine basic_initialization


!> Read the per-class equilibrium data written by the kinetic RE
!> drift-surface equilibrium solver (re_equilibrium.dat, format version 1).
!> Called by every MPI task (the file is small). Keep in sync with
!> re_eq_write_output in models/mod_re_kinetic_equilibrium.f90.
subroutine read_re_equilibrium_file(my_id)
  implicit none
  integer, intent(in) :: my_id
  integer, parameter  :: iunit = 441
  integer             :: ierr, s, k, idum
  real*8              :: rdum, cols(11)
  character(len=512)  :: line
  character(len=32)   :: key

  if (req_n_class .gt. 0) return   ! already read

  open(iunit, file=req_file_name, status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) "ERROR: cannot open '", req_file_name, "' needed by the"
    write(*,*) "       'equilibrium' RE initialization. Run the equilibrium"
    write(*,*) "       phase with re_kinetic_equilibrium=.true. first."
    stop 1
  endif

  ! --- header: keyword lines, '#' comments
  do
    read(iunit,'(A)',iostat=ierr) line
    if (ierr .ne. 0) then
      write(*,*) "ERROR: unexpected end of ", req_file_name
      stop 1
    endif
    if (index(adjustl(line), '#') .eq. 1) cycle
    read(line,*) key
    select case (trim(key))
    case ('n_class'); read(line,*) key, req_n_class
    case ('n_l');     read(line,*) key, req_n_l
    case ('taper');   read(line,*) key, req_edge_taper
    case ('I_RE');    read(line,*) key, req_I_RE
    case ('l_beam');  read(line,*) key, req_l_beam
    case ('l_beam_w'); read(line,*) key, req_l_beam_width
    case ('q_err', 'psi_bnd'); read(line,*) key, rdum
    case ('psi_axis'); read(line,*) key, req_psi_axis
    case ('dist_fmt'); read(line,*) key, req_dist_fmt
    case ('n_alpha');  read(line,*) key, req_n_alpha
    case ('av_zeff');  read(line,*) key, req_av_zeff
    case ('av_ztot');  read(line,*) key, req_av_ztot
    case ('av_eec');   read(line,*) key, req_av_eec
    case ('av_lnl');   read(line,*) key, req_av_lnl
    case ('n_atab');   read(line,*) key, req_n_atab
    case ('R_edge')
      read(line,*) key, rdum
      exit                          ! last header entry
    case default
      write(*,*) "ERROR: unexpected entry '", trim(key), "' in ", req_file_name
      stop 1
    end select
  enddo

  if ((req_n_class .le. 0) .or. (req_n_l .le. 1)) then
    write(*,*) "ERROR: invalid n_class / n_l in ", req_file_name
    stop 1
  endif

  allocate(req_cl_ekin(req_n_class), req_cl_xi(req_n_class), req_cl_w(req_n_class), &
           req_cl_alpha(req_n_class), req_cl_A_axis(req_n_class), req_cl_A_edge(req_n_class))
  allocate(req_nprof_l(req_n_l), req_nprof(req_n_l))

  ! --- class table (skip the comment line preceding it)
  s = 0
  do while (s .lt. req_n_class)
    read(iunit,'(A)',iostat=ierr) line
    if (ierr .ne. 0) then
      write(*,*) "ERROR: unexpected end of class table in ", req_file_name
      stop 1
    endif
    if (index(adjustl(line), '#') .eq. 1) cycle
    s = s + 1
    read(line,*) idum, cols
    req_cl_ekin(s)   = cols(1)
    req_cl_xi(s)     = cols(2)
    req_cl_w(s)      = cols(3)
    req_cl_alpha(s)  = cols(6)
    req_cl_A_axis(s) = cols(7)
    req_cl_A_edge(s) = cols(8)
  enddo

  ! --- Nprof table
  k = 0
  do while (k .lt. req_n_l)
    read(iunit,'(A)',iostat=ierr) line
    if (ierr .ne. 0) then
      write(*,*) "ERROR: unexpected end of Nprof table in ", req_file_name
      stop 1
    endif
    if (index(adjustl(line), '#') .eq. 1) cycle
    k = k + 1
    read(line,*) req_nprof_l(k), req_nprof(k)
  enddo

  ! --- alpha table (continuous formats only; absent for 'ekin_xi_w')
  if (req_n_atab .gt. 0) then
    allocate(req_atab_alpha(req_n_atab), req_atab_Aax(req_n_atab), &
             req_atab_Aed(req_n_atab))
    k = 0
    do while (k .lt. req_n_atab)
      read(iunit,'(A)',iostat=ierr) line
      if (ierr .ne. 0) then
        write(*,*) "ERROR: unexpected end of alpha table in ", req_file_name
        stop 1
      endif
      if (index(adjustl(line), '#') .eq. 1) cycle
      k = k + 1
      read(line,*) req_atab_alpha(k), req_atab_Aax(k), req_atab_Aed(k)
    enddo
  endif
  close(iunit)

  if (trim(req_dist_fmt) .ne. 'ekin_xi_w') then
    if (req_n_atab .le. 1) then
      write(*,*) "ERROR: continuous distribution '", trim(req_dist_fmt), &
                 "' needs the alpha table; re-run the equilibrium with a"
      write(*,*) "       version that writes format 2."
      stop 1
    endif
  endif

  if (my_id .eq. 0) then
    write(*,'(A,I4,A,I5,A)') "  read '"//req_file_name//"': ", req_n_class, &
      ' RE classes, Nprof table with ', req_n_l, ' points'
    if (trim(req_dist_fmt) .ne. 'ekin_xi_w') then
      write(*,'(A,A,A,I4,A)') "    continuous distribution '", &
        trim(req_dist_fmt), "', alpha table with ", req_n_atab, ' points'
      write(*,'(A,F7.3,A,F7.3,A,F8.3,A,F7.3)') '    Z_eff = ', req_av_zeff, &
        ', Z_tot = ', req_av_ztot, ', E/E_c = ', req_av_eec, &
        ', lnLambda = ', req_av_lnl
    endif
  endif

end subroutine read_re_equilibrium_file


!> Build the cumulative distribution of the avalanche MOMENTUM spectrum,
!>   dn/du ~ exp(-gamma/gamma_0),   u = p/(m_e c),  gamma = sqrt(1+u^2),
!> on a uniform grid in ln u, for inverse-CDF sampling.
!>
!> Inverse CDF rather than rejection, deliberately: the spectrum spans several
!> orders of magnitude, so rejection against any simple envelope would accept
!> ~1 in 1e3-1e4. Inverting a tabulated CDF costs one binary search per marker
!> and accepts everything. It is also the more general choice -- a tabulated
!> f(p) from an external code drops straight in, with no envelope to construct.
subroutine req_build_momentum_cdf()
  implicit none
  integer :: k
  real*8  :: u, gam, pdf_prev, pdf_k, g0

  g0 = sqrt(5.d0 + req_av_zeff) * req_av_lnl          ! gamma_0 = c_Z lnLambda
  req_u_lo = 1.d0
  req_u_hi = 6.d0 * g0                                 ! matches the equilibrium
  req_dlnu = (log(req_u_hi) - log(req_u_lo)) / dble(REQ_NCDF - 1)

  req_cdf(1) = 0.d0
  u   = req_u_lo
  gam = sqrt(1.d0 + u*u)
  pdf_prev = exp(-gam/g0) * u        ! d n / d(ln u) = u dn/du
  do k = 2, REQ_NCDF
    u   = exp(log(req_u_lo) + dble(k-1)*req_dlnu)
    gam = sqrt(1.d0 + u*u)
    pdf_k = exp(-gam/g0) * u
    req_cdf(k) = req_cdf(k-1) + 0.5d0*(pdf_prev + pdf_k)*req_dlnu
    pdf_prev = pdf_k
  enddo
  if (req_cdf(REQ_NCDF) .le. 0.d0) then
    write(*,*) 'ERROR: degenerate RE momentum spectrum (zero total weight)'
    stop 1
  endif
  req_cdf = req_cdf / req_cdf(REQ_NCDF)
  req_cdf_ready = .true.

end subroutine req_build_momentum_cdf


!> Invert the tabulated momentum CDF: return u with CDF(u) = t, t in [0,1].
!> Binary search plus linear interpolation in ln u.
pure function req_draw_u(t) result(u)
  implicit none
  real*8, intent(in) :: t
  real*8             :: u, tt, frac
  integer            :: lo, hi, mid

  tt = min(max(t, 0.d0), 1.d0)
  lo = 1;  hi = REQ_NCDF
  do while (hi - lo .gt. 1)
    mid = (lo + hi) / 2
    if (req_cdf(mid) .le. tt) then
      lo = mid
    else
      hi = mid
    endif
  enddo
  if (req_cdf(hi) .gt. req_cdf(lo)) then
    frac = (tt - req_cdf(lo)) / (req_cdf(hi) - req_cdf(lo))
  else
    frac = 0.d0
  endif
  u = exp(log(req_u_lo) + (dble(lo-1) + frac)*req_dlnu)
end function req_draw_u


!> Piecewise-linear evaluation of the common profile function Nprof at the
!> RAW label l with the SAME edge factor as the equilibrium solver
!> (mod_re_kinetic_equilibrium / re_eq_nprof_at): the C1 smoothstep
!> beam-edge envelope (current confined below req_l_beam, roll-off width
!> req_l_beam_width) times the C1 smoothstep wall taper beyond l = 1
!> (width req_edge_taper). All three widths are read from
!> re_equilibrium.dat; the table itself is raw.
pure function req_nprof_eval(l) result(nval)
  implicit none
  real*8, intent(in) :: l
  real*8             :: nval, x, dl, t
  integer            :: k
  x  = min(max(l, 0.d0), 1.d0)
  dl = req_nprof_l(2) - req_nprof_l(1)
  k  = min(int(x/dl) + 1, req_n_l - 1)
  nval = req_nprof(k) + (req_nprof(k+1) - req_nprof(k)) * (x - req_nprof_l(k)) / dl
  if (req_l_beam .lt. 1.d0) then
    t = (l - (req_l_beam - req_l_beam_width)) / max(req_l_beam_width, 1.d-12)
    t = min(max(t, 0.d0), 1.d0)
    nval = nval * (1.d0 - t*t*(3.d0 - 2.d0*t))
  endif
  if (l .gt. 1.d0) then
    t = min((l - 1.d0) / max(req_edge_taper, 1.d-12), 1.d0)
    nval = nval * (1.d0 - t*t*(3.d0 - 2.d0*t))
  endif
end function req_nprof_eval


!> Rejection-sampling density of the population being sampled (req_act_*):
!> the stationary drift-surface density n_s ~ Nprof(Ahat_s)/R, normalized
!> to [0,1] with req_sup_pdf. Expects var(1) = R, var(2) = psi (from
!> initialise_particles with variables=[-1, var_psi]).
pure function re_eq_marker_pdf(var) result(p)
  implicit none
  real*8, intent(in) :: var(2)
  real*8             :: p, lhat
  ! Label parameters come from the req_act_* scalars, not from a class index:
  ! the continuous path samples per ALPHA STRATUM, which has no class.
  lhat = (req_act_alpha*var(1) - var(2) - req_act_A_axis) &
         / (req_act_A_edge - req_act_A_axis)
  p = req_nprof_eval(lhat) / var(1) / req_sup_pdf
end function re_eq_marker_pdf


!> Linear interpolation of the A_axis / A_edge tables at an arbitrary alpha.
!> The grid is UNIFORM (written that way precisely so this is O(1) and exact
!> at alpha = 0), and A_axis(alpha) ~ -psi_axis + alpha R_axis is near-linear,
!> so linear interpolation is well matched to the data. Outside the table the
!> value is clamped -- |alpha| = |xi| p/e <= p/e = alpha_max by construction,
!> so this can only be hit by round-off at the very ends.
pure function req_atab_at(alpha, tab) result(v)
  implicit none
  real*8, intent(in) :: alpha, tab(:)
  real*8             :: v, x, frac
  integer            :: i
  x = (alpha - req_atab_alpha(1)) &
      / (req_atab_alpha(req_n_atab) - req_atab_alpha(1)) * dble(req_n_atab - 1)
  if (x .le. 0.d0) then
    v = tab(1)
  else if (x .ge. dble(req_n_atab - 1)) then
    v = tab(req_n_atab)
  else
    i    = int(x) + 1
    frac = x - dble(i - 1)
    v    = (1.d0 - frac)*tab(i) + frac*tab(i+1)
  endif
end function req_atab_at


!> Mean pitch of the Embreus pitch distribution at Lorentz factor gamma,
!> <xi> = [1 - e^{-2A}(1+2A)]/[A(1-e^{-2A})] - 1  ->  1/A - 1 for large A.
!> Mirrors re_eq_spectrum_classes in the equilibrium module.
pure function req_mean_pitch(gam) result(xi)
  implicit none
  real*8, intent(in) :: gam
  real*8             :: xi, A, e2
  A  = gam * (req_av_eec + 1.d0) / (req_av_ztot + 1.d0)
  e2 = exp(-2.d0*A)
  xi = (1.d0 - e2*(1.d0 + 2.d0*A)) / (A * (1.d0 - e2)) - 1.d0
end function req_mean_pitch


!> Draw the pitch at Lorentz factor gamma by EXACT inverse CDF.
!> The pitch pdf on s = 1 + xi in [0,2] is A e^{-As}/(1 - e^{-2A}), whose CDF
!> inverts in closed form:  s = -ln(1 - t (1 - e^{-2A})) / A,  t ~ U(0,1).
!> No rejection, so the multi-decade dynamic range of f never costs an
!> acceptance rate. At large A, e^{-2A} underflows to 0 and the expression
!> degenerates to -ln(1-t)/A, which is unbounded as t -> 1; s is clamped to
!> [0,2] for that reason (analytically s(1) = 2 exactly).
pure function req_draw_pitch(gam, t) result(xi)
  implicit none
  real*8, intent(in) :: gam, t
  real*8             :: xi, A, e2, s, tt
  A   = gam * (req_av_eec + 1.d0) / (req_av_ztot + 1.d0)
  e2  = exp(-2.d0*A)
  tt  = min(max(t, 0.d0), 1.d0 - 1.d-15)
  s   = -log(1.d0 - tt*(1.d0 - e2)) / A
  s   = min(max(s, 0.d0), 2.d0)
  xi  = s - 1.d0
end function req_draw_pitch


!> Initialize the RE markers of one group from the kinetic RE drift-surface
!> equilibrium: positions sampled per class from the stationary density
!> n_s ~ w_s Nprof(Ahat_s)/R (this IS the stationary density -- no further
!> relaxation correction is needed), momentum magnitude and pitch from the
!> class values, gyro-angle uniform. Classes are assigned deterministically
!> in the GLOBAL marker index space, so per-class marker counts and weight
!> ratios are consistent across MPI tasks without communication.
!>
!> WEIGHTS: the marker weights are normalized such that the toroidal current
!> carried by the markers, I = sum_p w_p q e v_phi,p / (2 pi R_p), equals
!> EXACTLY the RE current I_RE of the equilibrium (read from
!> re_equilibrium.dat). The namelist num_re is IGNORED for this
!> initialization: any mismatch between the marker current and the fluid
!> equilibrium current makes the resistive term eta*(j - j_RE) nonzero and
!> the coupled state decays away from the constructed equilibrium (observed
!> as an inboard drift of the current channel on the resistive time). The
!> implied physical RE count is reported in the log.
subroutine equilibrium_initialization(sim, group_num, rng)
  use mpi
  use phys_module,             only: part_group_configs, type_part_group_config
  use mod_particle_group_id,   only: matching_part_config_indices
  use mod_particle_allocation, only: calc_n_particles_per_mpi_array
  implicit none
  type(particle_sim), intent(inout) :: sim
  integer,            intent(in)    :: group_num
  class(type_rng),    intent(in)    :: rng

  type(type_part_group_config)      :: config
  integer, dimension(:), allocatable :: n_per_mpi
  integer :: n_global, n_local, i_glob_lo, s, j, i_lo, i_hi, ierr
  ! Population = class (discrete) or alpha stratum (continuous). Allocatable:
  ! the old fixed cls_glob_lo(0:1000) also imposed a hard 1000-class cap.
  integer, allocatable :: pop_glob_lo(:)
  integer :: n_pop, n_in_class
  logical :: req_continuous
  real*8  :: t_lo, t_hi, u_mid, gam_mid, xi_mid, u_j, gam_j, xi_j
  real*8  :: r_u, r_xi, xi_wide_max
  real*8, allocatable :: p_par_a(:), p_perp_a(:)
  real*8  :: cum_w, R_min, p_tot, gyro_angle
  real*8  :: psi, U, e1(3), e2(3)
  real*8, dimension(3) :: E_fld, B_fld, B_cart, B_norm
  real*8  :: weight_s
  real*8  :: I_loc, I_unit, w_factor, p_phi_c, gam, me_kg, phi_p

  config = part_group_configs(matching_part_config_indices(group_num))

  call read_re_equilibrium_file(sim%my_id)

  ! (the former hard cap of 1000 classes is gone: the population index array
  !  is allocatable, and the continuous path uses far more strata than that)

  ! --- global index layout: this task holds global indices
  !     i_glob_lo+1 .. i_glob_lo+n_local of n_global markers
  !     (same splitting as allocate_particles_for_sim)
  n_global  = int(sim%groups(group_num)%n_particles)
  n_per_mpi = calc_n_particles_per_mpi_array(n_global, sim%n_mpi)
  n_local   = n_per_mpi(sim%my_id+1)
  i_glob_lo = sum(n_per_mpi(1:sim%my_id))
  if (n_local .ne. size(sim%groups(group_num)%particles,1)) then
    write(*,*) 'ERROR: equilibrium_initialization: inconsistent local marker count'
    stop 1
  endif

  ! ================================================================
  ! Population layout. A "population" is a class (discrete table) or an
  ! ALPHA STRATUM (continuous distribution). Both are contiguous blocks of
  ! the global marker index, so per-population counts and weight ratios are
  ! identical on every MPI task without communication.
  !
  ! Continuous case: NATURAL sampling, g ~ f, with equal weights. Measured
  ! against the three moments the `rep` coupling actually deposits (P_par,
  ! P_perp, j_Phi), natural sampling minimises the WORST-CASE variance, and
  ! leaves j_Phi -- the moment stationarity depends on -- essentially
  ! noise-free, because v_par ~ +-c for every relativistic marker regardless
  ! of energy. An energy-priority tilt (g ~ gamma f) would make j_Phi ~1000x
  ! noisier to buy P_par. See doc/continuous_distribution_plan.md Sec. 3.1.
  !
  ! The strata are equal-PROBABILITY intervals of the momentum CDF, which is
  ! stratified sampling: every stratum gets the same number of markers and
  ! the same total weight, so all markers end up with equal weight.
  ! ================================================================
  req_continuous = (trim(req_dist_fmt) .ne. 'ekin_xi_w')

  if (req_continuous) then
    if (.not. req_cdf_ready) call req_build_momentum_cdf()
    ! ~1000 markers per stratum, bounded. The stratum alpha is used for the
    ! POSITION pdf only (the momentum and pitch of each marker are drawn
    ! individually), so the error is O(d alpha) across a stratum and this
    ! is a pure numerical choice -- deliberately not a namelist knob.
    n_pop = min(max(n_global/1000, 16), 4096)
    if (allocated(pop_glob_lo)) deallocate(pop_glob_lo)
    allocate(pop_glob_lo(0:n_pop))
    do s = 0, n_pop
      pop_glob_lo(s) = nint(dble(s)/dble(n_pop) * dble(n_global))
    enddo
  else
    n_pop = req_n_class
    if (allocated(pop_glob_lo)) deallocate(pop_glob_lo)
    allocate(pop_glob_lo(0:n_pop))
    pop_glob_lo(0) = 0
    cum_w = 0.d0
    do s = 1, req_n_class
      cum_w = cum_w + req_cl_w(s)
      pop_glob_lo(s) = nint(cum_w * dble(n_global))
    enddo
  endif
  pop_glob_lo(n_pop) = n_global

  ! --- rejection-sampling normalization: sup of Nprof/R over the domain
  R_min = minval(sim%fields%node_list%node(1:sim%fields%node_list%n_nodes)%x(1,1,1))

  xi_wide_max = 0.d0

  do s = 1, n_pop

    n_in_class = pop_glob_lo(s) - pop_glob_lo(s-1)
    if (n_in_class .le. 0) then
      if ((sim%my_id .eq. 0) .and. (.not. req_continuous)) &
        write(*,'(A,I4,A,ES10.2,A)') &
        '  WARNING: RE class ', s, ' (weight ', req_cl_w(s), &
        ') received no markers; its current is not represented'
      cycle
    endif

    ! --- label parameters of this population, and its per-marker weight
    if (req_continuous) then
      ! equal-probability stratum [t_lo, t_hi) of the momentum CDF
      t_lo = dble(s-1)/dble(n_pop)
      t_hi = dble(s)  /dble(n_pop)
      u_mid   = req_draw_u(0.5d0*(t_lo + t_hi))
      gam_mid = sqrt(1.d0 + u_mid*u_mid)
      xi_mid  = req_mean_pitch(gam_mid)
      ! alpha = xi p / e, with p = u m_e c  (gamma m v = p exactly)
      req_act_alpha  = xi_mid * u_mid * MASS_ELECTRON * SPEED_OF_LIGHT / EL_CHG
      req_act_A_axis = req_atab_at(req_act_alpha, req_atab_Aax)
      req_act_A_edge = req_atab_at(req_act_alpha, req_atab_Aed)
      weight_s = 1.d0 / dble(n_global)          ! natural sampling: equal weights
    else
      req_act_alpha  = req_cl_alpha(s)
      req_act_A_axis = req_cl_A_axis(s)
      req_act_A_edge = req_cl_A_edge(s)
      ! provisional weight = class fraction per marker (sum over all markers
      ! = 1); the global current-matching rescale follows below
      weight_s = req_cl_w(s) / dble(n_in_class)
    endif

    ! local slice of this population
    i_lo = max(pop_glob_lo(s-1) + 1, i_glob_lo + 1)          - i_glob_lo
    i_hi = min(pop_glob_lo(s),       i_glob_lo + n_local)    - i_glob_lo
    if (i_lo .gt. i_hi) cycle

    ! --- positions: rejection sampling from Nprof(Ahat)/R
    req_sup_pdf = maxval(req_nprof) / R_min
    call initialise_particles(sim%groups(group_num)%particles(i_lo:i_hi),  &
         sim%fields%node_list, sim%fields%element_list, rng,               &
         variables=[-1, var_psi], transform=re_eq_marker_pdf)

    ! --- momentum. Drawn SERIALLY into arrays rather than inside the OpenMP
    !     loop below: the draws must be reproducible and independent of the
    !     thread count, and random_number is not guaranteed thread-safe.
    if (allocated(p_par_a)) deallocate(p_par_a, p_perp_a)
    allocate(p_par_a(i_lo:i_hi), p_perp_a(i_lo:i_hi))
    if (req_continuous) then
      do j = i_lo, i_hi
        call random_number(r_u)
        call random_number(r_xi)
        ! stratified in the CDF: uniform WITHIN this stratum's interval
        u_j   = req_draw_u(t_lo + (t_hi - t_lo)*r_u)
        gam_j = sqrt(1.d0 + u_j*u_j)
        xi_j  = req_draw_pitch(gam_j, r_xi)
        xi_wide_max = max(xi_wide_max, abs(xi_j - xi_mid))
        p_tot = u_j * MASS_ELECTRON * SPEED_OF_LIGHT / ATOMIC_MASS_UNIT  ! [AMU*m/s]
        p_par_a(j)  = xi_j * p_tot
        p_perp_a(j) = sqrt(max(p_tot**2 - p_par_a(j)**2, 0.d0))
      enddo
    else
      p_tot = sqrt((req_cl_ekin(s)*EL_CHG/SPEED_OF_LIGHT &
                    + MASS_ELECTRON*SPEED_OF_LIGHT)**2   &
                   - (MASS_ELECTRON*SPEED_OF_LIGHT)**2) / ATOMIC_MASS_UNIT  ! [AMU*m/s]
      p_par_a(:)  = req_cl_xi(s) * p_tot
      p_perp_a(:) = sqrt(max(p_tot**2 - (req_cl_xi(s)*p_tot)**2, 0.d0))
    endif

    select type (particles => sim%groups(group_num)%particles)
    type is (particle_kinetic_relativistic)
      !$omp parallel do default(none) &
      !$omp private(E_fld, B_fld, psi, U, B_cart, B_norm, e1, e2, gyro_angle, j) &
      !$omp shared (sim, i_lo, i_hi, p_par_a, p_perp_a, weight_s)
      do j = i_lo, i_hi
        call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, &
                                    particles(j)%x(3), E_fld, B_fld, psi, U)
        B_cart = vector_cylindrical_to_cartesian(particles(j)%x(3), B_fld)
        B_norm = B_cart / norm2(B_cart)
        call get_orthonormals(B_norm, e1, e2)
        call random_number(gyro_angle)
        gyro_angle = gyro_angle * TWOPI
        particles(j)%p = p_par_a(j) * B_norm &
                       + p_perp_a(j)*(e1*cos(gyro_angle) + e2*sin(gyro_angle))
        particles(j)%q      = -1
        particles(j)%weight = weight_s
      end do
      !$omp end parallel do
    class default
      write(*,*) "ERROR: the 'equilibrium' RE initialization requires"
      write(*,*) "       type = 'particle_kinetic_relativistic'"
      stop 1
    end select

    if ((sim%my_id .eq. 0) .and. (.not. req_continuous)) then
      write(*,'(A,I4,A,I10,A,ES12.4)') '  RE class ', s,                   &
        ': global markers ', n_in_class, ', E_kin[eV] ', req_cl_ekin(s)
    endif

  enddo ! populations

  if (allocated(p_par_a)) deallocate(p_par_a, p_perp_a)

  if (req_continuous .and. (sim%my_id .eq. 0)) then
    write(*,'(A,I6,A,I10,A)') '  continuous RE sampling: ', n_pop, &
      ' alpha strata, ', n_global, ' markers (natural sampling, equal weights)'
    write(*,'(A,ES11.4,A,ES11.4)') '    momentum u = p/(m c) drawn on [', &
      req_u_lo, ' ,', req_u_hi
    ! The position pdf uses the STRATUM alpha, but a marker's own alpha is
    ! xi p/e with xi from the full pitch distribution. Within a stratum the
    ! momentum spread is O(1/n_pop) and negligible; the PITCH spread is not
    ! bounded by n_pop and is the real limit of this approximation. It is
    ! harmless when the pitch distribution is narrow (large A, i.e. strong
    ! field / low Z), which is the regime this is intended for.
    write(*,'(A,F8.4)') '    max |xi - <xi>| within a stratum: ', xi_wide_max
    if (xi_wide_max .gt. 0.1d0) then
      write(*,'(A)') '  WARNING: the pitch distribution is broad, so a marker''s own'
      write(*,'(A)') '           alpha departs significantly from its stratum alpha,'
      write(*,'(A)') '           which the POSITION sampling uses. Positions are then'
      write(*,'(A)') '           drawn on a drift surface that is not quite the'
      write(*,'(A)') '           marker''s own. Check j(psihat) against the equilibrium.'
    endif
  endif

  ! --- normalize the weights to the equilibrium RE current: compute the
  !     toroidal current carried by the markers at unit total weight,
  !     I_unit = sum_p w_p q e v_phi,p / (2 pi R_p), and rescale all weights
  !     by I_RE / I_unit
  I_loc = 0.d0
  me_kg = sim%groups(group_num)%mass * ATOMIC_MASS_UNIT
  select type (particles => sim%groups(group_num)%particles)
  type is (particle_kinetic_relativistic)
    do j = 1, n_local
      phi_p   = particles(j)%x(3)
      ! toroidal momentum component in JOREK's cylindrical convention:
      ! the cartesian mapping is x = R cos(phi), y = -R sin(phi) (clockwise
      ! phi, left-handed (R,Z,phi) ordering, cf. mod_coordinate_transforms),
      ! so e_phi = (-sin(phi), -cos(phi), 0)
      p_phi_c = -particles(j)%p(1)*sin(phi_p) - particles(j)%p(2)*cos(phi_p)  ! [AMU m/s]
      gam     = sqrt(1.d0 + (norm2(particles(j)%p)*ATOMIC_MASS_UNIT &
                             / (me_kg*SPEED_OF_LIGHT))**2)
      I_loc   = I_loc + particles(j)%weight * dble(particles(j)%q) * EL_CHG   &
                * (p_phi_c*ATOMIC_MASS_UNIT / (gam*me_kg))                    &
                / (TWOPI * particles(j)%x(1))
    enddo
  end select

  call MPI_Allreduce(I_loc, I_unit, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (abs(I_unit) .le. 0.d0) then
    write(*,*) 'ERROR: equilibrium_initialization: zero marker current at unit weight'
    stop 1
  endif
  w_factor = req_I_RE / I_unit
  if (w_factor .le. 0.d0) then
    write(*,*) 'ERROR: equilibrium_initialization: marker current has the opposite'
    write(*,*) '       sign to the equilibrium I_RE -- the pitch signs of the'
    write(*,*) '       distribution table are inconsistent with the equilibrium.'
    stop 1
  endif

  select type (particles => sim%groups(group_num)%particles)
  type is (particle_kinetic_relativistic)
    do j = 1, n_local
      particles(j)%weight = particles(j)%weight * w_factor
    enddo
  end select

  if (sim%my_id .eq. 0) then
    write(*,'(A,ES13.5,A)') '  weights normalized to the equilibrium RE current I_RE = ', &
      req_I_RE, ' A'
    write(*,'(A,ES13.5)')   '  implied physical RE count sum(weights)         = ', w_factor
    if (config%num_re .gt. 0.d0) then
      write(*,'(A,ES10.2,A)') '  NOTE: part_group_configs%num_re (= ', config%num_re, &
        ') is IGNORED by the equilibrium initialization'
    endif
  endif

end subroutine equilibrium_initialization

end module initialisers_RE