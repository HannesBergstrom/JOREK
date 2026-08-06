!> Kinetic runaway-electron (RE) drift-surface equilibrium.
!>
!> Given a discrete set of RE momentum/pitch classes s = 1..n (read from a
!> plain-text table) and a target safety-factor profile q_t(psihat_n), the
!> Grad-Shafranov source is built such that the RE current density lies on
!> drift surfaces, i.e. surfaces of constant canonical toroidal momentum
!>
!>     A_s(R,Z) = gamma_s m_e v_par,s R - e psi        (small-pitch form)
!>
!> instead of flux surfaces. The stationary guiding-centre density of each
!> class is n_s = w_s Nprof(Ahat_s)/R with ONE common normalized profile
!> function Nprof shared by all classes, where Ahat_s is the per-class
!> normalized invariant label (0 at the class drift axis, 1 at the outboard
!> midplane edge). The resulting toroidal current density
!> j_phi = -e sum_s v_par,s n_s enters the GS equation as
!>
!>     Delta* psi = mu0 e sum_s v_par,s w_s Nprof(Ahat_s(psi;R,Z))
!>
!> (the 1/R of n_s cancels the R of the GS right-hand side). Nprof is
!> determined by an outer profile-transplant iteration matching the target
!> q profile. Physics basis: V. Bandaru and M. Hoelzl, Phys. Plasmas 30,
!> 092508 (2023), Formulation B, Eqs. (3.23)-(3.27) and Sec. IV, generalized
!> from mono-energetic to multiple classes. The motivation is that at high RE
!> energy the drift-orbit shift is comparable to profile scale lengths near
!> the axis, so markers initialized on flux surfaces do not relax to the
!> intended q profile; the only quantity conserved per marker during an
!> axisymmetric relaxation is A_s, so the equilibrium is constructed in
!> invariant space from the start.
!>
!> A standalone reference implementation with the validation of the physics
!> and of the iteration scheme lives in util/re_equilibrium_prototype/
!> (see its README for the numerically-motivated choices: cumulative
!> transplant, ratio clamp, q-space convergence metric).
!>
!> All new functionality is inactive unless re_kinetic_equilibrium = .true.
!> Everything here works in the same units as the JOREK equilibrium solver
!> (psi, R, Z as in the restart files; mu0 explicit in the source), with the
!> class quantities (v_par, Nprof) kept in SI.
module mod_re_kinetic_equilibrium

use constants, only: MU_ZERO, EL_CHG, MASS_ELECTRON, SPEED_OF_LIGHT, PI

implicit none

private

! --- namelist parameters
public :: re_kinetic_equilibrium, re_eq_dist_file, re_eq_dist_format,          &
          re_eq_q_file, re_eq_match_mode, re_eq_transplant, re_eq_map_mode,    &
          re_eq_I_RE,                                                          &
          re_eq_xi_min, re_eq_alpha_out, re_eq_tol_q, re_eq_tol_q_soft,        &
          re_eq_edge_taper, re_eq_l_beam, re_eq_l_beam_width,                  &
          re_eq_ratio_clamp, re_eq_absorbing_edge, re_eq_op_lambda,             &
          re_eq_alpha_current,                                                 &
          re_eq_max_it_out, re_eq_n_l, re_eq_n_q_levels, re_eq_n_midplane,     &
          re_eq_finite_pitch
! --- driver interface (used by equilibrium.f90 and the GS element assembly)
public :: re_eq_init, re_eq_update_labels, re_eq_rescale_current,              &
          re_eq_source, re_eq_source_derivs, re_eq_outer_update,               &
          re_eq_write_output, re_eq_finalize, re_eq_done
! --- exposed for the standalone unit test (util/re_equilibrium_prototype)
public :: re_cl_alpha, re_cl_A_edge

! ------------------------------------------------------------------
! --- Namelist input parameters (registered in the model's in1 group)
! ------------------------------------------------------------------
logical            :: re_kinetic_equilibrium = .false. !< master switch: RE drift-surface equilibrium
character(len=256) :: re_eq_dist_file   = 'none'      !< RE momentum distribution table
character(len=32)  :: re_eq_dist_format = 'ekin_xi_w' !< table format (see re_eq_read_distribution)
character(len=256) :: re_eq_q_file      = 'none'      !< target q profile table: psihat_n, q_t
character(len=16)  :: re_eq_match_mode  = 'full_q'    !< 'full_q': match q_t incl. amplitude, I_RE is an output;
                                                      !< 'q_shape': match the shape at prescribed re_eq_I_RE
character(len=16)  :: re_eq_transplant  = 'cumulative'!< outer update variant: 'cumulative' (default) or 'pointwise'
character(len=16)  :: re_eq_map_mode    = 'midplane'  !< label map Ahat<->psihat: 'midplane' (default, the
                                                      !< 2-point midplane average, Eq. 8 of the doc) or
                                                      !< 'contour' (average psihat_n over the WHOLE drift
                                                      !< surface A_s=const by nodal kernel regression;
                                                      !< faithful in shaped geometry, where two midplane
                                                      !< points poorly represent the surface and floor the
                                                      !< achievable q match)
real*8             :: re_eq_I_RE        = 0.d0        !< prescribed RE current [A]: held exactly in
                                                      !< q_shape mode, and the target of the optional
                                                      !< total-current control (re_eq_alpha_current) in
                                                      !< full_q mode
real*8             :: re_eq_xi_min      = 0.9d0       !< minimum |pitch|; abort below (trapped REs out of scope)
real*8             :: re_eq_alpha_out   = 0.3d0       !< under-relaxation of the outer transplant update
real*8             :: re_eq_tol_q       = 1.d-3       !< outer convergence: max|q/q_t - 1|
real*8             :: re_eq_tol_q_soft  = 1.d-2       !< soft tolerance: a stagnated iteration with best error
                                                      !< below this is accepted with a warning (with one common
                                                      !< Nprof and strongly different class drift shifts, exactly
                                                      !< matching q_t can be outside the range of the ansatz)
real*8             :: re_eq_edge_taper  = 0.2d0       !< label width of the linear taper that removes the current
                                                      !< of drift surfaces leaving the domain (Ahat > 1): RE
                                                      !< orbits crossing the wall are lost. A hard clamp instead
                                                      !< (keeping Nprof(1) outside) creates an uncontrollable
                                                      !< halo current that destabilizes the outer iteration; a
                                                      !< too-narrow taper (< ~0.1) makes the transplant gain
                                                      !< vanish abruptly at the edge and the iteration treadmills.
real*8             :: re_eq_l_beam      = 1.d0        !< beam-edge label: the RE current occupies drift surfaces
                                                      !< with Ahat < re_eq_l_beam only, leaving a current-free
                                                      !< (vacuum) annulus between the beam edge and the wall.
                                                      !< For < 1, ALL current-carrying drift orbits are closed
                                                      !< (no scrape-off; the edge taper never activates). In the
                                                      !< zero-drift-orbit limit this label equals the normalized
                                                      !< poloidal flux of the beam edge. Default 1 = wall-limited
                                                      !< beam (previous behaviour). q_t is matched only inside
                                                      !< the beam; q in the vacuum annulus is an output.
real*8             :: re_eq_l_beam_width = 0.1d0      !< label width of the smoothstep roll-off of the beam
                                                      !< current inside [l_beam - width, l_beam]
logical            :: re_eq_absorbing_edge = .false.  !< force Nprof(Ahat = 1) = 0, i.e. treat the last
                                                      !< confined drift surface as an ABSORBING boundary
                                                      !< for the RE density. Off by default. Intended for
                                                      !< cases whose domain extends past the LCFS, where
                                                      !< Ahat = 1 lies in the open-field-line region and
                                                      !< the constraint therefore costs nothing where q_t
                                                      !< is demanded. On a domain that IS the LCFS (no
                                                      !< vacuum region, e.g. re_kin_equil_600) it deletes
                                                      !< current exactly where q_t needs it and is very
                                                      !< expensive -- do not enable there.
                                                      !< NOTE it makes re_eq_edge_taper INERT: the taper
                                                      !< only acts for lraw > 1, and re_eq_nprof_eval
                                                      !< clips to [0,1], so once Nprof(1) = 0 the source
                                                      !< already vanishes at and beyond Ahat = 1.
real*8             :: re_eq_alpha_current = 0.d0      !< under-relaxation of the total-current control:
                                                      !< after each outer update Nprof is scaled by
                                                      !< (|re_eq_I_RE| / |I_RE|)**re_eq_alpha_current.
                                                      !< 0 (default) = OFF; ~0.2 is a sensible start.
                                                      !< WHY a slow rescale and NOT a row in the operator
                                                      !< least squares: at FIXED psi the flux surfaces do
                                                      !< not move, so q ~ 1/I_enc rigidly -- a current row
                                                      !< is then the SAME equation as the outermost q row
                                                      !< with a different right-hand side, and the solve
                                                      !< just splits the difference (tried: it degrades
                                                      !< both criteria at every weight).
                                                      !< The size/current degeneracy q_a ~ a^2/I lives
                                                      !< ACROSS outer iterations, where psi re-solves and
                                                      !< the LCFS moves. Rescaling between iterations lets
                                                      !< the equilibrium slide ALONG that direction: the
                                                      !< GS solve expands the LCFS and q(psihat_n)
                                                      !< recovers, so the q match and the current control
                                                      !< separate instead of competing.
                                                      !< Only meaningful where the LCFS can move (a domain
                                                      !< larger than the plasma). Where the boundary IS
                                                      !< the LCFS (n_open = 0) the area cannot change,
                                                      !< there is no degenerate direction to slide along,
                                                      !< and this will simply fight the q match.
real*8             :: re_eq_op_lambda   = 1.d-2       !< smoothness regularization of the 'operator'
                                                      !< transplant variant (damped least squares on the
                                                      !< relative Nprof correction); unused otherwise
real*8             :: re_eq_ratio_clamp = 2.d0        !< per-iteration clamp of the transplant ratio
integer            :: re_eq_max_it_out  = 50          !< maximum outer iterations
integer            :: re_eq_n_l         = 101         !< number of points of the Nprof(l) table
integer            :: re_eq_n_q_levels  = 80          !< number of psihat levels of the q evaluation
integer            :: re_eq_n_midplane  = 400         !< number of midplane points of the label map
logical            :: re_eq_finite_pitch = .false.    !< use A_s with R*B_phi/B and mu-conserving v_par (not
                                                      !< yet implemented; the small-pitch default neglects
                                                      !< O((p_perp/p_par)^2 * dB/B))

! ------------------------------------------------------------------
! --- Per-class state (SI units; labels/axes in JOREK psi units)
! ------------------------------------------------------------------
integer             :: re_eq_n_class = 0
real*8, allocatable :: re_cl_ekin(:)    !< kinetic energy [eV]
real*8, allocatable :: re_cl_xi(:)      !< pitch p_par/p
real*8, allocatable :: re_cl_w(:)       !< weights, normalized to sum = 1
real*8, allocatable :: re_cl_gamma(:)   !< Lorentz factor
real*8, allocatable :: re_cl_vpar(:)    !< parallel velocity [m/s]
real*8, allocatable :: re_cl_alpha(:)   !< gamma m_e v_par / e [Wb/m]: A_s/e = alpha_s R - psi
real*8, allocatable :: re_cl_A_axis(:)  !< A_s/e at the class drift axis [Wb]
real*8, allocatable :: re_cl_A_edge(:)  !< A_s/e at the outboard midplane edge [Wb]
real*8, allocatable :: re_cl_R_axis(:)  !< R of the class drift axis [m]
real*8, allocatable :: re_cl_Z_axis(:)  !< Z of the class drift axis [m]
real*8, allocatable :: re_cl_edge_frac(:) !< carried-current fraction at Ahat > 0.95: how hard
                                          !< the beam edge sits against the loss boundary

!> common normalized profile function Nprof(l) [m^-2], uniform l grid in [0,1]
real*8, allocatable :: re_nprof_l(:), re_nprof(:)

!> target q profile table (with monotone-cubic slopes: a piecewise-LINEAR
!> target has kinks no smooth equilibrium q can match, which floors the
!> achievable q error at ~ curvature * spacing^2 of the table)
integer             :: re_eq_n_qt = 0
real*8, allocatable :: re_qt_psihat(:), re_qt_q(:), re_qt_slope(:)

!> outer iteration state
integer :: re_eq_outer_iter  = 0        !< current outer iteration (for logging)
logical :: re_eq_initialized = .false.  !< tables read and classes built
logical :: re_eq_labels_ready = .false. !< per-class labels valid for the present psi
real*8  :: re_eq_R_edge  = 0.d0         !< outboard midplane boundary radius [m]
real*8  :: re_eq_psi_bnd = 0.d0         !< boundary psi used in the labels
real*8  :: re_eq_q_err   = 1.d99        !< latest max|q/q_t - 1|
real*8  :: re_eq_I_now   = 0.d0         !< latest RE current [A]
!> best-iterate tracking / stagnation handling of the outer loop
real*8              :: re_eq_best_err_cur = 1.d99 !< best |I_RE/target - 1| so far; only used when the
                                                   !< total-current control is active
real*8              :: re_eq_best_err = 1.d99   !< best IN-BEAM max|q/q_t - 1| so far
                                                 !< (the controllable objective; the verdict
                                                 !<  uses the full-range error, see below)
real*8, allocatable :: re_eq_best_nprof(:)      !< Nprof of the best iterate
integer             :: re_eq_n_stall  = 0       !< outer iterations without improvement
logical :: re_eq_finishing     = .false.        !< best profile restored; final evaluation pass
logical :: re_eq_soft_accepted = .false.        !< finished above tol_q but below tol_q_soft
logical :: re_eq_reverted      = .false.        !< finishing: polish worsened q, best profile restored
logical :: re_eq_done          = .false.        !< finishing verdict rendered; the outer loop must stop
                                                !< (accepted OR given up -- read by the equilibrium driver)

integer, parameter :: RE_EQ_LOG_UNIT = 437  !< unit of re_eq_convergence.log

contains

!=======================================================================
!> Read the input tables, build the classes, write the startup log.
!> Called once (my_id == 0) at the beginning of the equilibrium solve.
subroutine re_eq_init(my_id)
  use phys_module, only: F0, R_geo, amin, FF_coef, T_coef, num_ffprime, num_T
  implicit none
  integer, intent(in) :: my_id
  integer :: s
  real*8  :: B0, d_s

  if (re_eq_initialized) return

  if (my_id .ne. 0) then
    write(*,*) 'ERROR: re_eq_init must only be called on MPI task 0'
    stop 1
  endif

  if (re_eq_finite_pitch) then
    write(*,*) 'ERROR: re_eq_finite_pitch=.true. is not implemented yet; the'
    write(*,*) '       small-pitch form A_s = gamma m v_par R - e psi is used.'
    stop 1
  endif

  if ((trim(re_eq_map_mode) .ne. 'midplane') .and. &
      (trim(re_eq_map_mode) .ne. 'contour')) then
    write(*,*) 'ERROR: re_eq_map_mode must be ''midplane'' or ''contour'', got: ', &
               trim(re_eq_map_mode)
    stop 1
  endif

  ! --- Guard against a NaN trap in the ANALYTIC profile evaluations: the
  ! --- analytic FFprime and temperature routines divide by the shape widths
  ! --- (FF_coef(4)/(8), T_coef(4)) even when the profile amplitude is zero,
  ! --- so a pure-RE equilibrium (FF_0 = FF_1 = 0) needs those FINITE. This
  ! --- only applies when the analytic form is used: with a numerical profile
  ! --- from a file (ffprime_file / T_file -> num_ffprime / num_T set true in
  ! --- read_num_profiles) the shape coefficients are unused, so skip the guard.
  if (.not. num_ffprime) then
    if ((FF_coef(4) .eq. 0.d0) .or. (FF_coef(8) .eq. 0.d0)) then
      write(*,*) 'ERROR: re_eq: FF_coef(4) and FF_coef(8) must be nonzero even for'
      write(*,*) '       a zero-amplitude FFprime profile (FF_0=FF_1=0): the analytic'
      write(*,*) '       profile evaluation divides by them and yields NaN otherwise.'
      write(*,*) '       Use e.g. FF_coef(4)=0.03, FF_coef(5)=5., FF_coef(8)=1.'
      write(*,*) '       (Not needed when FFprime is read from ffprime_file.)'
      stop 1
    endif
  endif
  if ((.not. num_T) .and. (T_coef(4) .eq. 0.d0)) then
    write(*,*) 'ERROR: re_eq: T_coef(4) must be nonzero for the analytic temperature'
    write(*,*) '       profile (see the FF_coef note above). Not needed when the'
    write(*,*) '       temperature is read from T_file.'
    stop 1
  endif

  if ((trim(re_eq_transplant) .ne. 'cumulative') .and. &
      (trim(re_eq_transplant) .ne. 'pointwise')  .and. &
      (trim(re_eq_transplant) .ne. 'operator')) then
    write(*,*) 'ERROR: re_eq_transplant must be ''cumulative'', ''pointwise'' or'
    write(*,*) '       ''operator'', got: ', trim(re_eq_transplant)
    stop 1
  endif

  if ((re_eq_alpha_current .gt. 0.d0) .and. (re_eq_I_RE .eq. 0.d0)) then
    write(*,*) 'ERROR: re_eq: re_eq_alpha_current > 0 requires a nonzero re_eq_I_RE'
    write(*,*) '       (the target of the total-current control).'
    stop 1
  endif
  if ((re_eq_alpha_current .gt. 0.d0) .and. (trim(re_eq_match_mode) .eq. 'q_shape')) then
    write(*,*) 'ERROR: re_eq: re_eq_alpha_current is for full_q. In q_shape the'
    write(*,*) '       current is already held exactly by re_eq_rescale_current,'
    write(*,*) '       and the q amplitude is free -- so the size/current'
    write(*,*) '       degeneracy returns and the control has nothing to pin.'
    stop 1
  endif

  call re_eq_read_distribution(my_id)
  call re_eq_read_q_target(my_id)
  call re_eq_init_nprof(my_id)

  ! --- startup log: the per-class drift parameter d_s = gamma m v_par/(e B0 a)
  !     is the single most useful sanity number (Delta_s ~ d_s q / (r/a))
  B0 = abs(F0) / R_geo
  write(*,*) '*******************************************************'
  write(*,*) '*   kinetic RE drift-surface equilibrium (re_eq)      *'
  write(*,*) '*******************************************************'
  write(*,'(A,I4)')     '   number of RE classes : ', re_eq_n_class
  write(*,'(A,A)')      '   match mode           : ', trim(re_eq_match_mode)
  write(*,'(A,A)')      '   transplant variant   : ', trim(re_eq_transplant)
  write(*,'(A,A)')      '   label map            : ', trim(re_eq_map_mode)
  if (re_eq_alpha_current .gt. 0.d0) &
    write(*,'(A,ES12.4,A,F8.3,A)') '   current control      : I_RE -> ', re_eq_I_RE, &
      ' A  (alpha ', re_eq_alpha_current, ')'
  if (re_eq_absorbing_edge) then
    write(*,'(A)')      '   absorbing edge       : ON  (Nprof(Ahat=1) = 0;'
    write(*,'(A)')      '                          re_eq_edge_taper has no effect)'
  endif
  write(*,*) '    s    E_kin[eV]      xi        weight     gamma      v_par[m/s]     d_s'
  do s = 1, re_eq_n_class
    d_s = re_cl_alpha(s) / (B0 * amin)
    write(*,'(I5,ES13.4,F10.5,ES13.4,F10.3,ES15.6,ES12.3)') &
      s, re_cl_ekin(s), re_cl_xi(s), re_cl_w(s), re_cl_gamma(s), re_cl_vpar(s), d_s
  enddo

  open(RE_EQ_LOG_UNIT, file='re_eq_convergence.log', action='write', status='replace')
  write(RE_EQ_LOG_UNIT,'(A)') '# outer  inner_iters  max|q/qt-1|   q_err_in_beam   I_RE[A]        I_err          max_edge_fraction'

  re_eq_outer_iter  = 0
  re_eq_initialized = .true.

end subroutine re_eq_init


!=======================================================================
!> Read the RE momentum-space distribution table. One reader per format;
!> adding a new table format means adding one case here.
!> Format 'ekin_xi_w': one class per line, columns
!>   E_kin [eV]   xi = p_par/p   weight (relative, normalized internally),
!> '#' starts a comment.
subroutine re_eq_read_distribution(my_id)
  use tr_module
  implicit none
  integer, intent(in) :: my_id
  integer, parameter  :: max_classes = 10000
  integer :: iunit, ierr, n, ipos, s
  real*8  :: cols(3), wsum, mec2_eV
  real*8  :: tmp_e(max_classes), tmp_xi(max_classes), tmp_w(max_classes)
  character(len=512) :: line

  select case (trim(re_eq_dist_format))
  case ('ekin_xi_w')
    ! fall through to the reader below
  case default
    write(*,*) 'ERROR: unknown re_eq_dist_format: ', trim(re_eq_dist_format)
    stop 1
  end select

  iunit = 438
  open(iunit, file=trim(re_eq_dist_file), status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) 'ERROR: cannot open re_eq_dist_file: ', trim(re_eq_dist_file)
    stop 1
  endif

  n = 0
  do
    read(iunit,'(A)',iostat=ierr) line
    if (ierr .ne. 0) exit
    ipos = index(line, '#')
    if (ipos .gt. 0) line = line(1:ipos-1)
    if (len_trim(line) .eq. 0) cycle
    n = n + 1
    if (n .gt. max_classes) then
      write(*,*) 'ERROR: more than ', max_classes, ' classes in ', trim(re_eq_dist_file)
      stop 1
    endif
    read(line,*,iostat=ierr) cols
    if (ierr .ne. 0) then
      write(*,*) 'ERROR: expected 3 columns (E_kin[eV] xi weight) in ', &
                 trim(re_eq_dist_file), ' line: ', trim(line)
      stop 1
    endif
    tmp_e(n) = cols(1); tmp_xi(n) = cols(2); tmp_w(n) = cols(3)
  enddo
  close(iunit)

  if (n .eq. 0) then
    write(*,*) 'ERROR: no distribution classes found in ', trim(re_eq_dist_file)
    stop 1
  endif

  re_eq_n_class = n
  call tr_allocate(re_cl_ekin,  1, n, "re_cl_ekin",  CAT_GRID)
  call tr_allocate(re_cl_xi,    1, n, "re_cl_xi",    CAT_GRID)
  call tr_allocate(re_cl_w,     1, n, "re_cl_w",     CAT_GRID)
  call tr_allocate(re_cl_gamma, 1, n, "re_cl_gamma", CAT_GRID)
  call tr_allocate(re_cl_vpar,  1, n, "re_cl_vpar",  CAT_GRID)
  call tr_allocate(re_cl_alpha, 1, n, "re_cl_alpha", CAT_GRID)
  call tr_allocate(re_cl_A_axis,1, n, "re_cl_A_axis",CAT_GRID)
  call tr_allocate(re_cl_A_edge,1, n, "re_cl_A_edge",CAT_GRID)
  call tr_allocate(re_cl_R_axis,1, n, "re_cl_R_axis",CAT_GRID)
  call tr_allocate(re_cl_Z_axis,1, n, "re_cl_Z_axis",CAT_GRID)
  call tr_allocate(re_cl_edge_frac, 1, n, "re_cl_edge_frac", CAT_GRID)

  re_cl_ekin(1:n) = tmp_e(1:n)
  re_cl_xi(1:n)   = tmp_xi(1:n)

  ! --- checks and derived quantities
  wsum = sum(tmp_w(1:n))
  if (wsum .le. 0.d0) then
    write(*,*) 'ERROR: distribution weights must be positive in ', trim(re_eq_dist_file)
    stop 1
  endif
  re_cl_w(1:n) = tmp_w(1:n) / wsum

  do s = 1, n
    if (abs(re_cl_xi(s)) .lt. re_eq_xi_min) then
      write(*,'(A,I4,A,F8.4,A,F8.4)') ' ERROR: RE class ', s, ' has |xi| = ', &
        abs(re_cl_xi(s)), ' < re_eq_xi_min = ', re_eq_xi_min
      write(*,*) '        Trapped-particle invariants are out of scope for the'
      write(*,*) '        drift-surface equilibrium (strongly passing REs assumed).'
      stop 1
    endif
    if (re_cl_ekin(s) .le. 0.d0) then
      write(*,*) 'ERROR: non-positive RE class energy in ', trim(re_eq_dist_file)
      stop 1
    endif
  enddo

  mec2_eV = MASS_ELECTRON * SPEED_OF_LIGHT**2 / EL_CHG
  do s = 1, n
    re_cl_gamma(s) = 1.d0 + re_cl_ekin(s) / mec2_eV
    re_cl_vpar(s)  = re_cl_xi(s) * SPEED_OF_LIGHT &
                     * sqrt(1.d0 - 1.d0/re_cl_gamma(s)**2)
    re_cl_alpha(s) = re_cl_gamma(s) * MASS_ELECTRON * re_cl_vpar(s) / EL_CHG
  enddo

  re_cl_A_axis = 0.d0;  re_cl_A_edge = 0.d0
  re_cl_R_axis = 0.d0;  re_cl_Z_axis = 0.d0
  re_cl_edge_frac = 0.d0

end subroutine re_eq_read_distribution


!=======================================================================
!> Read the target q profile table: two columns psihat_n, q_t.
subroutine re_eq_read_q_target(my_id)
  use tr_module
  implicit none
  integer, intent(in) :: my_id
  integer, parameter  :: max_rows = 10000
  integer :: iunit, ierr, n, ipos
  real*8  :: cols(2), tmp_p(max_rows), tmp_q(max_rows)
  character(len=512) :: line

  iunit = 438
  open(iunit, file=trim(re_eq_q_file), status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) then
    write(*,*) 'ERROR: cannot open re_eq_q_file: ', trim(re_eq_q_file)
    stop 1
  endif
  n = 0
  do
    read(iunit,'(A)',iostat=ierr) line
    if (ierr .ne. 0) exit
    ipos = index(line, '#')
    if (ipos .gt. 0) line = line(1:ipos-1)
    if (len_trim(line) .eq. 0) cycle
    n = n + 1
    if (n .gt. max_rows) then
      write(*,*) 'ERROR: more than ', max_rows, ' rows in ', trim(re_eq_q_file)
      stop 1
    endif
    read(line,*,iostat=ierr) cols
    if (ierr .ne. 0) then
      write(*,*) 'ERROR: expected 2 columns (psihat_n q) in ', trim(re_eq_q_file)
      stop 1
    endif
    tmp_p(n) = cols(1); tmp_q(n) = cols(2)
  enddo
  close(iunit)

  if (n .lt. 2) then
    write(*,*) 'ERROR: target q profile needs at least 2 rows: ', trim(re_eq_q_file)
    stop 1
  endif

  re_eq_n_qt = n
  call tr_allocate(re_qt_psihat, 1, n, "re_qt_psihat", CAT_GRID)
  call tr_allocate(re_qt_q,      1, n, "re_qt_q",      CAT_GRID)
  call tr_allocate(re_qt_slope,  1, n, "re_qt_slope",  CAT_GRID)
  re_qt_psihat(1:n) = tmp_p(1:n)
  re_qt_q(1:n)      = tmp_q(1:n)

  ! --- monotone-cubic slopes (Fritsch-Butland): harmonic mean of adjacent
  !     secants of the same sign, zero otherwise
  do ipos = 2, n-1
    cols(1) = (re_qt_q(ipos)   - re_qt_q(ipos-1)) / (re_qt_psihat(ipos)   - re_qt_psihat(ipos-1))
    cols(2) = (re_qt_q(ipos+1) - re_qt_q(ipos))   / (re_qt_psihat(ipos+1) - re_qt_psihat(ipos))
    if (cols(1)*cols(2) .gt. 0.d0) then
      re_qt_slope(ipos) = 2.d0*cols(1)*cols(2) / (cols(1) + cols(2))
    else
      re_qt_slope(ipos) = 0.d0
    endif
  enddo
  re_qt_slope(1) = (re_qt_q(2) - re_qt_q(1)) / (re_qt_psihat(2) - re_qt_psihat(1))
  re_qt_slope(n) = (re_qt_q(n) - re_qt_q(n-1)) / (re_qt_psihat(n) - re_qt_psihat(n-1))

end subroutine re_eq_read_q_target


!=======================================================================
!> Monotone-cubic (Hermite) evaluation of the target q at psihat, clipped
!> to the table range.
function re_eq_qt_eval(psihat) result(qval)
  implicit none
  real*8, intent(in) :: psihat
  real*8             :: qval, x, h, t, h00, h10, h01, h11
  integer            :: k
  x = min(max(psihat, re_qt_psihat(1)), re_qt_psihat(re_eq_n_qt))
  do k = 2, re_eq_n_qt
    if (x .le. re_qt_psihat(k)) exit
  enddo
  k = min(k, re_eq_n_qt)
  h = re_qt_psihat(k) - re_qt_psihat(k-1)
  t = (x - re_qt_psihat(k-1)) / h
  h00 = (1.d0 + 2.d0*t) * (1.d0 - t)**2
  h10 = t * (1.d0 - t)**2
  h01 = t*t * (3.d0 - 2.d0*t)
  h11 = t*t * (t - 1.d0)
  qval = h00*re_qt_q(k-1) + h10*h*re_qt_slope(k-1) &
       + h01*re_qt_q(k)   + h11*h*re_qt_slope(k)
end function re_eq_qt_eval


!=======================================================================
!> Piecewise-linear evaluation of the common profile function Nprof at the
!> label l, clipped to [0,1].
function re_eq_nprof_eval(l) result(nval)
  implicit none
  real*8, intent(in) :: l
  real*8             :: nval, x, dl
  integer            :: k
  x  = min(max(l, 0.d0), 1.d0)
  dl = re_nprof_l(2) - re_nprof_l(1)
  k  = min(int(x/dl) + 1, re_eq_n_l - 1)
  nval = re_nprof(k) + (re_nprof(k+1) - re_nprof(k)) * (x - re_nprof_l(k)) / dl
end function re_eq_nprof_eval


!=======================================================================
!> The multiplicative edge factor of the effective profile: the C1
!> smoothstep beam-edge envelope (current confined below re_eq_l_beam,
!> roll-off of width re_eq_l_beam_width) times the C1 smoothstep wall
!> taper (open drift orbits beyond lraw = 1 carry no current). Applied at
!> EVALUATION time, exactly once -- the stored Nprof table stays raw.
!> (Baking the envelope into the table is NOT idempotent: re-applying it
!> after every transplant update compounds to envelope^n and sharpens the
!> roll-off towards a step, which the transplant then fights.)
function re_eq_edge_factor(lraw) result(f)
  implicit none
  real*8, intent(in) :: lraw
  real*8             :: f, t
  f = 1.d0
  if (re_eq_l_beam .lt. 1.d0) then
    t = (lraw - (re_eq_l_beam - re_eq_l_beam_width)) &
        / max(re_eq_l_beam_width, 1.d-12)
    t = min(max(t, 0.d0), 1.d0)
    f = 1.d0 - t*t*(3.d0 - 2.d0*t)
  endif
  if (lraw .gt. 1.d0) then
    t = min((lraw - 1.d0) / max(re_eq_edge_taper, 1.d-12), 1.d0)
    f = f * (1.d0 - t*t*(3.d0 - 2.d0*t))
  endif
end function re_eq_edge_factor


!=======================================================================
!> d/dl of the edge factor (for the source derivatives).
function re_eq_edge_factor_deriv(lraw) result(df)
  implicit none
  real*8, intent(in) :: lraw
  real*8             :: df, t, w, f_env, f_tap, d_env, d_tap
  f_env = 1.d0;  d_env = 0.d0
  if (re_eq_l_beam .lt. 1.d0) then
    w = max(re_eq_l_beam_width, 1.d-12)
    t = (lraw - (re_eq_l_beam - w)) / w
    if ((t .gt. 0.d0) .and. (t .lt. 1.d0)) then
      f_env = 1.d0 - t*t*(3.d0 - 2.d0*t)
      d_env = -6.d0*t*(1.d0 - t) / w
    else if (t .ge. 1.d0) then
      f_env = 0.d0
    endif
  endif
  f_tap = 1.d0;  d_tap = 0.d0
  if (lraw .gt. 1.d0) then
    w = max(re_eq_edge_taper, 1.d-12)
    t = (lraw - 1.d0) / w
    if (t .lt. 1.d0) then
      f_tap = 1.d0 - t*t*(3.d0 - 2.d0*t)
      d_tap = -6.d0*t*(1.d0 - t) / w
    else
      f_tap = 0.d0
    endif
  endif
  df = d_env*f_tap + f_env*d_tap
end function re_eq_edge_factor_deriv


!=======================================================================
!> Nprof at the RAW (unclipped) label: the raw table times the edge factor
!> (beam envelope and wall taper). This function MUST be used wherever the
!> source density is evaluated (and its counterpart in
!> particles/initialisers/initialisers_RE.f90 kept in sync).
function re_eq_nprof_at(lraw) result(nval)
  implicit none
  real*8, intent(in) :: lraw
  real*8             :: nval
  nval = re_eq_nprof_eval(lraw) * re_eq_edge_factor(lraw)
end function re_eq_nprof_at


!=======================================================================
!> Initial guess of Nprof from the target q in cylindrical approximation:
!>   q = r B0 / (R0 Btheta)  =>  I(r) = 2 pi r^2 B0 / (mu0 R0 q),
!>   j(r) = I'/(2 pi r) = (B0/(mu0 R0)) (2/q - (r/q^2) dq/dr)
!> evaluated analytically (a finite-difference I'/r blows up at the axis),
!> with the label approximated by the cylindrical psihat(r). A crude guess
!> only costs a few extra outer iterations.
subroutine re_eq_init_nprof(my_id)
  use tr_module
  use phys_module, only: F0, R_geo, amin
  implicit none
  integer, intent(in) :: my_id
  integer, parameter  :: nr = 512
  integer :: i, k, ipass
  real*8  :: B0, r(nr), q(nr), ph(nr), Ienc(nr), Bth(nr), j(nr), dq_dr(nr)
  real*8  :: psi_pol, vbar, l, dl

  call tr_allocate(re_nprof_l, 1, re_eq_n_l, "re_nprof_l", CAT_GRID)
  call tr_allocate(re_nprof,   1, re_eq_n_l, "re_nprof",   CAT_GRID)
  call tr_allocate(re_eq_best_nprof, 1, re_eq_n_l, "re_eq_best_nprof", CAT_GRID)
  do k = 1, re_eq_n_l
    re_nprof_l(k) = dble(k-1) / dble(re_eq_n_l - 1)
  enddo
  re_eq_best_err      = 1.d99
  re_eq_best_err_cur  = 1.d99
  re_eq_n_stall       = 0
  re_eq_finishing     = .false.
  re_eq_soft_accepted = .false.
  re_eq_reverted      = .false.
  re_eq_done          = .false.

  B0 = abs(F0) / R_geo
  do i = 1, nr
    r(i)  = amin * dble(i) / dble(nr)
    ph(i) = (r(i)/amin)**2          ! first pass: psihat ~ (r/a)^2
  enddo

  do ipass = 1, 2
    do i = 1, nr
      ! use |q_t|: this cylindrical guess sets the current MAGNITUDE profile
      ! (Nprof >= 0), independent of the sign convention of q (which only
      ! reflects the current direction, i.e. the RE pitch sign). A negative
      ! q_t would otherwise give a negative j that is clamped to zero below.
      q(i)    = abs(re_eq_qt_eval(ph(i)))
      Ienc(i) = 2.d0*PI * r(i)**2 * B0 / (MU_ZERO * R_geo * q(i))
      Bth(i)  = MU_ZERO * Ienc(i) / (2.d0*PI * r(i))
    enddo
    ! cylindrical poloidal flux -> psihat(r) for the second pass
    psi_pol = 0.d0
    ph(1)   = 0.d0
    do i = 2, nr
      psi_pol = psi_pol + 0.5d0*(Bth(i)+Bth(i-1)) * (r(i)-r(i-1)) * R_geo
      ph(i)   = psi_pol
    enddo
    ph = ph / ph(nr)
  enddo

  do i = 2, nr-1
    dq_dr(i) = (q(i+1) - q(i-1)) / (r(i+1) - r(i-1))
  enddo
  dq_dr(1)  = dq_dr(2)
  dq_dr(nr) = dq_dr(nr-1)
  do i = 1, nr
    j(i) = B0 / (MU_ZERO * R_geo) * (2.d0/q(i) - r(i)*dq_dr(i)/q(i)**2)
    j(i) = max(j(i), 0.d0)
  enddo

  ! --- Nprof^0(l) ~ R0 j(r(l)) / (e vbar), l approximated by psihat(r)
  vbar = abs(sum(re_cl_w(1:re_eq_n_class) * re_cl_vpar(1:re_eq_n_class)))
  do k = 1, re_eq_n_l
    l = re_nprof_l(k)
    ! invert ph(r): find first i with ph(i) >= l
    ! (no compound condition: Fortran .and. does not short-circuit, and
    !  ph(i-1) must not be evaluated for i = 1)
    i = nr
    do while (i .gt. 1)
      if (ph(i-1) .lt. l) exit
      i = i - 1
    enddo
    if (i .eq. 1) then
      re_nprof(k) = j(1)
    else
      re_nprof(k) = j(i-1) + (j(i)-j(i-1)) * (l-ph(i-1)) / (ph(i)-ph(i-1))
    endif
    re_nprof(k) = re_nprof(k) * R_geo / (EL_CHG * vbar)
  enddo

  call re_eq_apply_beam_envelope()

end subroutine re_eq_init_nprof


!=======================================================================
!> A of the last drift surface of a class that does NOT reach the loss
!> boundary: the extremum of A_s = alpha_s R - psi over the BOUNDARY nodes.
!>
!> A drift orbit is lost when its surface reaches the WALL (limiter /
!> divertor target), not when it crosses a flux surface, so the critical
!> label is a property of the boundary contour and of nothing else. Written
!> this way the definition makes NO assumption about the boundary: the
!> previous form (alpha * maxval(R) - ES%psi_bnd) is the special case
!> psi = psi_bnd all along the boundary and is reproduced EXACTLY there.
!>
!> That special case does not hold in general. On a vessel-shaped fixed
!> boundary (psi prescribed and varying along the contour) the old form is
!> over-restrictive by ~ |alpha| * (R range) / |psi_bnd - psi_axis|, i.e. it
!> truncates the beam inside the true last confined surface by an amount
!> that grows LINEARLY with the class energy -- 19% of the label range at
!> 10 MeV on the JET-like limiter test case, where it leaves an unphysical
!> current-free gap between the beam and the wall around the whole contour.
!> It also forced the diverted case to treat the SEPARATRIX as the loss
!> boundary, whereas a confined orbit may legitimately excurse into the
!> scrape-off layer.
!>
!> sgn < 0 means A has a MINIMUM at the drift axis and increases outwards,
!> so the confined surfaces are those below the boundary minimum (and the
!> other way round for sgn > 0). It is the same branch selector that
!> re_eq_find_drift_axis uses, keyed on psi and NOT on the sign of alpha:
!> alpha and the sign of j_phi (hence of psi_axis - psi_bnd) both follow
!> sign(v_par), so the two flips cancel and either pitch sign is handled.
function re_eq_A_edge_bnd(alpha, node_list, bnd_node_list) result(A_edge)
  use data_structure
  use equil_info, only: ES
  use mod_model_settings, only: var_psi
  implicit none
  real*8,                      intent(in) :: alpha
  type (type_node_list),       intent(in) :: node_list
  type (type_bnd_node_list),   intent(in) :: bnd_node_list
  real*8  :: A_edge, A_i, sgn
  integer :: i, iv

  sgn = 1.d0
  if (ES%psi_axis .gt. ES%psi_bnd) sgn = -1.d0     ! psi max at axis -> A min at axis

  A_edge = -sgn * 1.d99
  do i = 1, bnd_node_list%n_bnd_nodes
    iv   = bnd_node_list%bnd_node(i)%index_jorek
    A_i  = alpha * node_list%node(iv)%x(1,1,1) - node_list%node(iv)%values(1,1,var_psi)
    if (sgn * A_i .gt. sgn * A_edge) A_edge = A_i
  enddo

end function re_eq_A_edge_bnd


!=======================================================================
!> Update the per-class invariant labels: locate the drift axis of every
!> class (the extremum of A_s, displaced from the magnetic axis) and the
!> loss-boundary value (re_eq_A_edge_bnd). MUST be called every time psi
!> changes during the Picard iteration -- both move as psi converges.
subroutine re_eq_update_labels(my_id, node_list, element_list, bnd_node_list)
  use data_structure
  use equil_info, only: ES
  implicit none
  integer,                    intent(in) :: my_id
  type (type_node_list),      intent(in) :: node_list
  type (type_element_list),   intent(in) :: element_list
  type (type_bnd_node_list),  intent(in) :: bnd_node_list
  integer :: s, ifail
  real*8  :: R0s, Z0s, R_ax, Z_ax, A_ax

  ! --- outboard midplane edge radius. No longer used for the A_edge label
  !     normalization (that now comes from the boundary contour itself, see
  !     re_eq_A_edge_bnd); retained because the midplane label map scans out
  !     to it and the startup report quotes it.
  if (ES%xpoint .and. (ES%LCFS_a .gt. 0.d0)) then
    re_eq_R_edge = ES%LCFS_Rgeo + ES%LCFS_a
  else
    re_eq_R_edge = maxval(node_list%node(1:node_list%n_nodes)%x(1,1,1))
  endif
  re_eq_psi_bnd = ES%psi_bnd

  do s = 1, re_eq_n_class
    ! start the search from the previous drift axis (or the magnetic axis)
    if (re_eq_labels_ready) then
      R0s = re_cl_R_axis(s);  Z0s = re_cl_Z_axis(s)
    else
      R0s = ES%R_axis;        Z0s = ES%Z_axis
    endif
    call re_eq_find_drift_axis(node_list, element_list, re_cl_alpha(s), &
                               R0s, Z0s, R_ax, Z_ax, A_ax, ifail)
    if (ifail .ne. 0) then
      write(*,'(A,I4,A)') ' WARNING: re_eq: drift-axis search failed for class ', s, &
                          ' -- keeping previous axis'
      if (.not. re_eq_labels_ready) then
        write(*,*) 'ERROR: re_eq: no valid drift axis for class ', s
        stop 1
      endif
    else
      re_cl_R_axis(s) = R_ax
      re_cl_Z_axis(s) = Z_ax
      re_cl_A_axis(s) = A_ax
    endif
    re_cl_A_edge(s) = re_eq_A_edge_bnd(re_cl_alpha(s), node_list, bnd_node_list)
    if (abs(re_cl_A_edge(s) - re_cl_A_axis(s)) .le. 0.d0) then
      write(*,*) 'ERROR: re_eq: degenerate label normalization for class ', s
      stop 1
    endif
  enddo
  re_eq_labels_ready = .true.

end subroutine re_eq_update_labels


!=======================================================================
!> Locate the drift axis of one class: the interior extremum of
!> A(R,Z) = alpha R - psi. Coarse node scan for the extremal node, then a
!> local least-squares paraboloid fit of the nodal A values around it;
!> the fit gives the sub-element extremum position and value analytically.
!> The extremum type of A is opposite to psi's (A ~ -psi near the axis;
!> the alpha R term only shifts the extremum).
!>
!> This deliberately uses ONLY nodal data: point location (find_RZ) is
!> unreliable in the whole neighbourhood of the degenerate polar-grid
!> centre -- exactly where the drift axes of weakly-shifted classes live.
subroutine re_eq_find_drift_axis(node_list, element_list, alpha, R0, Z0, &
                                 R_ax, Z_ax, A_ax, ifail)
  use data_structure
  use equil_info, only: ES
  use phys_module, only: amin
  use mod_model_settings, only: var_psi
  implicit none
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  real*8,                   intent(in)  :: alpha    !< gamma m v_par/e of the class [Wb/m]
  real*8,                   intent(in)  :: R0, Z0   !< unused (kept for interface stability)
  real*8,                   intent(out) :: R_ax, Z_ax
  real*8,                   intent(out) :: A_ax     !< A/e at the drift axis [Wb]
  integer,                  intent(out) :: ifail

  integer, parameter :: n_fit_min = 10
  real*8  :: sgn, A_node, A_best, R_scan, Z_scan, r_fit, dx, dy, d2
  real*8  :: M(6,6), rhs(6), c(6), row(6), det, ddx, ddy, phn_raw
  integer :: i, k, l, n_fit, i_pass

  ifail = 0

  ! --- coarse node scan: extremum of A over the grid nodes; minimum if psi
  !     has its maximum at the axis and vice versa
  sgn = 1.d0
  if (ES%psi_axis .gt. ES%psi_bnd) sgn = -1.d0    ! psi max at axis -> A min
  A_best = -1.d99
  R_scan = R0;  Z_scan = Z0
  do i = 1, node_list%n_nodes
    ! diverted grid: the drift axis is interior, so skip open-region nodes
    ! (scrape-off / private flux, beyond the separatrix) as candidates --
    ! at large outboard R the alpha*R term can otherwise make A extremal on
    ! an open field line. Non-X-point path unchanged (ES%xpoint=.false.).
    if (ES%xpoint) then
      phn_raw = (node_list%node(i)%values(1,1,var_psi) - ES%psi_axis) &
                / (ES%psi_bnd - ES%psi_axis)
      if (phn_raw .gt. 1.d0 + 1.d-6) cycle
    endif
    A_node = sgn * (alpha * node_list%node(i)%x(1,1,1) - node_list%node(i)%values(1,1,var_psi))
    if (A_node .gt. A_best) then
      A_best = A_node
      R_scan = node_list%node(i)%x(1,1,1)
      Z_scan = node_list%node(i)%x(1,1,2)
    endif
  enddo

  ! --- least-squares paraboloid A ~ c1 + c2 x + c3 y + c4 x^2 + c5 xy + c6 y^2
  !     over the nodes within r_fit of the scan extremum (x = R-R_scan,
  !     y = Z-Z_scan); grow the radius until enough nodes participate
  r_fit = 0.05d0 * amin
  do i_pass = 1, 8
    M = 0.d0;  rhs = 0.d0;  n_fit = 0
    do i = 1, node_list%n_nodes
      dx = node_list%node(i)%x(1,1,1) - R_scan
      dy = node_list%node(i)%x(1,1,2) - Z_scan
      d2 = dx*dx + dy*dy
      if (d2 .gt. r_fit*r_fit) cycle
      n_fit  = n_fit + 1
      A_node = alpha * node_list%node(i)%x(1,1,1) - node_list%node(i)%values(1,1,var_psi)
      row = (/ 1.d0, dx, dy, dx*dx, dx*dy, dy*dy /)
      do k = 1, 6
        do l = 1, 6
          M(k,l) = M(k,l) + row(k)*row(l)
        enddo
        rhs(k) = rhs(k) + row(k)*A_node
      enddo
    enddo
    if (n_fit .ge. n_fit_min) exit
    r_fit = 2.d0 * r_fit
  enddo
  if (n_fit .lt. n_fit_min) then
    ifail = 1
    return
  endif

  call re_eq_solve6(M, rhs, c, ifail)
  if (ifail .ne. 0) return

  ! --- extremum of the paraboloid: grad = 0
  det = 4.d0*c(4)*c(6) - c(5)*c(5)
  if (abs(det) .le. 1.d-30) then
    ! degenerate fit: keep the scan node itself
    ddx = 0.d0;  ddy = 0.d0
  else
    ddx = (-2.d0*c(6)*c(2) + c(5)*c(3)) / det
    ddy = (-2.d0*c(4)*c(3) + c(5)*c(2)) / det
  endif
  ! clamp the displacement to the fit region
  if (abs(ddx) .gt. r_fit) ddx = sign(r_fit, ddx)
  if (abs(ddy) .gt. r_fit) ddy = sign(r_fit, ddy)

  R_ax = R_scan + ddx
  Z_ax = Z_scan + ddy
  A_ax = c(1) + c(2)*ddx + c(3)*ddy + c(4)*ddx*ddx + c(5)*ddx*ddy + c(6)*ddy*ddy

end subroutine re_eq_find_drift_axis


!=======================================================================
!> Solve a 6x6 linear system by Gaussian elimination with partial
!> pivoting (normal equations of the paraboloid fit).
subroutine re_eq_solve6(A_in, b_in, x, ifail)
  implicit none
  real*8,  intent(in)  :: A_in(6,6), b_in(6)
  real*8,  intent(out) :: x(6)
  integer, intent(out) :: ifail
  real*8  :: A(6,6), b(6), piv, fac
  integer :: i, j, k, ip

  A = A_in;  b = b_in;  ifail = 0
  do k = 1, 6
    ip = k
    do i = k+1, 6
      if (abs(A(i,k)) .gt. abs(A(ip,k))) ip = i
    enddo
    if (abs(A(ip,k)) .le. 1.d-300) then
      ifail = 4
      return
    endif
    if (ip .ne. k) then
      do j = 1, 6
        piv = A(k,j);  A(k,j) = A(ip,j);  A(ip,j) = piv
      enddo
      piv = b(k);  b(k) = b(ip);  b(ip) = piv
    endif
    do i = k+1, 6
      fac = A(i,k) / A(k,k)
      A(i,k:6) = A(i,k:6) - fac * A(k,k:6)
      b(i)     = b(i)     - fac * b(k)
    enddo
  enddo
  do k = 6, 1, -1
    x(k) = (b(k) - sum(A(k,k+1:6)*x(k+1:6))) / A(k,k)
  enddo
end subroutine re_eq_solve6


!=======================================================================
!> psi and its (R,Z) gradient at an arbitrary point, via element search
!> plus finite-element interpolation.
!>
!> Robustness: find_RZ legitimately fails (ifail=99) at degenerate points
!> of the element mapping -- most notably the polar-grid centre node, which
!> is exactly where the magnetic axis of a circular fixed-boundary case
!> sits, and points exactly on the domain boundary. Those locations are
!> measure-zero, so on failure the evaluation is retried at small spatial
!> offsets (1e-5 a, then 1e-3 a); the interpolation error introduced is
!> negligible for the axis search and the label map.
subroutine re_eq_grad_psi(node_list, element_list, R, Z, psi, dpsi_dR, dpsi_dZ, ifail)
  use data_structure
  use mod_interp, only: interp_PRZ
  use mod_model_settings, only: var_psi
  use phys_module, only: amin
  implicit none
  interface
    subroutine find_RZ(node_list,element_list,R_find,Z_find,R_out,Z_out,ielm_out,s_out,t_out,ifail)
      use data_structure
      type (type_node_list)    :: node_list
      type (type_element_list) :: element_list
      real*8    :: R_find, Z_find, R_out, Z_out, s_out, t_out
      integer   :: ielm_out, ifail
    end subroutine find_RZ
  end interface
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  real*8,                   intent(in)  :: R, Z
  real*8,                   intent(out) :: psi, dpsi_dR, dpsi_dZ
  integer,                  intent(out) :: ifail

  integer :: i_elm, i_try
  real*8  :: R_out, Z_out, s, t, xjac, R_try, Z_try, eps
  real*8  :: P(1), P_s(1), P_t(1), P_phi(1)
  real*8  :: RR, R_s, R_t, ZZ, Z_s, Z_t
  ! offset pattern: the point itself, then 4 diagonal neighbours at two radii
  real*8, parameter :: off_R(9) = (/ 0.d0,  1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0,  1.d0, -1.d0 /)
  real*8, parameter :: off_Z(9) = (/ 0.d0,  1.d0,  1.d0, -1.d0, -1.d0,  1.d0,  1.d0, -1.d0, -1.d0 /)

  do i_try = 1, 9
    eps = 1.d-5 * amin
    if (i_try .ge. 6) eps = 1.d-3 * amin
    R_try = R + off_R(i_try) * eps
    Z_try = Z + off_Z(i_try) * eps
    call find_RZ(node_list, element_list, R_try, Z_try, R_out, Z_out, i_elm, s, t, ifail)
    if (ifail .eq. 0) exit
  enddo
  if (ifail .ne. 0) return

  call interp_PRZ(node_list, element_list, i_elm, [var_psi], 1, s, t, 0.d0, &
                  P, P_s, P_t, P_phi, RR, R_s, R_t, ZZ, Z_s, Z_t)

  xjac = R_s*Z_t - R_t*Z_s
  if (abs(xjac) .le. 1.d-30) then
    ifail = 3
    return
  endif
  psi     = P(1)
  dpsi_dR = (  P_s(1)*Z_t - P_t(1)*Z_s) / xjac
  dpsi_dZ = (- P_s(1)*R_t + P_t(1)*R_s) / xjac

end subroutine re_eq_grad_psi


!=======================================================================
!> The RE contribution to the Grad-Shafranov current variable
!> zj = Delta* psi (JOREK convention), at one point (psi, R):
!>   S_RE = mu0 e sum_s v_par,s w_s Nprof(Ahat_s)
!> Note there is no explicit R factor left (the 1/R of the class density
!> cancels the R of the GS right-hand side); the element assembly adds
!> S_RE / R to its rhs integrand. This routine accumulates no diagnostics;
!> the beam-edge bookkeeping is done in re_eq_total_current on the grid.
function re_eq_source(psi, R) result(S)
  implicit none
  real*8, intent(in) :: psi, R
  real*8             :: S, lhat
  integer            :: is
  S = 0.d0
  do is = 1, re_eq_n_class
    lhat = (re_cl_alpha(is)*R - psi - re_cl_A_axis(is)) &
           / (re_cl_A_edge(is) - re_cl_A_axis(is))
    S = S + re_cl_vpar(is) * re_cl_w(is) * re_eq_nprof_at(lhat)
  enddo
  S = MU_ZERO * EL_CHG * S
end function re_eq_source


!=======================================================================
!> The RE source and its psi- and R-derivatives (for filling the current
!> variable zj and its derivative degrees of freedom after the solve).
!> The derivative of the piecewise-linear Nprof is its slope in the
!> containing interval (zero outside [0,1]).
subroutine re_eq_source_derivs(psi, R, S, dS_dpsi, dS_dR)
  implicit none
  real*8, intent(in)  :: psi, R
  real*8, intent(out) :: S, dS_dpsi, dS_dR
  real*8  :: lhat, denom, dl, slope, cw, Nval
  integer :: is, k
  S = 0.d0;  dS_dpsi = 0.d0;  dS_dR = 0.d0
  dl = re_nprof_l(2) - re_nprof_l(1)
  do is = 1, re_eq_n_class
    denom = re_cl_A_edge(is) - re_cl_A_axis(is)
    lhat  = (re_cl_alpha(is)*R - psi - re_cl_A_axis(is)) / denom
    cw    = re_cl_vpar(is) * re_cl_w(is)
    Nval  = re_eq_nprof_at(lhat)
    S     = S + cw * Nval
    ! d/dl of (raw table) * (edge factor): product rule
    if ((lhat .gt. 0.d0) .and. (lhat .lt. 1.d0)) then
      k     = min(int(lhat/dl) + 1, re_eq_n_l - 1)
      slope = (re_nprof(k+1) - re_nprof(k)) / dl * re_eq_edge_factor(lhat)
    else
      slope = 0.d0
    endif
    slope = slope + re_eq_nprof_eval(lhat) * re_eq_edge_factor_deriv(lhat)
    dS_dpsi = dS_dpsi + cw * slope * (-1.d0/denom)
    dS_dR   = dS_dR   + cw * slope * (re_cl_alpha(is)/denom)
  enddo
  S       = MU_ZERO * EL_CHG * S
  dS_dpsi = MU_ZERO * EL_CHG * dS_dpsi
  dS_dR   = MU_ZERO * EL_CHG * dS_dR
end subroutine re_eq_source_derivs


!=======================================================================
!> Total RE current [A] carried by the source on the present psi:
!>   I_RE = int j_phi dA = -e sum_s v_par,s w_s int Nprof(Ahat_s)/R dA
!> evaluated by Gaussian quadrature over all elements. Also refreshes the
!> per-class beam-edge fractions (carried current at Ahat > 0.95).
subroutine re_eq_total_current(my_id, node_list, element_list, I_RE)
  use mod_parameters, only: n_vertex_max, n_degrees
  use data_structure
  use gauss
  use basis_at_gaussian
  use mod_model_settings, only: var_psi
  implicit none
  integer,                  intent(in)  :: my_id
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  real*8,                   intent(out) :: I_RE

  integer :: i, iv, kv, kf, ms, mt, s
  real*8  :: x_g, y_g, x_s, x_t, y_s, y_t, eq_g, eq_s, eq_t
  real*8  :: xjac, wst, lhat, Nval
  real*8  :: I_cl(re_eq_n_class), edge_cl(re_eq_n_class), tot_cl(re_eq_n_class)

  I_cl = 0.d0;  edge_cl = 0.d0;  tot_cl = 0.d0

  do i = 1, element_list%n_elements
    do ms = 1, n_gauss
      do mt = 1, n_gauss
        x_g = 0.d0; x_s = 0.d0; x_t = 0.d0
        y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
        eq_g = 0.d0
        do kv = 1, n_vertex_max
          iv = element_list%element(i)%vertex(kv)
          do kf = 1, n_degrees
            x_g = x_g + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
            y_g = y_g + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
            x_s = x_s + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
            x_t = x_t + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
            y_s = y_s + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
            y_t = y_t + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
            eq_g = eq_g + node_list%node(iv)%values(1,kf,var_psi) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          enddo
        enddo
        xjac = x_s*y_t - x_t*y_s
        wst  = wgauss(ms) * wgauss(mt) * abs(xjac)
        do s = 1, re_eq_n_class
          lhat = (re_cl_alpha(s)*x_g - eq_g - re_cl_A_axis(s)) &
                 / (re_cl_A_edge(s) - re_cl_A_axis(s))
          Nval      = re_eq_nprof_at(lhat)
          I_cl(s)   = I_cl(s) - EL_CHG * re_cl_vpar(s) * re_cl_w(s) * Nval / x_g * wst
          ! Beam-edge bookkeeping: how much of the CARRIED current sits in the
          ! outermost 5% of the label range, i.e. how hard the beam edge is
          ! against the loss boundary. This REPLACES the former lost-current
          ! fraction, which accumulated the CLIPPED table -- re_eq_nprof_eval
          ! returns Nprof(1) for lhat > 1 -- over every Gauss point outside the
          ! beam. On a domain that extends past the plasma that sum is
          ! dominated by vacuum volume and is not a current fraction at all
          ! (it reported 0.64 on the JET-like limiter case, whose beam is
          ! entirely confined). With A_edge taken from the loss boundary,
          ! everything on lhat <= 1 is confined by construction: nothing is
          ! lost, so the quantity worth reporting is the edge sharpness.
          tot_cl(s) = tot_cl(s) + abs(Nval) / x_g * wst
          if (lhat .gt. 0.95d0) edge_cl(s) = edge_cl(s) + abs(Nval) / x_g * wst
        enddo
      enddo
    enddo
  enddo

  do s = 1, re_eq_n_class
    if (tot_cl(s) .gt. 0.d0) then
      re_cl_edge_frac(s) = edge_cl(s) / tot_cl(s)
    else
      re_cl_edge_frac(s) = 0.d0
    endif
  enddo
  I_RE = sum(I_cl)
  re_eq_I_now = I_RE

end subroutine re_eq_total_current


!=======================================================================
!> Rescale Nprof so that the RE current equals the prescribed re_eq_I_RE
!> (q_shape mode; called every Picard iteration).
subroutine re_eq_rescale_current(my_id, node_list, element_list)
  use data_structure
  implicit none
  integer,                  intent(in) :: my_id
  type (type_node_list),    intent(in) :: node_list
  type (type_element_list), intent(in) :: element_list
  real*8 :: I_now
  call re_eq_total_current(my_id, node_list, element_list, I_now)
  if (abs(I_now) .gt. 0.d0) then
    re_nprof = re_nprof * (re_eq_I_RE / I_now)
  endif
end subroutine re_eq_rescale_current


!=======================================================================
!> The label map psihat_m(l) for a class with invariant slope alpha: the
!> drift surface with label l crosses the midplane (Z of the class drift
!> axis) at R_out and R_in; the map takes the average
!>   psihat_m(l) = 0.5 * [psihat(R_out(l)) + psihat(R_in(l))].
!> This is a separate, swappable modelling choice. Handles the near-axis
!> degeneracy where both crossings sit on the same side of the magnetic
!> axis (strong drift shift) by falling back to the outboard branch.
subroutine re_eq_label_map(my_id, node_list, element_list, alpha, A_edge_in, &
                           n_lmap, l_values, psihat_m)
  use data_structure
  use equil_info, only: ES
  use mod_model_settings, only: var_psi
  implicit none
  integer,                  intent(in)  :: my_id
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  real*8,                   intent(in)  :: alpha    !< gamma m v_par/e of the class [Wb/m]
  real*8,                   intent(in)  :: A_edge_in !< loss-boundary A of this class
                                                     !< (re_eq_A_edge_bnd; passed in so the
                                                     !< map cannot re-derive it from the old
                                                     !< flux-surface-boundary assumption)
  integer,                  intent(in)  :: n_lmap
  real*8,                   intent(in)  :: l_values(n_lmap)
  real*8,                   intent(out) :: psihat_m(n_lmap)

  integer :: i, k, l, ifail, i_ax, np
  real*8  :: alpha_eff, R_ax, Z_ax, A_ax, A_edge, dpsi
  real*8  :: R_lo, R_hi, dR, dum1, dum2, lv, ph_out, ph_in, R_out_l, R_in_l
  real*8, allocatable :: Rg(:), psig(:), lhatg(:), phg(:)
  ! contour-map (nodal kernel regression) workspace
  real*8, parameter   :: RE_EQ_MAP_BW = 0.03d0   ! label bandwidth of the kernel
  real*8, allocatable :: lhat_n(:), phn_n(:)
  real*8              :: wsum, wnum, dd, wk, phn_raw

  alpha_eff = alpha

  ! --- shared setup: drift axis (label normalization) and the edge value
  call re_eq_find_drift_axis(node_list, element_list, alpha_eff, &
                             ES%R_axis, ES%Z_axis, R_ax, Z_ax, A_ax, ifail)
  if (ifail .ne. 0) then
    write(*,*) 'ERROR: re_eq_label_map: class drift axis not found'
    stop 1
  endif
  A_edge = A_edge_in
  dpsi   = ES%psi_bnd - ES%psi_axis

  ! === Contour-average label map (re_eq_map_mode = 'contour') ================
  ! Average psihat_n over the WHOLE drift surface A_s = const rather than over
  ! its two midplane crossings (the midplane branch below is the 2-point
  ! special case, Eq. 8). Nadaraya-Watson kernel regression of the nodal
  ! psihat_n on the nodal label lhat: uses NO point location (robust like the
  ! drift-axis finder) and represents all poloidal angles of the surface,
  ! which two midplane points do not in shaped geometry. Closed surfaces only
  ! (lhat <= 1.2); psihat_m(l) is forced monotone.
  if (trim(re_eq_map_mode) .eq. 'contour') then
    allocate(lhat_n(node_list%n_nodes), phn_n(node_list%n_nodes))
    do i = 1, node_list%n_nodes
      lhat_n(i) = (alpha_eff*node_list%node(i)%x(1,1,1)                        &
                   - node_list%node(i)%values(1,1,var_psi) - A_ax)            &
                  / (A_edge - A_ax)
      phn_raw   = (node_list%node(i)%values(1,1,var_psi) - ES%psi_axis) / dpsi
      phn_n(i)  = min(max(phn_raw, 0.d0), 1.d0)
      ! diverted grid: exclude open-region nodes (the scrape-off layer and
      ! the private-flux region sit beyond the separatrix, phn_raw > 1) so
      ! they cannot pull the contour average towards the open field lines.
      ! Flagged via the lhat > 1.2 cut used in the binning loop below. Only
      ! when X-point -- a limiter grid has no open region, and gating on
      ! ES%xpoint keeps the non-X-point path bit-identical.
      if (ES%xpoint .and. (phn_raw .gt. 1.d0 + 1.d-6)) lhat_n(i) = 1.d99
    enddo
    do k = 1, n_lmap
      wsum = 0.d0;  wnum = 0.d0
      do i = 1, node_list%n_nodes
        if (lhat_n(i) .gt. 1.2d0) cycle
        dd = (lhat_n(i) - l_values(k)) / RE_EQ_MAP_BW
        wk = exp(-0.5d0*dd*dd)
        wnum = wnum + wk*phn_n(i)
        wsum = wsum + wk
      enddo
      if (wsum .gt. 0.d0) then
        psihat_m(k) = min(max(wnum/wsum, 0.d0), 1.d0)
      else
        psihat_m(k) = min(max(l_values(k), 0.d0), 1.d0)
      endif
    enddo
    do k = 2, n_lmap
      if (psihat_m(k) .lt. psihat_m(k-1)) psihat_m(k) = psihat_m(k-1)
    enddo
    deallocate(lhat_n, phn_n)
    return
  endif

  ! === Midplane-average label map (default, Eq. 8) ===========================
  np = re_eq_n_midplane
  allocate(Rg(np), psig(np), lhatg(np), phg(np))

  ! --- midplane scan at the Z of the effective drift axis. Point location
  !     can fail near the degenerate polar-grid centre; such interior points
  !     are filled by linear interpolation between their valid neighbours
  !     (failed points at the ends of the scan get the boundary psi).
  R_lo = minval(node_list%node(1:node_list%n_nodes)%x(1,1,1))
  R_hi = re_eq_R_edge
  dR   = (R_hi - R_lo) / dble(np + 1)
  dpsi = ES%psi_bnd - ES%psi_axis
  do i = 1, np
    Rg(i) = R_lo + dR * dble(i)
    call re_eq_grad_psi(node_list, element_list, Rg(i), Z_ax, psig(i), dum1, dum2, ifail)
    if (ifail .ne. 0) psig(i) = 1.d99          ! mark for the fill below
  enddo
  do i = 1, np
    if (psig(i) .lt. 1.d98) cycle
    ! nearest valid neighbours left and right
    k = i - 1
    do while ((k .ge. 1) .and. (psig(max(k,1)) .gt. 1.d98))
      k = k - 1
    enddo
    l = i + 1
    do while ((l .le. np) .and. (psig(min(l,np)) .gt. 1.d98))
      l = l + 1
    enddo
    if ((k .ge. 1) .and. (l .le. np)) then
      psig(i) = psig(k) + (psig(l) - psig(k)) * (Rg(i) - Rg(k)) / (Rg(l) - Rg(k))
    else if (k .ge. 1) then
      psig(i) = psig(k)
    else if (l .le. np) then
      psig(i) = psig(l)
    else
      psig(i) = ES%psi_bnd
    endif
  enddo
  do i = 1, np
    lhatg(i) = (alpha_eff*Rg(i) - psig(i) - A_ax) / (A_edge - A_ax)
    phg(i)   = min(max((psig(i) - ES%psi_axis) / dpsi, 0.d0), 1.d0)
  enddo

  ! --- index of the effective drift axis on the scan
  i_ax = 1
  do i = 2, np
    if (abs(Rg(i) - R_ax) .lt. abs(Rg(i_ax) - R_ax)) i_ax = i
  enddo

  do k = 1, n_lmap
    lv = min(max(l_values(k), 1.d-9), 1.d0)

    ! outboard branch: lhat increases from ~0 at the drift axis to 1 at the edge
    R_out_l = Rg(np);  ph_out = phg(np)
    do i = i_ax, np-1
      if ((lhatg(i) - lv) * (lhatg(i+1) - lv) .le. 0.d0) then
        ph_out = phg(i) + (phg(i+1)-phg(i)) * (lv - lhatg(i)) / (lhatg(i+1) - lhatg(i))
        exit
      endif
    enddo

    ! inboard branch; may not reach the label for strong drift shifts
    ph_in = ph_out
    do i = i_ax, 2, -1
      if ((lhatg(i) - lv) * (lhatg(i-1) - lv) .le. 0.d0) then
        ph_in = phg(i) + (phg(i-1)-phg(i)) * (lv - lhatg(i)) / (lhatg(i-1) - lhatg(i))
        exit
      endif
    enddo

    psihat_m(k) = 0.5d0 * (ph_out + ph_in)
  enddo

  deallocate(Rg, psig, lhatg, phg)

end subroutine re_eq_label_map


!=======================================================================
!> Linear response operator K(i,k) = dI(psihat_i)/dN_k : the toroidal
!> current that unit weight of the k-th Nprof basis function puts INSIDE
!> the flux surface psihat_i.
!>
!> This is the object the label map (re_eq_label_map) approximates by a
!> single psihat per label. The source is EXACTLY linear in the Nprof
!> coefficients at fixed psi, so K needs no linearization: it is assembled
!> by binning every Gauss point into (its label interval, its flux level).
!> The edge factor is folded in, so K describes the current actually
!> present in the equilibrium (beam envelope and wall cut included).
!>
!> Why it matters: a drift surface is NOT a flux surface. It spans a RANGE
!> of psihat (~0.18 at 10 MeV) and threads a given flux surface only
!> PARTIALLY. Collapsing that to one psihat per label is what makes the
!> implicit inverse ill-conditioned and what leaves the outermost labels
!> free to absorb an arbitrary amplitude. K represents the partial coverage
!> exactly.
subroutine re_eq_response_operator(node_list, element_list, n_lev, ph_lev, K, Ktot)
  use mod_parameters, only: n_vertex_max, n_degrees
  use data_structure
  use equil_info, only: ES
  use gauss
  use basis_at_gaussian
  use mod_model_settings, only: var_psi
  implicit none
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  integer,                  intent(in)  :: n_lev
  real*8,                   intent(in)  :: ph_lev(n_lev)
  real*8,                   intent(out) :: K(n_lev, re_eq_n_l)
  real*8,                   intent(out) :: Ktot(re_eq_n_l)   !< same integral with NO flux cut:
                                                             !< the TOTAL current per unit N_k

  integer :: i, iv, kv, kf, ms, mt, s, j, kk, jb
  real*8  :: x_g, y_g, x_s, x_t, y_s, y_t, eq_g, xjac, wst
  real*8  :: lhat, lc, dl, frac, coef, phn, dpsi
  real*8  :: D(n_lev+1, re_eq_n_l)

  D    = 0.d0
  dl   = re_nprof_l(2) - re_nprof_l(1)
  dpsi = ES%psi_bnd - ES%psi_axis

  do i = 1, element_list%n_elements
    do ms = 1, n_gauss
      do mt = 1, n_gauss
        x_g = 0.d0; x_s = 0.d0; x_t = 0.d0
        y_g = 0.d0; y_s = 0.d0; y_t = 0.d0
        eq_g = 0.d0
        do kv = 1, n_vertex_max
          iv = element_list%element(i)%vertex(kv)
          do kf = 1, n_degrees
            x_g = x_g + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
            y_g = y_g + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
            x_s = x_s + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
            x_t = x_t + node_list%node(iv)%x(1,kf,1) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
            y_s = y_s + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_s(kv,kf,ms,mt)
            y_t = y_t + node_list%node(iv)%x(1,kf,2) * element_list%element(i)%size(kv,kf) * H_t(kv,kf,ms,mt)
            eq_g = eq_g + node_list%node(iv)%values(1,kf,var_psi) * element_list%element(i)%size(kv,kf) * H(kv,kf,ms,mt)
          enddo
        enddo
        xjac = x_s*y_t - x_t*y_s
        wst  = wgauss(ms) * wgauss(mt) * abs(xjac)

        ! flux bin: this point lies inside every level with ph_lev >= phn
        phn = (eq_g - ES%psi_axis) / dpsi
        jb  = 1
        do j = 1, n_lev
          if (ph_lev(j) .gt. phn) exit
          jb = j + 1
        enddo
        ! jb > n_lev means the point is outside every evaluated surface: it
        ! still contributes to the TOTAL current, so bin it into the overflow
        ! row rather than dropping it.
        do s = 1, re_eq_n_class
          lhat = (re_cl_alpha(s)*x_g - eq_g - re_cl_A_axis(s)) &
                 / (re_cl_A_edge(s) - re_cl_A_axis(s))
          lc   = min(max(lhat, 0.d0), 1.d0)
          coef = -EL_CHG * re_cl_vpar(s) * re_cl_w(s) &
                 * re_eq_edge_factor(lhat) / x_g * wst
          kk   = min(int(lc/dl) + 1, re_eq_n_l - 1)
          frac = lc/dl - dble(kk-1)
          D(jb, kk  ) = D(jb, kk  ) + coef * (1.d0 - frac)
          D(jb, kk+1) = D(jb, kk+1) + coef * frac
        enddo
      enddo
    enddo
  enddo

  do kk = 1, re_eq_n_l
    K(1,kk) = D(1,kk)
    do j = 2, n_lev
      K(j,kk) = K(j-1,kk) + D(j,kk)
    enddo
    Ktot(kk) = K(n_lev,kk) + D(n_lev+1,kk)
  enddo

end subroutine re_eq_response_operator


!=======================================================================
!> Outer update via the response operator (re_eq_transplant = 'operator').
!>
!> q ~ 1/I_enc (exact in the cylindrical limit; the same relation the
!> cumulative transplant and re_eq_init_nprof already assume), so the
!> enclosed current must be rescaled by q/q_t at every level. With
!> I(psihat) = K N exactly, and writing the update multiplicatively as
!> N -> N (1 + u) so positivity is preserved and u is dimensionless:
!>
!>     M u = q/q_t - 1,      M(i,k) = K(i,k) N_k / I_i
!>
!> solved as damped least squares with a second-difference smoothness
!> penalty of weight re_eq_op_lambda. The ill-posedness that the label map
!> hid in its null space becomes one explicit regularization parameter.
!>
!> With re_eq_absorbing_edge the constraint N(1) = 0 (u(n_l) = -1) is
!> imposed INSIDE the solve, so the neighbouring labels relax into a
!> roll-off. That is the difference from the cumulative path, where the
!> same flag can only zero the last table point afterwards and leaves a
!> one-interval cliff (which the finite elements then render as an
!> element-scale ripple of a few percent of the peak current).
!>
!> Solved through the normal equations with LAPACK dgesv: the system is
!> only re_eq_n_l square, and the regularization keeps it well conditioned.
subroutine re_eq_operator_update(node_list, element_list, n_lev, ph_lev, q_lev, c_amp)
  use data_structure
  implicit none
  type (type_node_list),    intent(in) :: node_list
  type (type_element_list), intent(in) :: element_list
  integer,                  intent(in) :: n_lev
  real*8,                   intent(in) :: ph_lev(n_lev), q_lev(n_lev), c_amp

  integer :: i, k, n_good, info, ipiv(re_eq_n_l)
  real*8  :: Kop(n_lev, re_eq_n_l), M(n_lev, re_eq_n_l), Ivec(n_lev), rhs(n_lev)
  real*8  :: Ktot(re_eq_n_l)
  real*8  :: L(re_eq_n_l, re_eq_n_l), AtA(re_eq_n_l, re_eq_n_l), Atb(re_eq_n_l)
  real*8  :: u(re_eq_n_l), qt_at, wcon, umin, umax

  call re_eq_response_operator(node_list, element_list, n_lev, ph_lev, Kop, Ktot)

  Ivec = matmul(Kop, re_nprof(1:re_eq_n_l))
  M    = 0.d0;  rhs = 0.d0;  n_good = 0
  do i = 1, n_lev
    if (abs(Ivec(i)) .le. 0.d0) cycle
    n_good  = n_good + 1
    M(i,:)  = Kop(i,:) * re_nprof(1:re_eq_n_l) / Ivec(i)
    qt_at   = c_amp * re_eq_qt_eval(ph_lev(i))
    rhs(i)  = abs(q_lev(i)) / max(abs(qt_at), 1.d-30) - 1.d0
  enddo
  if (n_good .eq. 0) return

  ! second-difference smoothness operator (one-sided at the ends)
  L = 0.d0
  do k = 2, re_eq_n_l - 1
    L(k,k-1) =  1.d0;  L(k,k) = -2.d0;  L(k,k+1) = 1.d0
  enddo
  L(1,1) = -1.d0;  L(1,2) = 1.d0
  L(re_eq_n_l,re_eq_n_l-1) = 1.d0;  L(re_eq_n_l,re_eq_n_l) = -1.d0

  AtA = matmul(transpose(M), M) + re_eq_op_lambda * matmul(transpose(L), L)
  Atb = matmul(transpose(M), rhs)

  if (re_eq_absorbing_edge) then
    wcon = 1.d3 * max(maxval(abs(M)), 1.d0)
    AtA(re_eq_n_l,re_eq_n_l) = AtA(re_eq_n_l,re_eq_n_l) + wcon*wcon
    Atb(re_eq_n_l)           = Atb(re_eq_n_l)           - wcon*wcon
  endif

  call dgesv(re_eq_n_l, 1, AtA, re_eq_n_l, ipiv, Atb, re_eq_n_l, info)
  if (info .ne. 0) then
    write(*,'(A,I6,A)') ' WARNING: re_eq: operator solve failed (dgesv info = ', &
      info, '); Nprof left unchanged this iteration'
    return
  endif
  u = Atb

  ! same per-iteration bound as the transplant ratio clamp
  umin = 1.d0/re_eq_ratio_clamp - 1.d0
  umax = re_eq_ratio_clamp - 1.d0
  do k = 1, re_eq_n_l
    u(k) = min(max(u(k), umin), umax)
    re_nprof(k) = max(re_nprof(k) * (1.d0 + re_eq_alpha_out * u(k)), 0.d0)
  enddo

end subroutine re_eq_operator_update

!=======================================================================
!> One outer q-matching update. Takes the q profile computed on psihat
!> levels (from determine_q_profile on the converged psi), builds the
!> transplant ratio through the label map and updates Nprof. Sets
!> converged = .true. when max|q/q_t - 1| < re_eq_tol_q over the
!> label-controllable psihat range.
!>
!> Numerically-motivated choices (validated in the M0 prototype, see
!> util/re_equilibrium_prototype/README.md):
!>  - default 'cumulative' transplant: C(l)=int_0^l Nprof dl' is multiplied
!>    by ratio^alpha and differentiated back (cylindrical-exact Newton
!>    direction; the literal 'pointwise' update converges ~15x slower),
!>  - the ratio is clamped to [1/re_eq_ratio_clamp, re_eq_ratio_clamp],
!>  - the log-ratio is smoothed with a [1/4,1/2,1/4] kernel,
!>  - q_now and q_t are evaluated at the SAME clamped psihat argument,
!>  - convergence is measured directly in q space on the diagnostic levels.
subroutine re_eq_outer_update(my_id, node_list, element_list, n_lev, ph_lev, q_lev, &
                              n_inner, converged)
  use data_structure
  implicit none
  integer,                  intent(in)  :: my_id
  type (type_node_list),    intent(in)  :: node_list
  type (type_element_list), intent(in)  :: element_list
  integer,                  intent(in)  :: n_lev
  real*8,                   intent(in)  :: ph_lev(n_lev)  !< psihat_n of the q evaluation
  real*8,                   intent(in)  :: q_lev(n_lev)   !< q at those levels
  integer,                  intent(in)  :: n_inner        !< inner iterations used (for the log)
  logical,                  intent(out) :: converged

  integer :: k, i, s
  logical :: cur_active, improved
  real*8  :: phm(re_eq_n_l), phe, q_at, qt_at, ratio(re_eq_n_l), lr(re_eq_n_l)
  real*8  :: q_acc(re_eq_n_l), qt_acc(re_eq_n_l)
  real*8  :: cw(re_eq_n_class), cw_sum, ph_beam_cl(re_eq_n_class)
  real*8  :: c_amp, num, den, err, err_ctl, I_now, ph_ctl_max, ph_beam, l_eff, f_cur, err_cur
  real*8  :: C(re_eq_n_l), dl, qq

  re_eq_outer_iter = re_eq_outer_iter + 1

  call re_eq_total_current(my_id, node_list, element_list, I_now)

  ! --- per-class label maps, combined as a current-weighted geometric mean.
  !     A single effective-class map mis-models which psihat the update at
  !     label l controls when the class drift shifts differ strongly, and
  !     the outer iteration then oscillates near the edge (prototype
  !     finding, multi-class case with disparate energies).
  !     q and q_t are always evaluated at the SAME clamped argument.
  cw = abs(re_cl_w(1:re_eq_n_class) * re_cl_vpar(1:re_eq_n_class))
  cw_sum = sum(cw);  cw = cw / cw_sum
  lr = 0.d0;  q_acc = 0.d0;  qt_acc = 0.d0
  ph_ctl_max = ph_lev(1)
  do s = 1, re_eq_n_class
    call re_eq_label_map(my_id, node_list, element_list, re_cl_alpha(s), &
                         re_cl_A_edge(s), re_eq_n_l, re_nprof_l, phm)

    ! psihat of the controllable beam interior for this class: labels
    ! beyond re_eq_l_beam carry no current (vacuum annulus), and labels
    ! inside the envelope roll-off are equally uncontrollable (the envelope
    ! forces j -> 0 there). Ratio evaluations are clamped to the
    ! full-current region.
    if (re_eq_l_beam .ge. 1.d0) then
      l_eff = 1.d0
    else
      l_eff = max(re_eq_l_beam - re_eq_l_beam_width, 0.d0)
    endif
    ph_beam = phm(re_eq_n_l)
    do k = 2, re_eq_n_l
      if (re_nprof_l(k) .ge. l_eff) then
        ph_beam = phm(k-1) + (phm(k) - phm(k-1)) &
                  * (l_eff - re_nprof_l(k-1)) / (re_nprof_l(k) - re_nprof_l(k-1))
        exit
      endif
    enddo
    ph_beam_cl(s) = ph_beam

    do k = 1, re_eq_n_l
      phe   = min(max(phm(k), ph_lev(1)), min(ph_lev(n_lev), ph_beam))
      q_at  = re_eq_interp_q(n_lev, ph_lev, q_lev, phe)
      qt_at = re_eq_qt_eval(phe)
      ! transplant on the MAGNITUDE of q: the ratio |q|/|q_t| gives the
      ! correct update direction for either sign convention (q and q_t share
      ! a sign, enforced by the guard below). Clamping the raw signed values
      ! to 1e-30 would destroy a negative-q target (co-B_phi RE pitch).
      lr(k)     = lr(k)     + cw(s) * log(max(abs(q_at), 1.d-30) / max(abs(qt_at), 1.d-30))
      q_acc(k)  = q_acc(k)  + cw(s) * q_at
      qt_acc(k) = qt_acc(k) + cw(s) * qt_at
    enddo
    ph_ctl_max = max(ph_ctl_max, min(ph_beam, ph_lev(n_lev)))
  enddo

  ! --- sign-consistency guard: the equilibrium q and the target q_t must
  !     share a sign. The transplant can only match the MAGNITUDE |q| (the
  !     sign of q is not a free profile parameter). The sign of q is fixed by
  !     sign(q) = sign(j_phi) * sign(F0), i.e. by BOTH the RE current
  !     direction (the pitch xi in the distribution table) AND the toroidal
  !     field (sign of F0) -- so the pitch alone does not determine it. This
  !     guard checks the OUTCOME (sign of the solved q vs the target), so it
  !     is geometry- and F0-sign agnostic; it does not predict the right xi.
  !     On a mismatch no profile shaping can match q_t, so fail early and
  !     clearly instead of grinding to a large, irreducible error.
  if (re_eq_outer_iter .eq. 1) then
    if (sum(q_acc) * sum(qt_acc) .lt. 0.d0) then
      write(*,*) 'ERROR: re_eq: the equilibrium q has the OPPOSITE sign to the'
      write(*,*) '       target q_t. The sign of q is set by sign(current) *'
      write(*,*) '       sign(F0), i.e. by BOTH the RE pitch xi AND the sign of'
      write(*,*) '       F0 -- so the required xi depends on the sign of F0.'
      write(*,*) '       Make the two agree by flipping whichever is physically'
      write(*,*) '       correct for your scenario: the RE pitch xi, the sign of'
      write(*,*) '       F0, or the sign of the q_t table.'
      write(*,'(A,ES12.4,A,ES12.4,A)') '        (current-weighted mean q = ', &
        sum(q_acc), ', mean q_t = ', sum(qt_acc), ')'
      stop 1
    endif
  endif

  ! --- Where the beam edge actually sits in normalized flux. Reported for
  !     EVERY configuration (it used to be gated on l_beam < 1, so the
  !     default wall-limited case never showed it): with a drift-shifted
  !     label the Ahat = l_beam surface spans a RANGE of psihat_n and its
  !     mean can sit far inside the wall -- 0.90 on the 10 MeV hollow-q case
  !     -- which is exactly the flux region the outer loop then cannot
  !     control. ph_ctl_max is the top of the controllable range, so q above
  !     it is an OUTCOME, not a match.
  if (re_eq_outer_iter .eq. 1) then
    write(*,'(A,F7.4,A)') ' re_eq: beam-edge label l_beam = ', re_eq_l_beam, &
      '; per-class beam edge in normalized flux:'
    do s = 1, re_eq_n_class
      write(*,'(A,I4,A,F8.4,A,ES13.5)') '        class ', s, ':  psihat_n = ', &
        ph_beam_cl(s), '   A_edge [Wb] = ', re_cl_A_edge(s)
      if (ph_beam_cl(s) .gt. 0.98d0) &
        write(*,'(A)') '        WARNING: beam edge of this class is very close to the wall'
    enddo
    write(*,'(A,F8.4)') '        q is matched only below psihat_n = ', ph_ctl_max
  endif

  ! --- amplitude factor: 1 in full_q mode; least-squares shape amplitude in
  !     q_shape mode (the current is prescribed there, so only the shape of
  !     q can be matched; the achieved amplitude is reported)
  c_amp = 1.d0
  if (trim(re_eq_match_mode) .eq. 'q_shape') then
    num = sum(q_acc * qt_acc)
    den = sum(qt_acc * qt_acc)
    c_amp = num / den
    lr = lr - log(c_amp)
  endif

  ! Freeze the update factor beyond the control clamp to its value AT the
  ! clamp: labels in the beam-edge roll-off all map to (nearly) the same
  ! frozen psihat, and a kink in the update factor there differentiates
  ! (through the cumulative transplant) into a bump-dip pair in Nprof that
  ! the q-feedback cannot see and that accumulates over the outer
  ! iterations (observed as an overshoot/dip of the current density pinned
  ! at the beam edge, growing with FE resolution).
  if (re_eq_l_beam .lt. 1.d0) then
    l_eff = max(re_eq_l_beam - re_eq_l_beam_width, 0.d0)
    do k = 2, re_eq_n_l
      if (re_nprof_l(k) .gt. l_eff) then
        lr(k:re_eq_n_l) = lr(k-1)
        exit
      endif
    enddo
  endif

  do k = 1, re_eq_n_l
    ratio(k) = min(max(exp(lr(k)), 1.d0/re_eq_ratio_clamp), re_eq_ratio_clamp)
  enddo

  ! --- Convergence metric directly in q space, over the FULL evaluated
  !     psihat range. It used to be restricted to ph_lev(i) <= ph_ctl_max, on
  !     the reasoning that q above the beam-edge label is uncontrollable. That
  !     reasoning does not hold: ph_beam = phm(re_eq_n_l) and phm is forced
  !     monotone, so phm(k) <= ph_beam for every label and the clamp on the
  !     ratio evaluation above NEVER binds when l_beam = 1. The transplant is
  !     therefore already matching q across the whole label range, and
  !     restricting only the METRIC just hid the region where the match is
  !     worst -- reporting e.g. 7.9e-3 for a state whose true worst-case error
  !     was 9.1e-3, i.e. the difference between passing re_eq_tol_q and not.
  !     err_ctl is kept and reported alongside so the split stays visible.
  err     = 0.d0
  err_ctl = 0.d0
  do i = 1, n_lev
    qq  = q_lev(i) / (c_amp * re_eq_qt_eval(ph_lev(i)))
    err = max(err, abs(qq - 1.d0))
    if (ph_lev(i) .le. ph_ctl_max) err_ctl = max(err_ctl, abs(qq - 1.d0))
  enddo
  re_eq_q_err = err

  ! --- current-matching error, only when a target current is actually being
  !     enforced (re_eq_alpha_current > 0). Without it the q error alone is
  !     the criterion, exactly as before.
  cur_active = (re_eq_alpha_current .gt. 0.d0) .and. (re_eq_I_RE .ne. 0.d0)
  err_cur    = 0.d0
  if (cur_active) err_cur = abs(abs(I_now)/abs(re_eq_I_RE) - 1.d0)

  ! Both criteria must be met: a run that matched q but is still far from the
  ! requested current has not finished the job. Reuses re_eq_tol_q rather than
  ! introducing a second tolerance.
  converged   = (err .lt. re_eq_tol_q)
  if (cur_active) converged = converged .and. (err_cur .lt. re_eq_tol_q)

  ! --- The two errors have DIFFERENT roles and must not be conflated:
  !       err     (full range)  -> the VERDICT: convergence and the soft
  !                               tolerance, so re_eq_tol_q means what it says
  !                               and no run can claim success on the strength
  !                               of a region that was never reported.
  !       err_ctl (in the beam) -> the OBJECTIVE: best-iterate tracking,
  !                               stagnation and revert-on-worsening.
  !     Driving the best-iterate/stall logic with the full-range error makes
  !     the solver trade away accuracy where the transplant HAS leverage in
  !     order to chase a residual it cannot move (the flux region only
  !     partially covered by the outermost drift surfaces). Measured on the
  !     mono 60 MeV prototype case: doing so degraded the achievable in-beam
  !     match from 2.4e-3 to 1.4e-2 and made the tail diverge. The full-range
  !     error there floors at ~1.2e-2 whatever the solver does, and the right
  !     response is to REPORT it rather than chase it.
  !
  !     WHY it floors: the map from Nprof(Ahat) to the enclosed current
  !     I(psihat) is a SMOOTHING operator. A drift surface spans a range of
  !     psihat -- |alpha| (R_out - R_in) / |psi_bnd - psi_axis|, about 0.18 at
  !     10 MeV -- and threads a given flux surface only partially, so
  !     structure in the target finer than that width lies outside the
  !     operator's effective range. The width scales with |alpha|, i.e. with
  !     the class ENERGY, which is why 20 MeV matches to 1e-4 and 60 MeV
  !     cannot. It has NOTHING to do with how many classes there are: both of
  !     the cases above are MONO-ENERGETIC, so nothing is being shared and a
  !     per-class Nprof would change nothing. (The multi-class disparity noted
  !     below is a genuine but SEPARATE effect.) The failure is also local:
  !     the gross shape of a target varying over ~0.5 in psihat is matched
  !     fine; it is the EDGE, where only the outermost -- and most smeared --
  !     labels reach, that cannot be resolved.
  ! --- best-iterate tracking and stagnation detection: exactly matching q_t
  !     can be outside the range of the ansatz -- through the smearing above,
  !     and additionally, for a MULTI-class distribution, because one common
  !     Nprof must serve strongly different class maps. Keep the best profile; when no
  !     longer improving, restore it and finish (the caller runs one more
  !     inner Picard pass on the restored profile, in which this routine
  !     only re-evaluates the error and decides on soft acceptance).
  if (re_eq_finishing) then
    ! --- Finishing pass: render the final verdict on this candidate. The
    !     first candidate is the polished best profile; the edge polish is
    !     nearly q-invisible in the circular geometry it was tuned for, but
    !     with a less accurate label map (e.g. shaped plasmas) the structure
    !     it removes can carry real q information, so the polish may WORSEN
    !     the match. In that case revert to the unpolished best profile and
    !     re-converge once more, then accept that -- the polish must never
    !     lose a match the transplant had already achieved.
    if ((.not. re_eq_reverted) .and. (err_ctl .gt. re_eq_best_err)) then
      write(*,'(A,ES10.2,A,ES10.2,A)') &
        ' re_eq: the edge polish worsened the in-beam max|q/q_t-1| (', re_eq_best_err, &
        ' -> ', err_ctl, '); reverting to the unpolished best profile'
      re_nprof(1:re_eq_n_l) = re_eq_best_nprof(1:re_eq_n_l)
      call re_eq_apply_beam_envelope()
      re_eq_reverted = .true.
      write(RE_EQ_LOG_UNIT,'(I6,I8,5ES16.6)') re_eq_outer_iter, n_inner, err, err_ctl, I_now, err_cur, maxval(re_cl_edge_frac)
      call flush_it(RE_EQ_LOG_UNIT)
      return                          ! caller re-converges psi on the best profile
    endif
    if (.not. converged) then
      if ((err .lt. re_eq_tol_q_soft) .and. &
          ((.not. cur_active) .or. (err_cur .lt. re_eq_tol_q_soft))) then
        write(*,'(A)')        ' WARNING: re_eq: q matching stagnated above re_eq_tol_q;'
        write(*,'(A,ES10.2)') '          accepted at the soft tolerance with max|q/q_t-1| = ', err
        write(*,'(A)')        '          (a drift surface spans a RANGE of psihat, so Nprof cannot'
        write(*,'(A)')        '          resolve target structure finer than that width -- worst at the'
        write(*,'(A)')        '          edge, and growing with class energy; see the module header)'
        converged = .true.
        re_eq_soft_accepted = .true.
      endif
    endif
    ! Terminal: the finishing candidate has been judged (accepted or not).
    ! Nprof is frozen in this branch, so further outer iterations would only
    ! re-converge and re-evaluate the identical state -- stop the loop here.
    re_eq_done = .true.
    write(RE_EQ_LOG_UNIT,'(I6,I8,5ES16.6)') re_eq_outer_iter, n_inner, err, err_ctl, I_now, err_cur, maxval(re_cl_edge_frac)
    call flush_it(RE_EQ_LOG_UNIT)
    return
  endif
  ! Best-iterate tracking stays on the q error alone, so the finishing pass
  ! still restores the best q match. Stagnation, however, must not fire while
  ! the CURRENT is still improving: with the total-current control on, the
  ! early phase legitimately trades q error for current progress, and a
  ! q-only counter would abort the run in the middle of it.
  improved = (err_ctl .lt. 0.98d0 * re_eq_best_err)
  if (improved) then
    re_eq_best_err = err_ctl
    re_eq_best_nprof(1:re_eq_n_l) = re_nprof(1:re_eq_n_l)
  endif
  if (cur_active) then
    if (err_cur .lt. 0.98d0 * re_eq_best_err_cur) improved = .true.
    re_eq_best_err_cur = min(re_eq_best_err_cur, err_cur)
  endif
  if (improved) then
    re_eq_n_stall = 0
  else
    re_eq_n_stall = re_eq_n_stall + 1
  endif

  write(*,'(A,I4,A,ES11.3,A,ES13.5,A,ES10.2)') &
    ' re_eq outer ', re_eq_outer_iter, ':  max|q/q_t-1| = ', err, &
    '   I_RE [A] = ', I_now, '   max edge fraction = ', maxval(re_cl_edge_frac)
  write(*,'(A,ES11.3,A,F7.4,A)') &
    '                (within the beam-edge label range: ', err_ctl, &
    ', i.e. psihat_n <= ', ph_ctl_max, ')'
  if (cur_active) &
    write(*,'(A,ES11.3,A,ES12.4,A)') '                |I_RE/target - 1| = ', err_cur, &
      '   (target ', re_eq_I_RE, ' A)'
  if (trim(re_eq_match_mode) .eq. 'q_shape') &
    write(*,'(A,F10.5)') '                q amplitude (achieved/target) = ', c_amp
  if (maxval(re_cl_edge_frac) .gt. 2.d-1) &
    write(*,'(A,ES10.2,A)') ' WARNING: re_eq: ', maxval(re_cl_edge_frac), &
      ' of the current of the worst class is carried on the outermost 5% of'  // &
      ' the label range: the beam edge is hard against the loss boundary'

  write(RE_EQ_LOG_UNIT,'(I6,I8,5ES16.6)') re_eq_outer_iter, n_inner, err, err_ctl, I_now, err_cur, maxval(re_cl_edge_frac)
  call flush_it(RE_EQ_LOG_UNIT)

  if (converged .or. (re_eq_n_stall .ge. 15) .or. (re_eq_outer_iter .ge. re_eq_max_it_out)) then
    ! --- enter the finishing pass: take the best profile, remove the
    !     null-space ripple of Nprof near the controllability edge (nearly
    !     invisible to q, but it imprints element-scale-looking oscillations
    !     on the edge current density), re-converge psi once more, and give
    !     the final verdict on that state
    if (re_eq_best_err .lt. err_ctl) re_nprof(1:re_eq_n_l) = re_eq_best_nprof(1:re_eq_n_l)
    call re_eq_smooth_nprof()          ! includes the beam-edge table hygiene
    write(*,'(A,I4,A,ES10.2)') ' re_eq: finishing after ', re_eq_outer_iter, &
      ' outer iterations (best in-beam max|q/q_t-1| = ', min(re_eq_best_err, err_ctl)
    write(*,'(A)') '        ): applied the edge null-space polish to Nprof;'
    write(*,'(A)') '        final convergence is evaluated on the polished profile'
    re_eq_finishing = .true.
    converged = .false.
    return
  endif

  ! --- smooth the log-ratio (Nprof is smooth; point features in the measured
  !     ratio are q-evaluation artifacts the update must not chase)
  lr = log(ratio)
  do k = 2, re_eq_n_l - 1
    ratio(k) = exp(0.25d0*lr(k-1) + 0.5d0*lr(k) + 0.25d0*lr(k+1))
  enddo

  select case (trim(re_eq_transplant))
  case ('pointwise')
    do k = 1, re_eq_n_l
      re_nprof(k) = re_nprof(k) * ratio(k)**re_eq_alpha_out
    enddo
  case ('cumulative')
    dl = re_nprof_l(2) - re_nprof_l(1)
    C(1) = 0.d0
    do k = 2, re_eq_n_l
      C(k) = C(k-1) + 0.5d0*(re_nprof(k) + re_nprof(k-1)) * dl
    enddo
    do k = 1, re_eq_n_l
      C(k) = C(k) * ratio(k)**re_eq_alpha_out
    enddo
    do k = 2, re_eq_n_l - 1
      re_nprof(k) = (C(k+1) - C(k-1)) / (2.d0*dl)
    enddo
    re_nprof(1)         = (-1.5d0*C(1) + 2.d0*C(2) - 0.5d0*C(3)) / dl
    re_nprof(re_eq_n_l) = ( 1.5d0*C(re_eq_n_l) - 2.d0*C(re_eq_n_l-1) + 0.5d0*C(re_eq_n_l-2)) / dl
    re_nprof = max(re_nprof, 0.d0)
  case ('operator')
    call re_eq_operator_update(node_list, element_list, n_lev, ph_lev, q_lev, c_amp)
  case default
    write(*,*) 'ERROR: unknown re_eq_transplant: ', trim(re_eq_transplant)
    stop 1
  end select

  ! confine the updated profile to the beam (the cumulative update can
  ! regenerate small current beyond the beam edge when differentiating C)
  call re_eq_apply_beam_envelope()

  ! --- total-current control: a SLOW uniform rescale towards re_eq_I_RE,
  !     deliberately outside the q solve (see re_eq_alpha_current). I_now is
  !     the current of the state just evaluated; the same per-iteration bound
  !     as the transplant ratio keeps it from lurching.
  if ((re_eq_alpha_current .gt. 0.d0) .and. (abs(I_now) .gt. 0.d0)) then
    f_cur = (abs(re_eq_I_RE) / abs(I_now))**re_eq_alpha_current
    f_cur = min(max(f_cur, 1.d0/re_eq_ratio_clamp), re_eq_ratio_clamp)
    re_nprof(1:re_eq_n_l) = re_nprof(1:re_eq_n_l) * f_cur
    call re_eq_apply_beam_envelope()
  endif

end subroutine re_eq_outer_update


!=======================================================================
!> Apply the beam-edge envelope to the stored Nprof table: the current is
!> confined to labels below re_eq_l_beam with a C1 smoothstep roll-off of
!> width re_eq_l_beam_width, leaving a current-free (vacuum) annulus
!> between the beam edge and the wall. Applied to the TABLE itself so that
!> the GS source, the current integral, the output file and the marker
!> loader all inherit it consistently. Must be re-applied after every
!> operation that modifies Nprof (transplant update, smoothing): the
!> cumulative transplant and the polish would otherwise regenerate
!> current beyond the beam edge.
subroutine re_eq_apply_beam_envelope()
  ! table hygiene only: zero the RAW table strictly beyond the beam edge
  ! (where the edge factor vanishes anyway) so that the cumulative
  ! transplant carries no phantom current there. Idempotent by
  ! construction -- the smoothstep envelope itself is applied at
  ! EVALUATION time (re_eq_edge_factor), never to the table.
  implicit none
  integer :: k
  ! Absorbing edge: the last confined drift surface carries no density. This
  ! is the one place every profile-modifying path (transplant, edge polish,
  ! initial guess) passes through, so enforcing it here covers them all.
  ! In the cumulative transplant this is a post-hoc zero of the last table
  ! point: the update is explicit and pointwise, so the NEIGHBOURS do not
  ! adapt and the profile meets zero across one table interval rather than
  ! rolling off. Prototype (hollow 10 MeV): the jump at the beam edge falls
  ! from 0.94 to 0.28 of the profile peak. A true roll-off needs the
  ! constraint imposed inside a solve, i.e. the response-operator update.
  if (re_eq_absorbing_edge) re_nprof(re_eq_n_l) = 0.d0
  if (re_eq_l_beam .ge. 1.d0) return
  do k = 1, re_eq_n_l
    if (re_nprof_l(k) .gt. re_eq_l_beam) re_nprof(k) = 0.d0
  enddo
end subroutine re_eq_apply_beam_envelope


!=======================================================================
!> Remove the null-space ripple of Nprof: blended [1/4,1/2,1/4] smoothing,
!> full strength towards l = 1, off below l = 0.5. Near the controllability
!> edge of the transplant, short-wavelength structure in Nprof is nearly
!> invisible to q (the flux-surface average smears it), so iteration noise
!> accumulates there without feedback; it is physically insignificant but
!> imprints oscillations on the edge current density.
subroutine re_eq_smooth_nprof()
  implicit none
  integer :: k, ipass, n_pass
  real*8  :: Ns(re_eq_n_l), wb, w_lo, w_wd
  ! In vacuum-annulus mode (l_beam < 1) the polish is much stronger: the
  ! accumulated null-space structure at the beam edge has wavelengths of
  ! order 0.1 in the label (a bump-dip pair pinned at the control clamp),
  ! and the envelope re-imposes the physical edge shape after smoothing
  ! anyway. In wall-limited mode the milder validated setting is kept,
  ! since Nprof near l = 1 carries real edge current.
  if (re_eq_l_beam .ge. 1.d0) then
    n_pass = 6;   w_lo = 0.5d0;  w_wd = 0.3d0
  else
    ! 50 passes ~ smoothing width 0.05 in the label, matched to the
    ! lambda ~ 0.1 null-space undulation at the beam edge
    n_pass = 50;  w_lo = 0.45d0; w_wd = 0.25d0
  endif
  do ipass = 1, n_pass
    Ns = re_nprof(1:re_eq_n_l)
    do k = 2, re_eq_n_l - 1
      Ns(k) = 0.25d0*re_nprof(k-1) + 0.5d0*re_nprof(k) + 0.25d0*re_nprof(k+1)
    enddo
    do k = 1, re_eq_n_l
      wb = min(max((re_nprof_l(k) - w_lo)/w_wd, 0.d0), 1.d0)
      re_nprof(k) = (1.d0 - wb)*re_nprof(k) + wb*Ns(k)
    enddo
  enddo
  call re_eq_apply_beam_envelope()
end subroutine re_eq_smooth_nprof


!=======================================================================
!> Piecewise-linear interpolation of the q profile in psihat.
function re_eq_interp_q(n, ph, q, x) result(qval)
  implicit none
  integer, intent(in) :: n
  real*8,  intent(in) :: ph(n), q(n), x
  real*8              :: qval, xx
  integer             :: k
  xx = min(max(x, ph(1)), ph(n))
  do k = 2, n
    if (xx .le. ph(k)) exit
  enddo
  k = min(k, n)
  qval = q(k-1) + (q(k) - q(k-1)) * (xx - ph(k-1)) / (ph(k) - ph(k-1))
end function re_eq_interp_q


!=======================================================================
!> Write the per-class data, the Nprof table and a summary to
!> re_equilibrium.dat: everything the kinetic marker initializer needs to
!> sample the stationary density n_s ~ w_s Nprof(Ahat_s)/R. Plain text with
!> a versioned header (kept easily extensible; the reader lives in
!> particles/initialisers/initialisers_RE.f90 and must be kept in sync).
subroutine re_eq_write_output(my_id)
  implicit none
  integer, intent(in) :: my_id
  integer :: iunit, s, k

  iunit = 439
  open(iunit, file='re_equilibrium.dat', action='write', status='replace')
  write(iunit,'(A)') '# JOREK kinetic RE drift-surface equilibrium, format version 1'
  write(iunit,'(A)') '# A_s/e = alpha_s * R - psi;  Ahat_s = (A_s/e - A_axis)/(A_edge - A_axis)'
  write(iunit,'(A)') '# n_s(R,Z) = w_s * Nprof(Ahat_s) / R  [m^-3]; psi in the units of the restart file'
  write(iunit,'(A,I6)')     'n_class ', re_eq_n_class
  write(iunit,'(A,I6)')     'n_l     ', re_eq_n_l
  write(iunit,'(A,ES23.15)') 'I_RE    ', re_eq_I_now
  write(iunit,'(A,ES23.15)') 'q_err   ', re_eq_q_err
  write(iunit,'(A,ES23.15)') 'psi_bnd ', re_eq_psi_bnd
  write(iunit,'(A,ES23.15)') 'taper   ', re_eq_edge_taper
  write(iunit,'(A,ES23.15)') 'l_beam  ', re_eq_l_beam
  write(iunit,'(A,ES23.15)') 'l_beam_w', re_eq_l_beam_width
  write(iunit,'(A,ES23.15)') 'R_edge  ', re_eq_R_edge
  write(iunit,'(A)') '# classes: s  E_kin[eV]  xi  weight  gamma  v_par[m/s]  alpha[Wb/m]  A_axis[Wb]  A_edge[Wb]  R_axis[m]  Z_axis[m]  edge_fraction'
  do s = 1, re_eq_n_class
    write(iunit,'(I5,11ES23.15)') s, re_cl_ekin(s), re_cl_xi(s), re_cl_w(s),   &
      re_cl_gamma(s), re_cl_vpar(s), re_cl_alpha(s), re_cl_A_axis(s),          &
      re_cl_A_edge(s), re_cl_R_axis(s), re_cl_Z_axis(s), re_cl_edge_frac(s)
  enddo
  write(iunit,'(A)') '# nprof: l  Nprof(l) [m^-2]'
  do k = 1, re_eq_n_l
    write(iunit,'(2ES23.15)') re_nprof_l(k), re_nprof(k)
  enddo
  close(iunit)

  write(*,*) ' re_eq: wrote re_equilibrium.dat (per-class data + Nprof table)'

end subroutine re_eq_write_output


!=======================================================================
!> Close the convergence log (call at the end of the equilibrium phase).
subroutine re_eq_finalize(converged)
  implicit none
  logical, intent(in) :: converged
  close(RE_EQ_LOG_UNIT)
  if (.not. converged) then
    write(*,*) 'ERROR: re_eq: q-profile matching did NOT converge: max|q/q_t-1| = ', re_eq_q_err
    write(*,*) '       The residual history is in re_eq_convergence.log.'
    stop 1
  endif
end subroutine re_eq_finalize

end module mod_re_kinetic_equilibrium
