! Unit test of mod_re_kinetic_equilibrium against a synthetic polar grid
! replicating the re_kin_equil_600 100 keV case:
!   R_geo=10, amin=1, F0=3, psi = psi_ax*(1-(r/a)^2), psi_bnd=0.
! psi_ax > 0 (psi MAXIMUM at the axis) is the only sign consistent with the
! xi = -0.99 distribution below: Nprof >= 0 and v_par < 0 give a negative GS
! source, hence Delta* psi < 0 and psi peaked at the axis. The sign matters
! now that A_edge is taken from the boundary -- it selects which extremum of
! A = alpha R - psi bounds the confined region. Both signs are exercised.
! Runs with -ffpe-trap and bounds checks to catch any NaN/overflow in
! re_eq_init / re_eq_update_labels / re_eq_source / re_eq_source_derivs.
program test_re_eq_unit
  use data_structure
  use equil_info
  use phys_module
  use mod_re_kinetic_equilibrium
  implicit none

  type(type_node_list)      :: node_list
  type(type_element_list)   :: element_list
  type(type_bnd_node_list)  :: bnd_node_list
  integer, parameter :: n_ring = 50, n_pol_t = 16
  integer :: i, j, k, n, n_nan, i_sign
  real*8  :: psi_a
  real*8  :: r, th, R_c, Z_c, psi, S, dS_dpsi, dS_dR, S_min, S_max, A_edge_ref
  real*8, parameter :: a = 1.d0, R0 = 10.d0, psi_ax = 0.1d0

  ! --- input files
  open(20, file='re_distribution.dat', status='replace')
  write(20,'(A)') '1.0e5  -0.99  1.0'
  close(20)
  open(20, file='qprofile_target.dat', status='replace')
  write(20,'(A)') '0.00  1.30'
  write(20,'(A)') '0.10  1.33'
  write(20,'(A)') '0.20  1.42'
  write(20,'(A)') '0.30  1.56'
  write(20,'(A)') '0.40  1.74'
  write(20,'(A)') '0.50  1.95'
  write(20,'(A)') '0.60  2.19'
  write(20,'(A)') '0.70  2.46'
  write(20,'(A)') '0.80  2.76'
  write(20,'(A)') '0.90  3.09'
  write(20,'(A)') '1.00  3.45'
  close(20)

  n = 1 + n_ring * n_pol_t
  allocate(bnd_node_list%bnd_node(n_pol_t))
  n_nan = 0
  S_min = 1.d99;  S_max = -1.d99

  ! ===== both psi branches: psi_ax > 0 is the physical one for xi < 0 (see
  ! ===== the header); psi_ax < 0 is kept as an extra NaN/branch trap.
  do i_sign = 1, 2
  psi_a = merge(psi_ax, -psi_ax, i_sign .eq. 1)
  write(*,'(A,ES11.3)') ' --- synthetic grid with psi_axis = ', psi_a

  ! --- synthetic polar node set: centre + rings
  node_list%n_nodes = n
  if (.not. allocated(node_list%node)) allocate(node_list%node(n))
  node_list%node(1)%x(1,:,:) = 0.d0
  node_list%node(1)%x(1,1,1) = R0
  node_list%node(1)%x(1,1,2) = 0.d0
  node_list%node(1)%values = 0.d0
  node_list%node(1)%values(1,1,1) = psi_a
  k = 1
  do i = 1, n_ring
    r = a * dble(i) / dble(n_ring)
    do j = 1, n_pol_t
      th = 2.d0 * 3.14159265358979d0 * dble(j-1) / dble(n_pol_t)
      k = k + 1
      node_list%node(k)%x = 0.d0
      node_list%node(k)%values = 0.d0
      node_list%node(k)%x(1,1,1) = R0 + r*cos(th)
      node_list%node(k)%x(1,1,2) = r*sin(th)
      node_list%node(k)%values(1,1,1) = psi_a * (1.d0 - (r/a)**2)
    enddo
  enddo
  element_list%n_elements = 0

  ! --- boundary nodes = the outermost ring (psi = psi_bnd there, so this
  !     synthetic grid is the flux-surface special case and the new
  !     loss-boundary A_edge must reproduce the old alpha*maxR - psi_bnd)
  bnd_node_list%n_bnd_nodes = n_pol_t
  do j = 1, n_pol_t
    bnd_node_list%bnd_node(j)%index_jorek = 1 + (n_ring-1)*n_pol_t + j
  enddo

  ! --- equilibrium state + geometry as in the failing run
  ES%psi_axis = psi_a
  ES%psi_bnd  = 0.d0
  ES%R_axis   = R0
  ES%Z_axis   = 0.d0
  F0    = 3.d0
  R_geo = R0
  amin  = a

  ! --- module inputs
  re_kinetic_equilibrium = .true.
  re_eq_dist_file = 're_distribution.dat'
  re_eq_q_file    = 'qprofile_target.dat'

  call re_eq_init(0)          ! guarded: the tables are read once
  write(*,*) 'init OK'

  call re_eq_update_labels(0, node_list, element_list, bnd_node_list)
  write(*,*) 'labels OK'

  ! --- INVARIANCE CHECK for the loss-boundary A_edge. psi = psi_bnd all along
  !     this boundary, i.e. the flux-surface special case, so the boundary
  !     extremum of alpha*R - psi must reproduce the previous
  !     alpha * maxval(R) - psi_bnd EXACTLY -- that is what keeps the existing
  !     circular/limiter cases bit-identical. It holds for the physical sign
  !     pairing only: sign(j_phi) = -sign(v_par) = -sign(alpha) forces psi to
  !     peak at the axis when alpha < 0, and then the bounding boundary point
  !     is the one at maximum R. In the flipped (unphysical) branch the
  !     bounding point is at MINIMUM R and the two definitions legitimately
  !     differ -- which is precisely the error the old formula would make.
  A_edge_ref = re_cl_alpha(1) * maxval(node_list%node(1:n)%x(1,1,1)) - ES%psi_bnd
  write(*,'(A,2ES16.8)') ' A_edge (boundary / old formula) : ', re_cl_A_edge(1), A_edge_ref
  if (i_sign .eq. 1) then
    if (abs(re_cl_A_edge(1) - A_edge_ref) .gt. 1.d-13 * max(abs(A_edge_ref), 1.d-30)) then
      write(*,*) 'UNIT TEST FAILED: A_edge not invariant on a flux-surface boundary'
      stop 1
    endif
    write(*,*) 'A_edge invariance OK (physical branch)'
  else
    A_edge_ref = re_cl_alpha(1) * minval(node_list%node(1:n)%x(1,1,1)) - ES%psi_bnd
    if (abs(re_cl_A_edge(1) - A_edge_ref) .gt. 1.d-13 * max(abs(A_edge_ref), 1.d-30)) then
      write(*,*) 'UNIT TEST FAILED: A_edge is not the boundary extremum'
      stop 1
    endif
    write(*,*) 'A_edge = boundary extremum at MIN R, as expected in this branch'
  endif

  ! --- evaluate the GS source over the grid and at off-node psi values
  do k = 1, n
    R_c = node_list%node(k)%x(1,1,1)
    psi = node_list%node(k)%values(1,1,1)
    S = re_eq_source(psi, R_c)
    if (S .ne. S) n_nan = n_nan + 1
    S_min = min(S_min, S);  S_max = max(S_max, S)
    call re_eq_source_derivs(psi, R_c, S, dS_dpsi, dS_dR)
    if (S .ne. S .or. dS_dpsi .ne. dS_dpsi .or. dS_dR .ne. dS_dR) n_nan = n_nan + 1
  enddo
  ! off-node samples (Gauss-point-like psi/R combinations)
  do k = 1, 1000
    R_c = 9.0d0 + 2.0d0 * dble(k-1) / 999.d0
    psi = psi_a * (1.d0 - ((R_c - R0)/a)**2) * 0.97d0
    S = re_eq_source(psi, R_c)
    if (S .ne. S) n_nan = n_nan + 1
    S_min = min(S_min, S);  S_max = max(S_max, S)
  enddo

  enddo   ! i_sign

  write(*,'(A,I6)')     ' NaN count      : ', n_nan
  write(*,'(A,2ES14.5)') ' source min/max : ', S_min, S_max
  if (n_nan .eq. 0) then
    write(*,*) 'UNIT TEST PASSED'
  else
    write(*,*) 'UNIT TEST FAILED'
    stop 1
  endif
end program test_re_eq_unit
