! Unit test of mod_re_kinetic_equilibrium against a synthetic polar grid
! replicating the re_kin_equil_600 100 keV case:
!   R_geo=10, amin=1, F0=3, psi = -0.1*(1-(r/a)^2), psi_bnd=0.
! Runs with -ffpe-trap and bounds checks to catch any NaN/overflow in
! re_eq_init / re_eq_update_labels / re_eq_source / re_eq_source_derivs.
program test_re_eq_unit
  use data_structure
  use equil_info
  use phys_module
  use mod_re_kinetic_equilibrium
  implicit none

  type(type_node_list)    :: node_list
  type(type_element_list) :: element_list
  integer, parameter :: n_ring = 50, n_pol_t = 16
  integer :: i, j, k, n, n_nan
  real*8  :: r, th, R_c, Z_c, psi, S, dS_dpsi, dS_dR, S_min, S_max
  real*8, parameter :: a = 1.d0, R0 = 10.d0, psi_ax = -0.1d0

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

  ! --- synthetic polar node set: centre + rings
  n = 1 + n_ring * n_pol_t
  node_list%n_nodes = n
  allocate(node_list%node(n))
  node_list%node(1)%x(1,:,:) = 0.d0
  node_list%node(1)%x(1,1,1) = R0
  node_list%node(1)%x(1,1,2) = 0.d0
  node_list%node(1)%values = 0.d0
  node_list%node(1)%values(1,1,1) = psi_ax
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
      node_list%node(k)%values(1,1,1) = psi_ax * (1.d0 - (r/a)**2)
    enddo
  enddo
  element_list%n_elements = 0

  ! --- equilibrium state + geometry as in the failing run
  ES%psi_axis = psi_ax
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

  call re_eq_init(0)
  write(*,*) 'init OK'

  call re_eq_update_labels(0, node_list, element_list)
  write(*,*) 'labels OK'

  ! --- evaluate the GS source over the grid and at off-node psi values
  n_nan = 0
  S_min = 1.d99;  S_max = -1.d99
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
    psi = psi_ax * (1.d0 - ((R_c - R0)/a)**2) * 0.97d0
    S = re_eq_source(psi, R_c)
    if (S .ne. S) n_nan = n_nan + 1
    S_min = min(S_min, S);  S_max = max(S_max, S)
  enddo

  write(*,'(A,I6)')     ' NaN count      : ', n_nan
  write(*,'(A,2ES14.5)') ' source min/max : ', S_min, S_max
  if (n_nan .eq. 0) then
    write(*,*) 'UNIT TEST PASSED'
  else
    write(*,*) 'UNIT TEST FAILED'
    stop 1
  endif
end program test_re_eq_unit
