! Minimal stub modules to syntax-check mod_re_kinetic_equilibrium.f90
! outside the full JOREK build (interfaces mirror the real modules).

module constants
  implicit none
  real*8, parameter :: PI             = 3.1415926535897932385d0
  real*8, parameter :: TWOPI          = 2.d0*PI
  real*8, parameter :: MU_ZERO        = 4.d-7*PI
  real*8, parameter :: EL_CHG         = 1.602176565d-19
  real*8, parameter :: MASS_ELECTRON  = 9.10938291d-31
  real*8, parameter :: SPEED_OF_LIGHT = 2.997924580105029d+8
end module constants

module tr_module
  implicit none
  integer, parameter :: CAT_GRID = 2
  interface tr_allocate
    module procedure tr_allocate1d_d
  end interface
  interface tr_deallocate
    module procedure tr_deallocate1d_d
  end interface
contains
  subroutine tr_allocate1d_d(array1d, begin_dim1, end_dim1, var_name, ocategory, pzeroing)
    real*8, dimension(:), allocatable :: array1d
    integer, intent(in) :: begin_dim1, end_dim1
    character*(*), intent(in) :: var_name
    integer, optional, intent(in) :: ocategory
    logical, optional, intent(in) :: pzeroing
    allocate(array1d(begin_dim1:end_dim1))
    array1d = 0.d0
  end subroutine tr_allocate1d_d
  subroutine tr_deallocate1d_d(array1d, var_name, ocategory)
    real*8, dimension(:), allocatable :: array1d
    character*(*), intent(in) :: var_name
    integer, optional, intent(in) :: ocategory
    deallocate(array1d)
  end subroutine tr_deallocate1d_d
end module tr_module

module mod_parameters
  implicit none
  integer, parameter :: n_vertex_max = 4
  integer, parameter :: n_degrees    = 4
  integer, parameter :: n_tor        = 1
end module mod_parameters

module data_structure
  use mod_parameters
  implicit none
  type type_node
    real*8 :: x(1,n_degrees,2)
    real*8 :: values(1,n_degrees,10)
  end type type_node
  type type_node_list
    integer :: n_nodes
    type(type_node), allocatable :: node(:)
  end type type_node_list
  type type_element
    integer :: vertex(n_vertex_max)
    real*8  :: size(n_vertex_max,n_degrees)
  end type type_element
  type type_element_list
    integer :: n_elements
    type(type_element), allocatable :: element(:)
  end type type_element_list
end module data_structure

module phys_module
  implicit none
  real*8 :: F0    = 3.d0
  real*8 :: R_geo = 10.d0
  real*8 :: amin  = 1.d0
end module phys_module

module equil_info
  implicit none
  type type_equil_state
    real*8 :: R_axis, Z_axis, psi_axis, psi_bnd
    real*8 :: psi_xpoint(2), Z_xpoint(2)
  end type type_equil_state
  type(type_equil_state) :: ES
end module equil_info

module mod_model_settings
  implicit none
  integer, parameter :: var_psi = 1
end module mod_model_settings

module gauss
  implicit none
  integer, parameter :: n_gauss = 4
  real*8 :: wgauss(n_gauss)
end module gauss

module basis_at_gaussian
  use mod_parameters
  use gauss
  implicit none
  real*8 :: H(n_vertex_max,n_degrees,n_gauss,n_gauss)
  real*8 :: H_s(n_vertex_max,n_degrees,n_gauss,n_gauss)
  real*8 :: H_t(n_vertex_max,n_degrees,n_gauss,n_gauss)
end module basis_at_gaussian

module mod_interp
  implicit none
  interface interp_PRZ
    module procedure interp_PRZ_1
  end interface
contains
  pure subroutine interp_PRZ_1(node_list, element_list, i_elm, i_v, n_v, s, t, phi, &
                               P, P_s, P_t, P_phi, R, R_s, R_t, Z, Z_s, Z_t, deltas)
    use data_structure
    type (type_node_list),    intent(in)  :: node_list
    type (type_element_list), intent(in)  :: element_list
    integer,                  intent(in)  :: i_elm
    integer,                  intent(in)  :: n_v, i_v(n_v)
    real*8,                   intent(in)  :: s, t, phi
    real*8,                   intent(out) :: P(n_v), P_s(n_v), P_t(n_v), P_phi(n_v)
    real*8,                   intent(out) :: R, R_s, R_t, Z, Z_s, Z_t
    logical, optional, intent(in)         :: deltas
    P = 0.d0; P_s = 0.d0; P_t = 0.d0; P_phi = 0.d0
    R = 0.d0; R_s = 0.d0; R_t = 0.d0; Z = 0.d0; Z_s = 0.d0; Z_t = 0.d0
  end subroutine interp_PRZ_1
end module mod_interp

subroutine find_RZ(node_list, element_list, R_find, Z_find, R_out, Z_out, ielm_out, s_out, t_out, ifail)
  use data_structure
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list
  real*8    :: R_find, Z_find, R_out, Z_out, s_out, t_out
  integer   :: ielm_out, ifail
  R_out = R_find; Z_out = Z_find; ielm_out = 1; s_out = 0.5d0; t_out = 0.5d0; ifail = 0
end subroutine find_RZ

subroutine flush_it(unit)
  integer, intent(in) :: unit
  flush(unit)
end subroutine flush_it
