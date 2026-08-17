!> Implements the free boundary equilibrium.
module vacuum_equilibrium
  
  use constants
  use vacuum
  use phys_module, only: resistive_wall
  
  implicit none
  
  
  
  contains
  
  
  
  !> Reads the external fields and the poloidal field coils from the STARWALL output
  subroutine import_external_fields(filename, my_id)
    
    use vacuum_response
    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: my_id
    ! --- Local variables
    integer, parameter   :: filehandle = 60
    integer              :: file_version, n_bnd_elems, n_bnd_nodes, dim(2), err  !n_coils already defined in vacuum module
    integer              :: i_start_coil, i_end_coil                                 !Indices for SW coils 
    character(len=512)   :: comment
    
    if ( sr%n_tor == 0 ) return
    if ( sr%i_tor(1) /= 1 ) return ! external fields not necessary in this case
    
    if (sr%ncoil /= 0 ) starwall_equil_coils = .true.  ! Use STARWALL PF coils if provided
    
    ! --- Decide wheter the coils will be given with COIL_FIELD or STARWALL
    if (.not. starwall_equil_coils) then
    
      if ( my_id == 0 ) then
        
        write(*,*) ''            
        write(*,*) '****************************************'
        write(*,*) '* Using COIL_FIELD program for PF coils *'
        write(*,*) '****************************************'
        write(*,*) ''
        write(*,*) 'WARNING: coil currents should be constant in time. If you want proper coil currents time evolution'
        write(*,*) 'please use STARWALL PF coils'
        write(*,*) ''
        
        ! --- Read data from coil field file (only mpi proc 0 reads the file)
        open(filehandle, file=trim(filename), form='formatted', status='old', action='read', iostat=err)
        if ( err /= 0 ) then
          write(*,*) 'ERROR: Could not open external field file ',trim(filename),'.'
          stop
        end if
        read(filehandle,'(a)') comment
        
        file_version = read_intparam(filehandle, 'file_version')
        if ( file_version > 1 ) then
          write(*,*) 'ERROR: COIL data file version ', file_version, ' is not supported.'
          stop
        end if
        
        n_coils     = read_intparam(filehandle, 'n_coils')
        n_bnd_nodes = read_intparam(filehandle, 'n_bnd_nodes')
        n_bnd_elems = read_intparam(filehandle, 'n_bnd_elems')
        
        dim         = (/ 2*n_bnd_nodes, n_coils /)
        
        call read_array(filehandle, 'B_t', dim, float2d=bext_tan)
        call read_array(filehandle, 'B_n', dim, float2d=bext_nor)
        call read_array(filehandle, 'Psi', dim, float2d=bext_psi)

        ! --- From STARWALL file version 6, the coil currents follow the JOREK sign convention
        ! --- and the fields created by the old COIL_FIELD must be reversed
        if (sr%file_version >= 6) then
          bext_tan = -bext_tan
          bext_nor = -bext_nor
          bext_psi = -bext_psi
        endif
        
        32 format(3x,77('-'))
        33 format(3x,a,i8)
        34 format(3x,a,1x,a)
        35 format(3x,a,es20.12)
        write(*,*)
        write(*,32)
        write(*,33) 'EXTERNAL FIELD INFORMATION:'
        write(*,32)
        write(*,34) 'filename    =', trim(filename)
        write(*,33) 'n_coils     =', n_coils
        write(*,33) 'n_bnd_nodes =', n_bnd_nodes
        write(*,33) 'n_bnd_elems =', n_bnd_elems
        if ( vacuum_debug ) then
          write(*,35) 'sum(bext_tan) =', sum(bext_tan)
          write(*,35) 'sum(bext_nor) =', sum(bext_nor)
          write(*,35) 'sum(bext_psi) =', sum(bext_psi)
        end if
        write(*,32)
        write(*,*)
        
        call check_coil_curr_time_trace_input(n_coils) ! check if the user has introduced non existing coils 
        
        if ( .not. allocated(I_coils) ) then
          allocate( I_coils(n_coils) )
          I_coils(1:n_coils) =  pf_coils(1:n_coils)%current 
          write(*,*) 'I_coils allocated '               
        end if
        
      end if
      
      ! --- Broadcast to other MPI procs.
      call MPI_bcast(n_coils,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_bcast(n_bnd_nodes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_bcast(n_bnd_elems, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_bcast(dim,         2, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      
      if ( my_id /= 0 ) then
        if ( allocated(bext_tan) ) deallocate(bext_tan)
        allocate( bext_tan(dim(1),dim(2)) )
        
        if ( allocated(bext_nor) ) deallocate(bext_nor)
        allocate( bext_nor(dim(1),dim(2)) )
        
        if ( allocated(bext_psi) ) deallocate(bext_psi)
        allocate( bext_psi(dim(1),dim(2)) )        
      end if
      
      call MPI_bcast(bext_tan,        dim(1)*dim(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
      call MPI_bcast(bext_nor,        dim(1)*dim(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
      call MPI_bcast(bext_psi,        dim(1)*dim(2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
      
      n_pf_coils = n_coils
    
    else    !STARWALL coils
      
      if ( my_id == 0 ) then
      
        call check_coil_curr_time_trace_input(sr%ncoil)   ! check if the user has introduced non existing coils
      
        if ( .not. resistive_wall ) then
          write(*,*) 'WARNING: ideal wall with equilibrium with starwall_coils is not ready to use yet'
          stop
        end if
      
        write(*,*) ''            
        write(*,*) '***************************************'
        write(*,*) '* Using STARWALL equilibrium PF coils *'
        write(*,*) '***************************************'
        write(*,*) ''
       
        if ( .not. allocated(I_coils) ) then
          allocate( I_coils(sr%ncoil) )
          I_coils(:)                =  0.d0
          i_start_coil = sr%ind_start_pol_coils
          i_end_coil   = i_start_coil + sr%n_pol_coils - 1
          I_coils(i_start_coil:i_end_coil) =  pf_coils(1:sr%n_pol_coils)%current
          
          i_start_coil = sr%ind_start_rmp_coils
          i_end_coil   = i_start_coil + sr%n_rmp_coils - 1
          I_coils(i_start_coil:i_end_coil) =  rmp_coils(1:sr%n_rmp_coils)%current 
          n_coils                   =  sr%ncoil
          write(*,*) 'I_coils allocated '            
        endif
      endif
      
      n_pf_coils = sr%n_pol_coils
  
    endif   !End choice of STARWALL or COIL_FIELD coils
    
    call MPI_bcast(n_coils,                       1, MPI_INTEGER,          0, MPI_COMM_WORLD, err)
    if ( my_id /= 0 ) then        
      if ( allocated(I_coils) ) deallocate(I_coils)
      allocate( I_coils(n_coils) )
    end if 
    call MPI_bcast(I_coils,                 n_coils, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)
    call MPI_bcast(starwall_equil_coils,          1,          MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
    
  end subroutine import_external_fields
  
  
  
  
  !> Calculates the matrix contribution of the boundary integral to the Grad-Shafranov equation
  !! using the vacuum response from STARWALL
  subroutine vacuum_equil(my_id,node_list, bnd_node_list, bnd_elm_list, psi_axis, psi_bnd, a_mat, rhs_vec)
    
    use mod_parameters
    use data_structure
    use gauss
    use basis_at_gaussian
    use vacuum_response
    use constants
    use mpi_mod
    
    implicit none
    
    ! --- Routine parameters
    integer,                      intent(in) :: my_id
    type (type_node_list),        intent(in) :: node_list
    type (type_bnd_node_list),    intent(in) :: bnd_node_list
    type (type_bnd_element_list), intent(in) :: bnd_elm_list
    real*8,                       intent(in) :: psi_axis
    real*8,                       intent(in) :: psi_bnd
    type(type_SP_MATRIX) :: a_mat
    type(type_RHS) :: rhs_vec    
    
    ! --- Local variables
    type (type_bnd_element) :: bndelem_m
    integer :: m_bndelem, l_vertex, l_dof, l_node, l_dir, l_node_bnd, l_index, ms
    integer :: i_vertex, i_dof, i_node, i_dir, i_node_bnd, i_index, i_resp, i_resp_st
    integer :: j_node_bnd, j_dof, j_node, j_dir, j_index, j_resp, ilarge, n_c
    integer :: i, j, i_start_pf, i_end_pf, offset=1
    real*8  :: size_l, dA, testfunc_l, size_i, basfunc_i
    real*8  :: x(n_gauss), y(n_gauss), x_s(n_gauss), y_s(n_gauss)
    real*8  :: common_prefactor, psi_coil_j, B_tan_coil_i, B_tan_coil_i_loc, psi_0_j
    real*8, allocatable :: potentials_real(:)
    integer :: ierr, step

    call equilibrium_VFB(my_id) 
    
    if (starwall_equil_coils) then
      i_start_pf = sr%ind_start_pol_coils
      i_end_pf   = i_start_pf + sr%n_pol_coils -1
      if ( .not. allocated(wall_curr) )       allocate( wall_curr(n_wall_curr) ) 
      if ( .not. allocated (potentials_real)) allocate(potentials_real(sr%ncoil))      
      wall_curr       = 0.d0
      potentials_real = 0.d0
      potentials_real(i_start_pf:i_end_pf) =  I_coils(i_start_pf:i_end_pf) * mu_zero

      offset = sr%ind_start_coils
      do i = 1, n_wall_curr
        if ( (i>=sr%s_ww_inv%ind_start) .and. (i<=sr%s_ww_inv%ind_end) ) then
          wall_curr(i) = & 
          sum(sr%s_ww_inv%loc_mat(i-my_id*sr%s_ww_inv%step,offset:offset+sr%ncoil -1) *potentials_real(:))
        endif
      end do    
        call MPI_ALLReduce(MPI_IN_PLACE, wall_curr, size(wall_curr),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    endif
    
    ilarge = a_mat%nnz
  
    do m_bndelem = 1, bnd_elm_list%n_bnd_elements
      bndelem_m = bnd_elm_list%bnd_element(m_bndelem)
      
      ! --- Determine the values of R,s and Z,s at the Gaussian points.
      call det_coord_bnd(bndelem_m, node_list, R=x, Z=y, R_S=x_s, Z_S=y_s)

      ! --- Select a test function (weak form equation must hold for every test function)
      do l_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
        do l_dof = 1, 2 ! (loop over node dofs)
          l_node     = bndelem_m%vertex(l_vertex)
          l_dir      = bndelem_m%direction(l_vertex,l_dof)
          l_node_bnd = bndelem_m%bnd_vertex(l_vertex)
          l_index    = node_list%node(l_node)%index(l_dir)
          size_l     = bndelem_m%size(l_vertex,l_dof)
          
          ! --- Integration in s-direction is carried out as Gauss quadrature,
          !       int ds ... = sum_ms wgauss(ms) ...,
          !     where wgauss(ms) denotes the Gaussian weights.
          do ms = 1, n_gauss
            testfunc_l = H1(l_vertex,l_dof,ms) * size_l
            dA         = sqrt(x_s(ms)**2 + y_s(ms)**2) ! Integration factor from definition of dA
            
            ! --- Sum over boundary dofs at which response is calculated
            do i_vertex = 1, 2 ! (loop over nodes in element m_bndelem)
              do i_dof = 1, 2 ! (loop over node dofs)
                i_node     = bndelem_m%vertex(i_vertex)
                i_dir      = bndelem_m%direction(i_vertex,i_dof)
                i_node_bnd = bndelem_m%bnd_vertex(i_vertex)
                i_index    = node_list%node(i_node)%index(i_dir)
                size_i     = bndelem_m%size(i_vertex,i_dof)
                i_resp     = bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof)
                basfunc_i  = H1(i_vertex,i_dof,ms) *size_i
                i_resp_st  = (bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) - 1)*sr%n_tor0 &
                  + bnd_node_list%bnd_node(i_node_bnd)%index_starwall(i_dof)-bnd_node_list%bnd_node(i_node_bnd)%index_starwall(1) + 1

                common_prefactor = wgauss(ms) * dA * testfunc_l * basfunc_i
 
                if (starwall_equil_coils) then
                  B_tan_coil_i     = 0.d0
                  B_tan_coil_i_loc = 0.d0
                  if (i_resp_st>=sr%a_ey%ind_start .AND. i_resp_st<=sr%a_ey%ind_end) then
                    B_tan_coil_i_loc  = - sum (sr%a_ey%loc_mat(i_resp_st-my_id*sr%a_ey%step,:) * wall_curr(:) )
                  endif
                  call MPI_Reduce(B_tan_coil_i_loc, B_tan_coil_i, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
                else
                  B_tan_coil_i         =   sum ( I_coils(:) * bext_tan(i_resp,:) )  
                endif
                
                if (my_id == 0) then
                  rhs_vec%val(l_index) = rhs_vec%val(l_index) + common_prefactor * B_tan_coil_i
                  
                  ! --- Sum over boundary dofs contributing to the response
                  do j_node_bnd = 1, bnd_node_list%n_bnd_nodes ! (loop over boundary nodes)
                    do j_dof = 1, 2 ! (loop over node dofs)
                      j_node     = bnd_node_list%bnd_node(j_node_bnd)%index_jorek
                      j_dir      = bnd_node_list%bnd_node(j_node_bnd)%direction(j_dof)
                      j_index    = node_list%node(j_node)%index(j_dir)
                      j_resp     = bnd_node_list%bnd_node(j_node_bnd)%index_starwall(j_dof)
                     
                      if (.not. starwall_equil_coils) then 
                        psi_coil_j = sum( I_coils(:) * bext_psi(j_resp,:) )
                      else
                        psi_coil_j = 0.d0
                      endif
                     
                      psi_0_j    = node_list%node(j_node)%values(1,j_dir,1)
                      
                      if (j_dir == 1) then   ! shift the flux calue by a global constant, if chosen wisely, it helps a lot with convergence
                        psi_0_j = psi_0_j - psi_offset_freeb 
                      endif
                      
                      ilarge                 = ilarge + 1
                      a_mat%irn(ilarge)  = l_index
                      a_mat%jcn(ilarge)  = j_index
                      a_mat%val(ilarge)    = common_prefactor * response_m_eq(i_resp,j_resp)
                      rhs_vec%val(l_index) = rhs_vec%val(l_index)                               &
                        + common_prefactor * response_m_eq(i_resp,j_resp) * psi_coil_j              &
                        - common_prefactor * response_m_eq(i_resp,j_resp) * psi_0_j
                    end do
                  end do
                end if ! my_id == 0
                
              end do
            end do
            
          end do
          
        end do
      end do
      
    end do
    
  if ( my_id ==0 ) a_mat%nnz = ilarge   ! update the size of the matrix
    
  end subroutine vacuum_equil
  
  
  
  subroutine equilibrium_VFB(my_id)
    use mpi_mod
   implicit none
  
   integer, intent(in) :: my_id

   integer  :: i, ierr
   
   if (my_id == 0) write(*,*) ' vertical_FB = ', vertical_FB
   if (my_id == 0) write(*,*) ' radial_FB = ', radial_FB
   
   if ((my_id == 0) .and. (re_coil_ctl .ne. 0.d0)) &
     write(*,*) ' re_coil_ctl = ', re_coil_ctl

   do i=1, n_pf_coils
     ! Kinetic-RE elongation actuator: an ADDITIVE shaping current, applied to
     ! every PF coil including those carrying no position feedback (those are
     ! set once at allocation and never revisited below). It composes with the
     ! position feedbacks rather than competing for coils -- a coil may serve
     ! only one FEEDBACK channel (see the abort below), but this is a separate
     ! additive term. Guarded on /= 0 so a run without it is bit-identical.
     if (re_coil_ctl .ne. 0.d0) &
       I_coils(i) = pf_coils(i)%current + re_eq_coil_amp(i) * re_coil_ctl
     if( abs(vert_FB_amp(i)) .gt. 1.d-6 ) then
       I_coils(i) =  pf_coils(i)%current * (1 + vert_FB_amp(i) * vertical_FB ) &
                  +  re_eq_coil_amp(i) * re_coil_ctl
       if (my_id == 0) write(*,'(a,I7,a,1es12.4)') 'FB coil ==> I_coil(', i, ') = ', I_coils(i)
     endif
     if( abs(rad_FB_amp(i)) .gt. 1.d-6 ) then
       I_coils(i) =  pf_coils(i)%current * (1 + rad_FB_amp(i) * radial_FB ) &
                  +  re_eq_coil_amp(i) * re_coil_ctl
       if (my_id == 0) write(*,'(a,I7,a,1es12.4)') 'FB coil ==> I_coil(', i, ') = ', I_coils(i)
     endif
     if (( abs(vert_FB_amp(i)) .gt. 1.d-6 ) .and. (abs(rad_FB_amp(i)) .gt. 1.d-6 ))  then
       write(*,*) 'Error: You cannot use the same coil for radial and vertical feedback'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
   enddo
    
  end subroutine equilibrium_VFB




  !**********************************************************************
  !* routines borrowed from EQUAL (WZ)                                  *
  !**********************************************************************
  
  !> Find the PF coil-current combination that best reproduces a SHAPING
  !> perturbation of the boundary flux, i.e. the direction that scales the
  !> variation of psi over the boundary about its mean.
  !>
  !> Why this direction. In fixed boundary, scaling the variation of
  !> Psi_boundary about its mean by mu = 1.05 produced a 10 MeV equilibrium
  !> matching BOTH the target q profile and the prescribed current -- it is the
  !> one perturbation known to supply the missing degree of freedom, with an
  !> effective gain dc_glob/d(mu-1) ~ 0.64. A UNIFORM scale of all coil
  !> currents is a different and far weaker knob: measured on the 10 MeV
  !> free-boundary case, a 40% swing of every coil bought 2.6% of q amplitude
  !> (gain ~0.053, and saturating), so c_glob = 1 sat at a coil scale of ~2.0.
  !> The mean of psi over the boundary is pure gauge -- adding a constant to
  !> psi everywhere changes nothing, which is exactly why psi_offset_freeb is
  !> free to exist -- so it is removed from the target here.
  !>
  !> The solve. Each coil contributes bext_psi(j,i) to the boundary flux at
  !> response DOF j per unit current, so the least-squares problem is
  !>
  !>     min_dI  || A dI - d ||^2 + lambda ||dI||^2,
  !>     A(j,i) = bext_psi(j,i),   d(j) = psi_b(j) - mean(psi_b)
  !>
  !> restricted to the VALUE DOFs (j_dof = 1). The derivative DOFs are
  !> deliberately excluded: they carry different units from the values, and
  !> mixing them in one unweighted least squares would silently weight the two
  !> against each other by whatever the mesh scaling happens to be. For a
  !> smooth field the value pattern determines the derivative pattern anyway.
  !>
  !> Normalization: d is the FULL mean-free flux, i.e. the perturbation for
  !> mu - 1 = 1. The returned dI_shape is therefore "per unit (mu-1)", directly
  !> comparable to the mu = 1.05 experiment (which corresponds to 0.05*dI_shape).
  !>
  !> rel_residual is the number that decides whether this approach can work at
  !> all: it is ||A dI - d|| / ||d||, the fraction of the wanted shaping
  !> perturbation the coil set CANNOT produce. Near 0 means the coils span the
  !> direction; near 1 means they cannot make this shape and no control built
  !> on them will supply the missing degree of freedom.
  subroutine re_coil_shaping_direction(my_id, node_list, bnd_node_list, lambda, &
                                       dI_shape, rel_residual, ifail)
    use data_structure
    implicit none

    integer,                   intent(in)  :: my_id
    type (type_node_list),     intent(in)  :: node_list
    type (type_bnd_node_list), intent(in)  :: bnd_node_list
    real*8,                    intent(in)  :: lambda        !< Tikhonov weight, relative
    real*8,                    intent(out) :: dI_shape(:)   !< size n_coils
    real*8,                    intent(out) :: rel_residual
    integer,                   intent(out) :: ifail

    integer :: n_c, n_p, j_bnd, j_node, j_dir, j_resp, i, k, info
    integer, allocatable :: ipiv(:)
    real*8,  allocatable :: A(:,:), d(:), AtA(:,:), Atd(:), r(:)
    real*8               :: psi_mean, dnorm, rnorm, tr

    ifail = 0
    dI_shape     = 0.d0
    rel_residual = 1.d0
    if (my_id .ne. 0) return

    if (.not. allocated(bext_psi)) then
      write(*,*) 'ERROR: re_coil_shaping_direction: bext_psi is not allocated;'
      write(*,*) '       this needs the STARWALL/coil field file (free boundary).'
      ifail = 1;  return
    endif

    ! With STARWALL-supplied equilibrium coils the solver does NOT use
    ! bext_psi -- vacuum_equil takes the coil field through wall_curr and
    ! sr%a_ey instead, and sets psi_coil_j = 0. Fitting against bext_psi there
    ! would return a direction unrelated to what the solve actually sees, so
    ! refuse rather than answer confidently and wrongly.
    if (starwall_equil_coils) then
      write(*,*) 'ERROR: re_coil_shaping_direction: starwall_equil_coils = .true.,'
      write(*,*) '       so the coil field reaches the boundary through the STARWALL'
      write(*,*) '       response (wall_curr / a_ey), not through bext_psi. This fit'
      write(*,*) '       would not describe the operative coil-to-boundary map.'
      ifail = 5;  return
    endif

    n_c = size(bext_psi, 2)
    n_p = bnd_node_list%n_bnd_nodes
    if (size(dI_shape) .lt. n_c) then
      write(*,*) 'ERROR: re_coil_shaping_direction: dI_shape too small,', &
                 size(dI_shape), ' <', n_c
      ifail = 2;  return
    endif
    if (n_p .lt. n_c) then
      write(*,*) 'WARNING: re_coil_shaping_direction: fewer boundary points (', n_p, &
                 ') than coils (', n_c, '); the fit is underdetermined and only the'
      write(*,*) '         regularization makes it unique.'
    endif

    allocate( A(n_p, n_c), d(n_p), AtA(n_c, n_c), Atd(n_c), r(n_p), ipiv(n_c) )

    ! --- assemble on the VALUE dofs of the boundary nodes
    do j_bnd = 1, n_p
      j_node = bnd_node_list%bnd_node(j_bnd)%index_jorek
      j_dir  = bnd_node_list%bnd_node(j_bnd)%direction(1)
      j_resp = bnd_node_list%bnd_node(j_bnd)%index_starwall(1)
      d(j_bnd)   = node_list%node(j_node)%values(1, j_dir, 1)
      A(j_bnd,:) = bext_psi(j_resp, :)
    enddo

    ! --- strip the gauge (mean) component from BOTH sides: the coils can
    !     produce a uniform flux offset, and it does nothing, so leaving it in
    !     would let the fit spend coil current on an unobservable direction
    psi_mean = sum(d) / dble(n_p)
    d = d - psi_mean
    do i = 1, n_c
      A(:,i) = A(:,i) - sum(A(:,i)) / dble(n_p)
    enddo

    dnorm = sqrt(sum(d*d))
    if (dnorm .le. 0.d0) then
      write(*,*) 'ERROR: re_coil_shaping_direction: the boundary flux has no'
      write(*,*) '       variation about its mean; nothing to match.'
      ifail = 3;  deallocate(A,d,AtA,Atd,r,ipiv);  return
    endif

    ! --- regularized normal equations. n_coils is small (tens), so the cost is
    !     irrelevant; the regularization is here for CONDITIONING -- coil flux
    !     patterns at the boundary are smooth and strongly overlapping, so A
    !     is near rank deficient and an unregularized solve would hand back a
    !     huge cancelling current set that fits the same shape.
    AtA = matmul(transpose(A), A)
    Atd = matmul(transpose(A), d)
    tr  = 0.d0
    do i = 1, n_c
      tr = tr + AtA(i,i)
    enddo
    do i = 1, n_c
      AtA(i,i) = AtA(i,i) + lambda * tr / dble(n_c)
    enddo

    call dgesv(n_c, 1, AtA, n_c, ipiv, Atd, n_c, info)
    if (info .ne. 0) then
      write(*,*) 'ERROR: re_coil_shaping_direction: dgesv failed, info =', info
      ifail = 4;  deallocate(A,d,AtA,Atd,r,ipiv);  return
    endif
    dI_shape(1:n_c) = Atd(1:n_c)

    r     = matmul(A, dI_shape(1:n_c)) - d
    rnorm = sqrt(sum(r*r))
    rel_residual = rnorm / dnorm

    write(*,*)
    write(*,*) '------------------------------------------------------------'
    write(*,*) ' RE coil shaping direction (per unit mu-1)'
    write(*,*) '------------------------------------------------------------'
    write(*,'(A,I6,A,I4)')   '   boundary value dofs = ', n_p, '     coils = ', n_c
    write(*,'(A,ES12.4)')    '   |mean-free psi_b|   = ', dnorm
    write(*,'(A,F9.4)')      '   relative residual   = ', rel_residual
    if (rel_residual .gt. 0.5d0) then
      write(*,*) '   WARNING: the coil set can reproduce less than half of this'
      write(*,*) '            shaping perturbation. A control built on it will be'
      write(*,*) '            weak for the same reason the uniform scale was.'
    endif
    ! The number that decides usability: how big a change of the EXISTING coil
    ! currents does this direction ask for at the amplitude known to work?
    ! mu = 1.05 closed a 3.2% q offset in fixed boundary, and we need ~6.7% at
    ! 10 MeV, so read the last column as roughly half the eventual demand.
    write(*,'(A)') '   coil    dI per unit (mu-1) [A]      at mu=1.05 [A]   as % of present'
    do i = 1, n_c
      if ((i .le. n_pf_coils) .and. (abs(pf_coils(min(i,n_pf_coils))%current) .gt. 1.d-30)) then
        write(*,'(I7,2ES24.6,F16.2)') i, dI_shape(i), 0.05d0*dI_shape(i), &
          100.d0 * 0.05d0*dI_shape(i) / pf_coils(i)%current
      else
        write(*,'(I7,2ES24.6,A16)') i, dI_shape(i), 0.05d0*dI_shape(i), '     n/a'
      endif
    enddo

    open(77, file='re_coil_shaping.txt', status='replace', action='write')
    write(77,'(A)')          '# RE coil shaping direction: dI per unit (mu-1)'
    write(77,'(A,ES14.6)')   '# relative_residual = ', rel_residual
    write(77,'(A,I6)')       '# n_coils = ', n_c
    write(77,'(A)')          '#  coil        dI_per_unit_mu           dI_at_mu_1.05'
    do i = 1, n_c
      write(77,'(I7,2ES24.10)') i, dI_shape(i), 0.05d0*dI_shape(i)
    enddo
    close(77)
    write(*,*) '   written to re_coil_shaping.txt'
    write(*,*) '------------------------------------------------------------'

    deallocate(A, d, AtA, Atd, r, ipiv)

  end subroutine re_coil_shaping_direction


  subroutine pfcoils(R,Z,br,bz,psi)
  
  implicit none
  
  integer :: n_coils, i
  real*8  :: R, Z, br, bz, psi
  real*8,allocatable :: R_coils(:),Z_coils(:), I_coils(:)
  real*8,allocatable :: br_out(:), bz_out(:), psi_out(:), R_out(:), Z_out(:)
  
  !real*8 :: R_coils(4),Z_coils(4), I_coils(4)
  !real*8 :: br_out(4), bz_out(4), psi_out(4), R_out(4), Z_out(4)
  
  stop 'pfcoils programmed stop'
  
  n_coils = 4
  
  allocate(R_coils(n_coils),Z_coils(n_coils),I_coils(n_coils))
  allocate(R_out(n_coils),  Z_out(n_coils),  br_out(n_coils), bz_out(n_coils),psi_out(n_coils))
  
  R_coils = (/  9.d0, 9.d0, 11.d0, 11.d0 /)
  Z_coils = (/ -1.d0, +1.d0,  -1.d0, +1.d0 /) 
  I_coils = (/  1d0, 1.d0, -1.d0, -1.d0 /)
  I_coils = -0.005 * I_coils    
  
  !R_coils = (/  9.d0, 11.d0, 11.d0, 9.d0 /)
  !Z_coils = (/ -1.d0, -1.d0,  1.d0, 1.d0 /) 
  !I_coils = (/  1d0, -0.7d0, -0.7d0, 1.d0 /)
  
  Z_out(1:n_coils) = - Z + Z_coils(1:n_coils) 
  R_out(1:n_coils) = R 
  
  call brbzv(R_out,R_coils,Z_out,br_out,bz_out,n_coils)
  
  call psicalv(R_out,R_coils,Z_out,psi_out,n_coils)
  
  br  = 0.d0
  bz  = 0.d0
  psi = 0.d0
  
  do i=1, n_coils
    br  = br  + br_out(i) * I_coils(i)
    bz  = bz  + bz_out(i) * I_coils(i)
    psi = psi + psi_out(i)* I_coils(i)
  enddo
  
  return
  end subroutine pfcoils
  
  subroutine brbzv(a1,r1,z1,br,bz,n)                                
  !********************************************************************** 
  !**                                                                  ** 
  !**     MAIN PROGRAM:  MHD FITTING CODE                              ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !**     SUBPROGRAM DESCRIPTION:                                      ** 
  !**          psical computes mutual inductance/2/pi between two      ** 
  !**          circular filaments of radii aa1 and r1 and              ** 
  !**          separation of z1, for mks units multiply returned       ** 
  !**          value by 2.0e-07.                                       ** 
  !**                                                                  ** 
  !**     CALLING ARGUMENTS:                                           ** 
  !**       aa1.............first filament radius                      ** 
  !**       r1..............second filament radius                     ** 
  !**       z1..............vertical separation                        ** 
  !**                                                                  ** 
  !**     REFERENCES:                                                  ** 
  !**          (1) f.w. mcclain and b.b. brown, ga technologies        ** 
  !**              report ga-a14490 (1977).                            ** 
  !**                                                                  ** 
  !**     RECORD OF MODIFICATION:                                      ** 
  !**          26/04/83..........first created                         ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !**                                                                  ** 
  !********************************************************************** 
  implicit none
  integer  ::  n
  real*8   ::  a1(1:n),r1(1:n),z1(1:n),br(1:n),bz(1:n)
  
  real*8            :: z2,r2z2,den,sden,xnuma,xnum,x1,xalog,cay,ee,brkt,a2                               
  ! elliptic functions approximation (see Handbook of Mathematical Functions p.591
  real*8, parameter :: ak1=1.38629436112d0, ak2=0.09666344259d0, ak3=0.03590092383d0, ak4=0.03742563713d0, ak5=0.01451196212d0 
  real*8, parameter :: bk1=0.5d0,           bk2=0.12498593597d0, bk3=0.06880248576d0, bk4=0.03328355346d0, bk5=0.00441787012d0
                                                                   
  real*8, parameter :: ae1=0.44325141463d0, ae2=0.06260601220d0, ae3=0.04757383546d0, ae4=0.01736506451d0, &  
                       be1=0.24998368310d0, be2=0.09200180037d0, be3=0.04069697526d0, be4=0.00526449639d0                   
  integer  ::  i                                                                
  	                          
  if (n>0) then
                                                                     
     do i=1,n                        
  
         z2=z1(i)*z1(i)            
         r2z2=z2+r1(i)*r1(i)       
         den=(a1(i)+r1(i))**2+z2   
         sden=sqrt(den)      
         xnuma=(a1(i)-r1(i))**2+z2 
         xnum=max(xnuma,1d-10*den) 
         x1=xnum/den               
                                   
         xalog=-log(x1)            
                                                                       
         cay=((((x1*bk5+bk4)*x1+bk3)*x1+bk2)*x1+bk1)*xalog  &            
             +(((x1*ak5+ak4)*x1+ak3)*x1+ak2)*x1+ak1                     
                                                                       
         ee=(((be4*x1+be3)*x1+be2)*x1+be1)*x1*xalog          &           
            +(((x1*ae4+ae3)*x1+ae2)*x1+ae1)*x1+1.0d0                      
                                                                       
         brkt=ee/xnum                                                   
  !---------------------------------------------------------------------- 
  !--    br  computation                                                -- 
  !---------------------------------------------------------------------- 
         a2=a1(i)*a1(i)                                                 
         br(i)=((a2+r2z2)*brkt-cay)*z1(i)/(r1(i)*sden)                  
  !---------------------------------------------------------------------- 
  !--    bz  computation                                                -- 
  !---------------------------------------------------------------------- 
         bz(i)=((a2-r2z2)*brkt+cay)/sden                                
  
     end do          
  
  endif                                                          
                                                                         
  return                                                            
  end subroutine brbzv
  
  
  subroutine psicalv(a1,r1,z1,psi,n)
  !**********************************************************************
  !**                                                                  **
  !**     MAIN PROGRAM:  MHD FITTING CODE                              **
  !**                                                                  **
  !**                                                                  **
  !**     SUBPROGRAM DESCRIPTION:                                      **
  !**          psical computes mutual inductance/2/pi between two      **
  !**          circular filaments of radii aa1 and r1 and              **
  !**          separation of z1, for mks units multiply returned       **
  !**          value by 2.0e-07.                                       **
  !**                                                                  **
  !**     CALLING ARGUMENTS:                                           **
  !**       aa1.............first filament radius                      **
  !**       r1..............second filament radius                     **
  !**       z1..............vertical separation                        **
  !**                                                                  **
  !**     REFERENCES:                                                  **
  !**          (1) f.w. mcclain and b.b. brown, ga technologies        **
  !**              report ga-a14490 (1977).                            **
  !**                                                                  **
  !**     RECORD OF MODIFICATION:                                      **
  !**          26/04/83..........first created                         **
  !**                                                                  **
  !**                                                                  **
  !**                                                                  **
  !**********************************************************************
  implicit none
  integer  ::  n
  real*8   ::  a1(1:n),r1(1:n),z1(1:n),psi(1:n)
        
  real*8            :: z2,r2z2,xk,den,sden,xnuma,xnum,x1,xalog,cay,ee,a2 
                                
  real*8, parameter :: ak1=1.38629436112d0, ak2=0.09666344259d0, ak3=0.03590092383d0, ak4=0.03742563713d0, ak5=0.01451196212d0 
  real*8, parameter :: bk1=0.5d0,           bk2=0.12498593597d0, bk3=0.06880248576d0, bk4=0.03328355346d0, bk5=0.00441787012d0
                                                                   
  real*8, parameter :: ae1=0.44325141463, ae2=0.06260601220, ae3=0.04757383546, ae4=0.01736506451, &  
                       be1=0.24998368310, be2=0.09200180037, be3=0.04069697526, be4=0.00526449639                   
  integer  :: i
  
  psi=0.d0
  if(n>0) then
     do i=1,n
        z2=z1(i)*z1(i)
        r2z2=z2+r1(i)*r1(i)
        den=(a1(i)+r1(i))**2+z2
        xk=4.d0*a1(i)*r1(i)/den
        sden=sqrt(den)
        xnuma=(a1(i)-r1(i))**2+z2
  !WZ   xnum=max(xnuma,1d-10*den)
        xnum=max(xnuma,1d-20*den)
        x1=xnum/den
  
        xalog=-log(x1)
  
        cay=((((x1*bk5+bk4)*x1+bk3)*x1+bk2)*x1+bk1)*xalog &
            +(((x1*ak5+ak4)*x1+ak3)*x1+ak2)*x1+ak1
  
        ee=(((be4*x1+be3)*x1+be2)*x1+be1)*x1*xalog &
          +(((x1*ae4+ae3)*x1+ae2)*x1+ae1)*x1+1.d0
  !----------------------------------------------------------------------
  !--   psi computation                                                --
  !----------------------------------------------------------------------
        psi(i)= sden*((1.d0-0.5d0*xk)*cay-ee)
     end do
  endif
  
  return
  end subroutine psicalv
  
  
  
    
end module vacuum_equilibrium
