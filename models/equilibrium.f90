subroutine equilibrium(my_id,node_list,element_list,bnd_node_list,bnd_elm_list,xpoint2,xcase2, nice_q)
!-----------------------------------------------------------------------
! Solve the Grad-Shafranov equation to determine the plasma equilibrium
!   both freeboundary and fixed boundary solutions
!-----------------------------------------------------------------------
use tr_module 
use mod_parameters
use data_structure
use phys_module
use mod_poiss
use mod_iterate2area
use mod_plasma_response
use equil_info
use vacuum
use vacuum_equilibrium, only: re_coil_shaping_direction
use mpi_mod
use mod_interp, only: interp
use mod_F_profile
use mod_re_kinetic_equilibrium
implicit none

          
! --- Routine parameters
integer,                      intent(in)    :: my_id
type (type_node_list),        intent(inout) :: node_list
type (type_element_list),     intent(inout) :: element_list
type (type_bnd_node_list),    intent(inout) :: bnd_node_list
type (type_bnd_element_list), intent(inout) :: bnd_elm_list
logical,                      intent(in)    :: xpoint2
integer,                      intent(in)    :: xcase2
logical,                      intent(in)    :: nice_q

! --- Local variables.
type (type_surface_list) :: surface_list, sep_list
integer    :: ierr, n_iter, iter, i, in, mm, i_elm_axis, i_elm_xpoint(2), i_elm_lim, ifail, i_elm
real*8     :: amplitude, psi, psi_bnd
real*8     :: zn,  dn_dpsi,  dn_dpsi2,  dn_dz,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi2_dz,  dn_dpsi_dz2
real*8     :: zT,  dT_dpsi,  dT_dpsi2,  dT_dz,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi2_dz,  dT_dpsi_dz2
real*8     :: zTi, dTi_dpsi, dTi_dpsi2, dTi_dz, dTi_dz2, dTi_dpsi_dz, dTi_dpsi3, dTi_dpsi2_dz, dTi_dpsi_dz2
real*8     :: zTe, dTe_dpsi, dTe_dpsi2, dTe_dz, dTe_dz2, dTe_dpsi_dz, dTe_dpsi3, dTe_dpsi2_dz, dTe_dpsi_dz2
real*8     :: Ti_prof, Te_prof
real*8     :: zFFprime,dFFprime_dpsi,dFFprime_dz, dFFprime_dpsi_dz, dFFprime_dz2, dFFprime_dpsi2
real*8     :: F_prof, dF_dpsi, dF_dz, dF_dpsi2, dF_dz2, dF_dpsi_dz
real*8     :: xx, x_s, x_t, x_st, x_ss, x_tt, yy, y_s, y_t, y_st, y_ss, y_tt
real*8     :: R_axis, Z_axis, s_axis, t_axis, psi_axis,R, Z, BigR, T0, BigR_s, T0_s
real*8     :: R_lim, Z_lim, s_lim, t_lim, psi_lim, R_out, Z_out, s_out, t_out
real*8     :: R_xpoint(2),Z_xpoint(2),s_xpoint(2),t_xpoint(2), psi_xpoint(2)
real*8     :: zjz, dj_dpsi, dj_dR, dj_dZ, dj_dR_dZ, dj_dR_DR, dj_dZ_dZ, dj_dpsi2, dj_dR_dpsi, dj_dZ_dpsi, psi_n
real*8     :: ps0_s, ps0_t, p_s, p_t, p_ss, p_st, p_tt 
real*8     :: zj0_s, zj0_t, equil_error, equil_value, ps0_x, ps0_y, Z_s, Z_t, xjac, direction, Btot
real*8     :: current_tot, current_int, diff, R_xpoint2(2), Z_xpoint2(2)
real*8     :: sigmas(22), dZ_axis, dR_axis, Z_axis_int, Z_axis_old, R_axis_old, R_axis_int, area_ref
integer    :: n_grids(12)
logical    :: freeboundary_equil2
real*8     :: T_prof, T_0_old, FF_0_old, T_1_old, FF_1_old
real*8, allocatable     :: T_profile(:)
real*8     :: density_prof
real*8, allocatable     :: density_profile(:)
integer    :: nj
real*8     :: rr,ww, drr_dR, drr_dZ, drr_dR2, drr_dZ2, drr_dRdZ

! --- Kinetic RE drift-surface equilibrium (re_kinetic_equilibrium)
integer    :: iter_outer, n_outer_eq, n_outer_fb, i_lev, n_lev_q
logical    :: re_eq_converged
real*8     :: S_re, dS_re_dpsi, dS_re_dR, ph_top
real*8     :: coil_dI_shape(MAX_COILS), coil_shape_res
type (type_surface_list) :: surface_list_q
real*8, allocatable      :: q_lev(:), rad_lev(:), ph_lev(:)

if (my_id .eq. 0) then
  write(*,*) '***************************************'
  write(*,*) '*           equilibrium               *'
  write(*,*) '***************************************'
  write(*,*) '   freeboundary_equil : ',freeboundary_equil
  write(*,*) '   X-point      : ',xpoint2
  write(*,*) '   Xcase        : ',xcase2

  if ((newton_GS_fixbnd .or. newton_GS_freebnd) .and. (use_pastix_eq)) then
    write(*,*) ' '
    write(*,*) ' WARNING: PASTIX 5 IS NOT EFFICIENT FOR THE GRAD-SHAFRANOV SOLVER'
    write(*,*) '           WITH THE NEWTON METHOD. PLEASE USE PASTIX 6,          '
    write(*,*) '           MUMPS OR STRUMPACK INSTEAD. For example               '
    write(*,*) '           (add use_mumps_eq=.t. to namelist and  USE_MUMPS = 1  '
    write(*,*) '           in Makefile.inc)                                      '
    write(*,*) ' '
  endif

  if ((newton_GS_fixbnd .or. newton_GS_freebnd) .and. (.not. xpoint2)) then
    write(*,*) ' '
    write(*,*) ' WARNING: THE NEWTON METHOD FOR THE GRAD-SHAFRANOV SOLVER DOES   '
    write(*,*) '           NOT TAKE EFFECT FOR LIMITER PLASMAS (XPOINT=.F.)      '
    write(*,*) '           AND PICARD ITERATIONS ARE RECOVERED                   '
    write(*,*) '           FURTHER DEVELOPMENTS ARE NEEDED FOR LIMITER PLASMAS   '
    write(*,*) ' '
  endif

endif

freeboundary_equil2 = freeboundary_equil
freeboundary_equil  = .false.

!------------------------------------ fixed boundary equilibrium
n_iter       = 200
psi_bnd      = 0.d0
ES%psi_bnd   = 0.d0
Z_xpoint(1)  = -99.d0
Z_xpoint(2)  = +99.d0
R_xpoint(:)  = R_geo
vertical_FB  = 0.d0
i_elm_xpoint = 0
current_tot  = 0.

! --- Kinetic RE drift-surface equilibrium: nested iteration. The inner
! --- (Picard) loop below solves GS with the per-class invariant labels
! --- refreshed every iteration; the outer loop matches the target q profile
! --- by a transplant update of the common profile function Nprof.
! --- See mod_re_kinetic_equilibrium for the physics and references.
n_outer_eq      = 1
re_eq_converged = .true.
if (re_kinetic_equilibrium) then
  if (newton_GS_fixbnd .or. newton_GS_freebnd) then
    if (my_id == 0) then
      write(*,*) 'ERROR: re_kinetic_equilibrium requires PICARD iterations'
      write(*,*) '       (newton_GS_fixbnd=.f., newton_GS_freebnd=.f.).'
      write(*,*) '       The Newton branches would need the Jacobian of the RE'
      write(*,*) '       source, which is not implemented.'
      write(*,*) '       Fixed-boundary (limiter and diverted) and free-boundary'
      write(*,*) '       Picard equilibria ARE supported.'
    endif
    stop 1
  endif
  if (freeboundary_equil .and. (my_id == 0)) then
    write(*,*) ' re_eq: FREE-BOUNDARY equilibrium. The RE labels are refreshed and'
    write(*,*) '        the prescribed current held inside the free-boundary loop;'
    write(*,*) '        the FF''/p'' current feedback is bypassed (it scales FF_0 and'
    write(*,*) '        T_0, which are zero for a pure-RE equilibrium, so it cannot'
    write(*,*) '        control I_RE). The plasma SIZE is then set by the coils --'
    write(*,*) '        see the doc: q_t and I_RE can only both be matched if the'
    write(*,*) '        external field is free to adjust.'
  endif
  if (xpoint2 .and. (my_id == 0)) then
    write(*,*) ' re_eq: DIVERTED (X-point) equilibrium: labels normalized against'
    write(*,*) '        the separatrix (ES%psi_bnd) and the outboard LCFS radius;'
    write(*,*) '        open-region (scrape-off / private-flux) nodes are excluded,'
    write(*,*) '        and q matching is capped below the separatrix.'
  endif
  n_outer_eq      = re_eq_max_it_out + 1   ! +1: final evaluation pass on the
                                           ! restored best profile after stagnation
  re_eq_converged = .false.
  if (my_id == 0) call re_eq_init(my_id)
endif

if (my_id == 0) then

  do iter_outer = 1, n_outer_eq

  do iter = 1, n_iter


    call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
    call print_equil_state(.true.)

    ! --- refresh the per-class drift-surface labels for the present psi
    ! --- (the per-class drift axes move while psi converges), and hold the
    ! --- prescribed RE current whenever one is requested.
    ! --- Gated on re_eq_I_RE rather than on the match mode: holding the
    ! --- current EXACTLY here removes the Nprof amplitude from the outer
    ! --- optimization entirely, leaving it only the shape to work with. That
    ! --- is why q_shape converges so easily -- a uniform rescale changes no
    ! --- shape, and the shape match changes no current, so the two never
    ! --- compete. In full_q the amplitude is otherwise pulled by BOTH the q
    ! --- amplitude and the current target, which is badly conditioned; with
    ! --- the current pinned here, absolute q is still matched, but through
    ! --- the shape (which moves the LCFS, hence q_a ~ a^2 B0 / I).
    if (re_kinetic_equilibrium) then
      call re_eq_update_labels(my_id, node_list, element_list, bnd_node_list)
      if ((trim(re_eq_match_mode) .eq. 'q_shape') .or. (re_eq_I_RE .ne. 0.d0)) &
        call re_eq_rescale_current(my_id, node_list, element_list)
    endif

    if ((ES%ifail_axis .ne. 0) .and. (iter .le. 5)) then
      call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
      call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
      write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
      ES%R_axis     = R_out;    ES%Z_axis = Z_out;  ES%psi_axis   = psi_axis;   
      ES%s_axis     = s_out;    ES%t_axis = t_out;  ES%i_elm_axis = i_elm;
      ES%ifail_axis = ifail   
    endif
    
    if (xpoint2) then
      if (ES%ifail_xpoint == 0) then ! (otherwise, keep the values of the previous iteration as a reasonable guess)
        ES%psi_bnd  = ES%psi_xpoint(1)
        if( (xcase2 .eq. UPPER_XPOINT) .or. ((xcase2 .eq. DOUBLE_NULL) .and. (abs(ES%psi_xpoint(2)-ES%psi_axis) .lt. abs(ES%psi_xpoint(1)-ES%psi_axis))) ) then
          ES%psi_bnd = ES%psi_xpoint(2)
        endif
        psi_bnd     = ES%psi_bnd
        R_xpoint(1) = ES%R_xpoint(1)
        Z_xpoint(1) = ES%Z_xpoint(1)
        R_xpoint(2) = ES%R_xpoint(2)
        Z_xpoint(2) = ES%Z_xpoint(2)
        if(xcase2 .eq. LOWER_XPOINT) ES%Z_xpoint(2) = +99.d0
        if(xcase2 .eq. UPPER_XPOINT) ES%Z_xpoint(1) = -99.d0
      else
        ES%R_xpoint = R_xpoint
        ES%Z_xpoint = Z_xpoint
        ES%psi_bnd  = psi_bnd
        if (freeboundary_equil) then
          ES%Z_xpoint(1) = -99.d0
          ES%Z_xpoint(2) = +99.d0
        endif
      endif
    else
      ES%psi_bnd = ES%psi_lim
    endif

    if (.not. xpoint) then
      if ( (ES%Z_lim .gt. ES%Z_xpoint(1)) .and. (ES%Z_lim .lt. ES%Z_xpoint(2)) ) then
        if (n_limiter /= 0) then   ! else n_limiter = 0 and psi_bnd is set to 0
          ES%psi_bnd = ES%psi_lim
          write(*,'(A,3f8.3)') ' LIMITER PLASMA ',ES%psi_lim,ES%R_lim,ES%Z_lim
        endif
      endif
    endif
  
    if(xcase2 .eq. LOWER_XPOINT) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND  : ',ES%psi_axis,ES%psi_bnd,ES%Z_xpoint(1),ES%ifail_xpoint
    if(xcase2 .eq. UPPER_XPOINT) write(*,'(A,3es14.6,i3)') ' PSI_AXIS, PSI_BND  : ',ES%psi_axis,ES%psi_bnd,ES%Z_xpoint(2),ES%ifail_xpoint

    write(*,'(A,1f14.8)')                       ' PSI_BND - PSI_AXIS : ', ES%psi_bnd-ES%psi_axis 

    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1

    if ( (my_id == 0) .and. forceSDN .and. iter .gt. 2) then
      if (abs(ES%psi_xpoint(1)-ES%psi_xpoint(2)) .ge. SDN_threshold) then
        ! --- Project psi to enforce up/down symmetry
        call Poisson(0,0,node_list,element_list,bnd_node_list,bnd_elm_list, var_psi,var_psi,1, &
                     0.0,1.0,.true.,xcase,ES%Z_xpoint,.false.,.false.,1)
        call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
      end if
    end if

    diff = 0.d0
    do i=1, node_list%n_nodes
      diff = diff + abs(node_list%node(i)%deltas(1,1,1))
    enddo  
    diff = diff / float(node_list%n_nodes)

    ! Error handling, really is no point continuing if diff is NaN
    if (ISNAN(diff)) then
      write(*,*)'Equilibrium diff is NaN - stop here'
      stop
    end if

    write(*,'(A,I4,A,ES10.3)') ' Iteration ', iter, ': diff=', diff
    
    if ( (iter > 1) .and. (diff < equil_accuracy) ) then
      write(*,'(A,I4,A)') ' Fixed boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if ( iter == n_iter) then
      write(*,'(A,ES10.3)') ' WARNING: Fixed boundary equilibrium not fully converged: diff=', diff
      exit
    end if

  enddo

  ! --- Kinetic RE drift-surface equilibrium: outer q-matching update.
  ! --- Compute q(psihat_n) from the converged psi by flux-surface
  ! --- integration (the standard machinery), then transplant-update Nprof.
  if (re_kinetic_equilibrium) then

    call re_eq_q_transplant(iter)

    if (re_eq_converged) then
      write(*,'(A,I4,A)') ' re_eq: q-profile matching converged after ', iter_outer, ' outer iterations'
    else if (re_eq_done) then
      write(*,'(A)') ' re_eq: q-profile matching stopped without reaching the tolerance'
    endif
  endif

  if (re_eq_converged .or. re_eq_done) exit

  enddo ! iter_outer

  if (re_kinetic_equilibrium) then
    if (allocated(surface_list_q%psi_values)) &
      call tr_deallocate(surface_list_q%psi_values,"surface_list_q%psi_values",CAT_GRID)
    if (allocated(ph_lev)) then
      call tr_deallocate(ph_lev, "ph_lev", CAT_GRID)
      call tr_deallocate(q_lev,  "q_lev",  CAT_GRID)
      call tr_deallocate(rad_lev,"rad_lev",CAT_GRID)
    endif
    ! Hand-off (re_equilibrium.dat) and the final verdict are DEFERRED when a
    ! free-boundary solve follows: that solve moves the plasma boundary, so the
    ! per-class labels (A_axis, A_edge, R_axis) and I_RE all change. Writing
    ! here would describe the fixed-boundary equilibrium while the restart
    ! holds the free-boundary one, and the marker loader would sample on stale
    ! labels. re_eq_finalize also closes the convergence log, which we still
    ! want open through the free-boundary iterations.
    if (.not. freeboundary_equil2) then
      call re_eq_write_output(my_id, node_list, element_list, bnd_node_list)
      call re_eq_finalize(re_eq_converged)
    else
      write(*,'(A)') ' re_eq: fixed-boundary phase done; hand-off deferred until'
      write(*,'(A)') '        after the free-boundary solve (the labels move with'
      write(*,'(A)') '        the boundary).'
      if (.not. re_eq_converged) then
        write(*,*) 'ERROR: re_eq: the fixed-boundary q matching did NOT converge;'
        write(*,*) '       it provides the Nprof the free-boundary solve starts'
        write(*,*) '       from, so there is no point continuing.'
        call re_eq_finalize(re_eq_converged)
      endif
    endif
  endif

end if ! my_id == 0

!--------------------------------------- freeboundary equilibrium
freeboundary_equil = freeboundary_equil2

current_int = 0.d0; Z_axis_int = 0.d0; R_axis_int = 0.d0
 
T_0_old = T_0;  FF_0_old = FF_0;  T_1_old = T_1;  FF_1_old = FF_1

if (freeboundary_equil) then

  if (my_id == 0) then

    write(*,*)
    write(*,*) '------------------------------------------------------'
    write(*,*) '--- Iterative solution of freeboundary equilibrium ---'
    write(*,*) '------------------------------------------------------'
    write(*,*)

    ! Take target delta_psi from fixed boundary if not specified
    if ((delta_psi_GS >= 10000.d0) .and. newton_GS_freebnd) then
      write(*,*) ' '
      write(*,*) ' Taking target delta_psi_GS=psi_bnd-psi_axis from fixed boundary equilibrium'
      write(*,*) ' as it has not been specified in the input file'
      write(*,*) ' '
      delta_psi_GS = ES%psi_bnd - ES%psi_axis
    endif
    
    ! Target current and axis for Picard iterations
    if (current_ref .gt. 1.d20) then    !choose fix bnd equilibrium final current in case of non specification of target current
      call integral_current(node_list,element_list,ES%psi_axis,ES%psi_bnd, xpoint2, xcase2, ES%Z_xpoint, current_ref)
    endif
   
    if (Z_axis_ref .gt. 1.d20) then     !choose fix bnd equilibrium final Zaxis in case of non specification of target Zaxis
      Z_axis_ref = ES%Z_axis
    endif
    
    ! Target poloidal cross section area for limiter plasmas
    if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
      n_limiter = 0    ! Use the full domain to search psibnd enclosing given area
      call area_inside_flux_contour(node_list,element_list, xpoint2, xcase2, ES%psi_bnd, area_ref, ES%R_lim, ES%Z_lim)
      write(*,*) ' The reference area from fixed boundaray is = ', area_ref
    endif
  
  end if ! my_id == 0

  ! === Stage B: outer q-matching loop AROUND the free-boundary Picard =====
  ! The free-boundary solve replaces the fixed-boundary one as the "inner"
  ! solve; the transplant update is otherwise identical, so the same
  ! re_eq_q_transplant runs at the end of each pass.
  n_outer_fb = 1
  if (re_kinetic_equilibrium) then
    n_outer_fb = re_eq_max_it_out + 1   ! +1: final evaluation pass, as in phase 1
    if (my_id == 0) then
      ! Re-arm the outer loop. The fixed-boundary phase has already run to a
      ! verdict, so without this re_eq_done is set and the loop below would
      ! exit on its first pass having done nothing. See re_eq_restart_outer
      ! for why the best-iterate records and the Broyden history must go too.
      call re_eq_restart_outer()
      re_eq_converged = .false.
      write(*,*)
      write(*,'(A,I4,A)') ' re_eq: free-boundary q matching enabled, up to ', &
        re_eq_max_it_out, ' outer iterations around the free-boundary solve.'
      write(*,'(A)')      '        Nprof starts from the converged fixed-boundary profile.'
      if (re_eq_lcfs_a .gt. 0.d0) then
        ! Arm it HERE, not from re_eq_lcfs_a alone: the fixed-boundary phase
        ! above has its boundary frozen, so the LCFS is prescribed and nothing
        ! can move it. Including the size in the convergence test there would
        ! give Phase 1 a criterion it cannot meet.
        re_eq_size_active = .true.
        write(*,*)
        write(*,'(A,F9.5,A)') '        SIZE CONTROL ON: driving the LCFS minor radius to ', &
          re_eq_lcfs_a, ' m'
        write(*,'(A)')       '        NOTE: R_axis_ref from the namelist is only the STARTING value.'
        write(*,'(A,F9.5,A)') '              It is a controlled variable from here on (start ', &
          R_axis_ref, ', clamped to +/-25% of the'
        write(*,'(A)')       '              minor radius about it), and the magnetic axis will'
        write(*,'(A)')       '              settle wherever the drift shift puts it -- which is the'
        write(*,'(A)')       '              point: at fixed LCFS the axis SHOULD move with energy.'
        if (R_axis_ref .le. 0.d0) then
          write(*,*) 'ERROR: re_eq: the size control trims R_axis_ref, but a negative'
          write(*,*) '       R_axis_ref switches the radial feedback off entirely'
          write(*,*) '       (equilibrium.f90: if (R_axis_ref<0) radial_FB=0), so there'
          write(*,*) '       is no actuator. Set R_axis_ref to a sensible starting value.'
          stop 1
        endif
        if (re_eq_lcfs_kappa .gt. 0.d0) then
          write(*,'(A,F9.5)') '        Elongation channel ON: target kappa = ', re_eq_lcfs_kappa
          write(*,'(A)')      '        Actuator: additive coil current along re_eq_coil_amp.'
        else
          write(*,'(A)') '        Elongation is NOT controlled: only the minor radius is held.'
          write(*,'(A)') '        The shape is then free, and it moves with RE energy -- the'
          write(*,'(A)') '        vertical field needed to hold the beam radially grows with'
          write(*,'(A)') '        energy and squeezes the plasma taller (measured at matched a:'
          write(*,'(A)') '        kappa 1.187 at 100 keV against 1.207 at 10 MeV, i.e. 2% of q'
          write(*,'(A)') '        amplitude). Set re_eq_lcfs_kappa to hold it.'
        endif
        if (trim(re_eq_match_mode) .ne. 'q_shape') then
          write(*,'(A)') ' WARNING: re_eq: the size control fixes the LCFS while re_eq_I_RE fixes'
          write(*,'(A)') '          the current, and those two together DETERMINE the q amplitude.'
          write(*,'(A)') '          full_q then asks the transplant to match an amplitude it has no'
          write(*,'(A)') '          freedom left to change. Use re_eq_match_mode = ''q_shape'' and'
          write(*,'(A)') '          read the amplitude as a consistency check.'
        endif
        if (xpoint2) then
          write(*,'(A)') ' WARNING: re_eq: the size control is untested for DIVERTED boundaries.'
          write(*,'(A)') '          It works by compressing the plasma against the inboard'
          write(*,'(A)') '          limiter, so the size and the radial position are one knob;'
          write(*,'(A)') '          with an X-point the separatrix moves with the coils instead'
          write(*,'(A)') '          and the response may be weak or reversed. The secant measures'
          write(*,'(A)') '          it either way and will report if there is no response.'
        endif
      endif
      if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
        write(*,'(A)') ' WARNING: re_eq: freeb_equil_iterate_area forces the plasma area back'
        write(*,'(A)') '          to its FIXED-boundary value every iteration, which removes'
        write(*,'(A)') '          exactly the freedom this loop needs -- the plasma size is'
        write(*,'(A)') '          how the coils reconcile q_t with I_RE. Expect the q'
        write(*,'(A)') '          amplitude to stall. Turn it off for RE free-boundary runs.'
      endif
    endif
  endif

  do iter_outer = 1, n_outer_fb

  do iter=1, n_iter_freeb

    if (my_id == 0) then

      write(*,*)
      write(*,'(1x,a,i5,a)') '>>> ITERATION', iter, ' <<<'
 
      call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
      call print_equil_state(.true.)

      if ((ES%ifail_axis .ne. 0) .and. (iter .le. 5)) then
        call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
        call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
        write(*,'(A,3f10.5)')  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
        ES%R_axis     = R_out;    ES%Z_axis = Z_out;  ES%psi_axis   = psi_axis;   
        ES%s_axis     = s_out;    ES%t_axis = t_out;  ES%i_elm_axis = i_elm;
        ES%ifail_axis = ifail   
      endif
      
      write(10,'(i6,9e20.12)') iter, current_tot, ES%R_axis, ES%Z_axis, ES%psi_bnd-ES%psi_axis
      
      ES%psi_bnd = 0.d0
   
      if (xpoint2) then
        if (ES%ifail_xpoint .ne. 1) then      
          ES%psi_bnd  = ES%psi_xpoint(1)
          if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (abs(ES%psi_xpoint(2)-ES%psi_axis) .lt. abs(ES%psi_xpoint(1)-ES%psi_axis))) ) then
            ES%psi_bnd = ES%psi_xpoint(2)
          endif
          if(xcase2 .eq. LOWER_XPOINT) ES%Z_xpoint(2) = +99.d0
          if(xcase2 .eq. UPPER_XPOINT) ES%Z_xpoint(1) = -99.d0
        else
          ES%Z_xpoint(1) = -99.d0 
          ES%Z_xpoint(2) = +99.d0
        endif
      endif
  
      if (.not. xpoint2) then
        if ( (ES%Z_lim .gt. ES%Z_xpoint(1)) .and. (ES%Z_lim .lt. ES%Z_xpoint(2)) ) then
          call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
          if (ES%axis_is_psi_minimum) then
            ES%psi_bnd = min(ES%psi_lim,ES%psi_bnd)
          else
            ES%psi_bnd = max(ES%psi_lim,ES%psi_bnd)
          endif
          write(*,'(A,4f8.3)') ' LIMITER PLASMA ',ES%psi_lim, ES%psi_bnd, ES%R_lim,ES%Z_lim
        endif
      endif

      if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
        call iterate2area(node_list,element_list, ES%psi_axis, ES%psi_lim, xpoint2, xcase2, area_ref, ES%psi_bnd)
      endif
      
      write(*,'(A,1f8.3)') ' Psi_bnd = ', ES%psi_bnd   

      ! --- Kinetic RE drift-surface equilibrium: the per-class labels must be
      ! --- refreshed here for the same reason as in the fixed-boundary Picard
      ! --- (the drift axes and the loss-boundary A_edge move as psi and the
      ! --- plasma boundary converge), and the prescribed current held exactly
      ! --- when one is requested. Analogue of the fixed-boundary block above.
      if (re_kinetic_equilibrium) then
        call re_eq_update_labels(my_id, node_list, element_list, bnd_node_list)
        if ((trim(re_eq_match_mode) .eq. 'q_shape') .or. (re_eq_I_RE .ne. 0.d0)) &
          call re_eq_rescale_current(my_id, node_list, element_list)
      endif

      ! Calculate current feedback
      call integral_current(node_list,element_list,ES%psi_axis, ES%psi_bnd, xpoint2, xcase2, ES%Z_xpoint, current_tot)
  
      current_int = current_int + (current_tot-current_ref)
      
      ! The feedback below controls the total current by scaling FF' and p'.
      ! For a pure-RE equilibrium FF_0 = FF_1 = 0 and T_0 ~ 0, so it is a NO-OP
      ! that cannot control I_RE -- the current comes from Nprof and is held by
      ! re_eq_rescale_current above. Freeze the factor rather than let it wind
      ! up on a current error it has no authority over.
      if (re_kinetic_equilibrium) then
        current_FB_fact = 1.d0
        current_int     = 0.d0
      else if ((mod(iter,n_feedback_current) .eq. 0) .and. (.not. newton_GS_freebnd)) then
        current_FB_fact  = current_FB_fact * (1. - FB_Ip_position * (current_tot-current_ref)/current_ref &
                                                 - FB_Ip_integral *  current_int/current_ref   )
      else if ( cte_current_FB_fact > -1.d90 ) then
        current_FB_fact  = cte_current_FB_fact
      endif
      
      !-------------- Multiplying FF' and p' profiles by the same factor to scale total current -------------------------
      FF_0 = FF_0_old * current_FB_fact   
      FF_1 = FF_1_old * current_FB_fact      
        
      T_0  = T_0_old  * current_FB_fact    
      T_1  = T_1_old  * current_FB_fact
      !------------------------------------------------------------------------------------------------------------------
      
      write(*,'(A,1e12.4)') 'Current Feedback factor = ',  current_FB_fact
      
      !Vertical feedback - needed for vertically unstable plasmas        
      Z_axis_int = Z_axis_int + (ES%Z_axis - Z_axis_ref)
      R_axis_int = R_axis_int + (ES%R_axis - R_axis_ref)
      if (iter .eq. 1) then
        dZ_axis = 0.d0
        dR_axis = 0.d0
      else
        dZ_axis = ES%Z_axis - Z_axis_old
        dR_axis = ES%R_axis - R_axis_old
      end if

    
      if ((mod(iter,n_feedback_vertical) .eq. 0) .and. (iter .ge. start_VFB) .and. (.not. newton_GS_freebnd) ) then
        vertical_FB = FB_Zaxis_position   * (ES%Z_axis-Z_axis_ref) &   ! vertical_FB is used in vacuum_equilibrium.f90 to modify the coils current
                    + FB_Zaxis_integral   * Z_axis_int          &   
                    + FB_Zaxis_derivative * dZ_axis
        radial_FB = FB_Zaxis_position   * (ES%R_axis-R_axis_ref) &   ! radial_FB is used in vacuum_equilibrium.f90 to modify the coils current
                    + FB_Zaxis_integral   * R_axis_int          &   
                    + FB_Zaxis_derivative * dR_axis

      endif
        
      Z_axis_old = ES%Z_axis
      R_axis_old = ES%R_axis
       
    end if ! my_id == 0
    
    if (R_axis_ref<0) radial_FB=0.d0
    call MPI_bcast(vertical_FB, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
    call MPI_bcast(radial_FB, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    ! --- Iterate equation
    call poisson(my_id,-1,node_list,element_list,bnd_node_list,bnd_elm_list,3,1,1, &
                 ES%psi_axis,ES%psi_bnd,xpoint2,xcase2,ES%Z_xpoint,freeboundary_equil,refinement,iter)   !----------- for GS use -1
  
  !  call boundary_check
   
    if (my_id == 0) then
      diff = 0.d0
      do i=1, node_list%n_nodes
        diff = diff + abs(node_list%node(i)%deltas(1,1,1))
      enddo  
      diff = diff / float(node_list%n_nodes)
    
      write(*,'(A,i5,e14.6)') ' iteration, diff : ',iter,diff
    end if ! my_id == 0
  
    call MPI_bcast(diff, 1, MPI_DOUBLE_PRECISION,  0, MPI_COMM_WORLD,ierr)
  
    if ( (iter > 1) .and. (diff < equil_accuracy_freeb) ) then
      if (my_id == 0) write(*,'(A,I4,A)') ' Free boundary equilibrium converged: after', iter, ' iterations'
      exit
    else if (iter == n_iter_freeb) then
      if (my_id == 0) write(*,'(A,ES10.3)') ' WARNING: Free boundary equilibrium not fully converged: diff=', diff
      exit
    end if
  
  enddo ! iter (free-boundary Picard)

  ! --- outer q-matching update on the converged free-boundary equilibrium.
  !     Rank 0 owns Nprof and the GS assembly, so only the VERDICT has to be
  !     broadcast -- the loop below encloses poisson, whose vacuum_equil call
  !     is collective, so every rank must leave it on the same iteration.
  if (re_kinetic_equilibrium) then
    if (my_id == 0) then
      call re_eq_q_transplant(iter)
      ! Stage C: the radial setpoint carries the plasma SIZE, the transplant
      ! carries the q shape, re_eq_rescale_current carries the current -- three
      ! controls on disjoint subspaces. Updated AFTER the transplant so it acts
      ! on the LCFS just measured, and only while the loop is still running:
      ! once the verdict is in, Nprof is frozen and moving the boundary would
      ! invalidate it.
      !
      ! NOT during the finishing pass. That pass freezes Nprof and exists only
      ! to judge whether the edge polish helped, by comparing the q error
      ! before and after; moving the boundary underneath it makes the size
      ! control's effect look like the polish's. Observed exactly that on the
      ! first 10 MeV run: one good size step (a 0.71711 -> 0.70516) was
      ! attributed to the polish, which then reverted BOTH, and the run
      ! finished at the size it started with. The size control must have
      ! finished its work before the endgame begins -- which it will, now that
      ! the size error is part of the convergence test.
      if (re_eq_size_active .and. (.not. re_eq_converged) &
          .and. (.not. re_eq_done) .and. (.not. re_eq_finishing)) &
        call re_eq_lcfs_update(R_axis_ref, re_coil_ctl)
      if (re_eq_converged) then
        write(*,'(A,I4,A)') ' re_eq: free-boundary q-profile matching converged after ', &
          iter_outer, ' outer iterations'
      else if (re_eq_done) then
        write(*,'(A)') ' re_eq: free-boundary q-profile matching stopped without reaching the tolerance'
      endif
    endif
    call MPI_bcast(re_eq_converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(re_eq_done,      1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(R_axis_ref,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_bcast(re_coil_ctl,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (re_eq_converged .or. re_eq_done) exit
  else
    exit                        ! no q matching: one free-boundary solve only
  endif

  enddo ! iter_outer (free-boundary q matching)

  if (freeb_equil_iterate_area .and. (.not. xpoint2)) then
    n_limiter = 1  ! set found limiter (defined inside iterate2area)
  endif

else
  
  psi_offset_freeb = 0.d0
  
endif

if (my_id == 0) then
  ! Update psi axis and boundary with new values from the last iteration of equilibrium solvers 
  call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
        
  write(10,'(i6,9e20.12)') iter, current_tot, R_axis, Z_axis, psi_bnd-psi_axis
  
  if ((ifail .ne. 0) .and. (iter .le. 5)) then
    call find_RZ(node_list,element_list,R_geo,Z_geo,R_out,Z_out,i_elm,s_out,t_out,ifail)
    call interp(node_list,element_list,i_elm,1,1,s_out,t_out,psi_axis,P_s,P_t,P_st,P_ss,P_tt)
    write(*,*)  ' changed magnetic axis to :  ', R_out,Z_out,psi_axis
  endif
  
  psi_bnd = 0.d0
  
  if (xpoint2) then
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint,Z_xpoint,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    if (ifail .ne. 1) then      
      psi_bnd  = psi_xpoint(1)
      if( (xcase2 .eq. 2) .or. ((xcase2 .eq. 3) .and. (abs(psi_xpoint(2)-psi_axis) .lt. abs(psi_xpoint(1)-psi_axis))) ) then
        psi_bnd = psi_xpoint(2)
      endif
      if(xcase2 .eq. LOWER_XPOINT) Z_xpoint(2) = +99.d0
      if(xcase2 .eq. UPPER_XPOINT) Z_xpoint(1) = -99.d0
    else
      Z_xpoint(1) = -99.d0 
      Z_xpoint(2) = +99.d0
    endif
  endif
  if (.not. xpoint2) then
    call find_limiter(my_id,node_list,element_list,bnd_elm_list,psi_lim,R_lim,Z_lim)
    if ( (Z_lim .gt. Z_xpoint(1)) .and. (Z_lim .lt. Z_xpoint(2)) ) then
      call is_axis_psi_mininum(node_list, element_list, bnd_elm_list)
      if (ES%axis_is_psi_minimum) then
        psi_bnd = min(psi_lim,psi_bnd)
      else
        psi_bnd = max(psi_lim,psi_bnd)
      endif
      write(*,'(A,4f8.3)') ' LIMITER PLASMA ',psi_lim, psi_bnd, R_lim,Z_lim
    endif
  endif

  if (freeboundary_equil .and. freeb_equil_iterate_area .and. (.not. xpoint2)) then
    call iterate2area(node_list,element_list, psi_axis, psi_lim, xpoint2, xcase2, area_ref, psi_bnd)
  endif
  
  !------------------------------- end of equilibrium, start filling data
  psi_axis = psi_axis - psi_offset_freeb
  psi_bnd  = psi_bnd  - psi_offset_freeb

  ! --- Kinetic RE: the labels must follow that shift. A = alpha R - psi, so
  !     psi -> psi - offset maps A -> A + offset exactly; without this the zj
  !     fill below (which shifts each node's psi and then evaluates the source)
  !     and the marker hand-off would both combine the new psi with the old
  !     labels, landing at Ahat + offset/(A_edge - A_axis) -- clipped past the
  !     beam edge, so no RE current at all. The shift is zero in the
  !     fixed-boundary branch, which is why this only appears in free boundary.
  ! --- Kinetic RE: which coil combination reproduces a SHAPING perturbation of
  !     the boundary flux, i.e. the mu-direction that was shown in fixed
  !     boundary to make q_t and I_RE simultaneously reachable. Diagnostic
  !     only -- it changes no state, just measures whether this coil set spans
  !     that direction, which is what the uniform scale demonstrably does not
  !     (40% of every coil bought 2.6% of q amplitude at 10 MeV). Computed on
  !     the CONVERGED boundary flux, so it must sit after the solve.
  if (re_kinetic_equilibrium .and. freeboundary_equil) then
    call re_coil_shaping_direction(my_id, node_list, bnd_node_list, 1.d-6, &
                                   coil_dI_shape, coil_shape_res, ifail)
  endif

  if (re_kinetic_equilibrium) then
    if (psi_offset_freeb .ne. 0.d0) call re_eq_shift_labels(psi_offset_freeb)
    ! Hand-off written HERE: psi is now in its final form, the same one the
    ! restart carries, so re_equilibrium.dat and the restart agree.
    if (freeboundary_equil) then
      call re_eq_write_output(my_id, node_list, element_list, bnd_node_list)
      call re_eq_finalize(re_eq_converged)
      ! q-evaluation workspace, re-allocated by re_eq_q_transplant during the
      ! free-boundary outer loop after the fixed-boundary phase released it
      if (allocated(surface_list_q%psi_values)) &
        call tr_deallocate(surface_list_q%psi_values,"surface_list_q%psi_values",CAT_GRID)
      if (allocated(ph_lev)) then
        call tr_deallocate(ph_lev, "ph_lev", CAT_GRID)
        call tr_deallocate(q_lev,  "q_lev",  CAT_GRID)
        call tr_deallocate(rad_lev,"rad_lev",CAT_GRID)
      endif
    endif
  endif

  ! --- This fills in the data for the current variable "zj" (for R-MHD only)
#ifndef fullmhd

  do i=1,node_list%n_nodes
    node_list%node(i)%values(1,1,1) = node_list%node(i)%values(1,1,1) - psi_offset_freeb
    psi = node_list%node(i)%values(1,1,1)
    R   = node_list%node(i)%x(1,1,1)
    Z   = node_list%node(i)%x(1,1,2)
  
    call density(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
                                                               dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
  
    if (with_TiTe) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      zT  	= zTi	       + zTe
      dT_dpsi	= dTi_dpsi     + dTe_dpsi
      dT_dpsi2	= dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3	= dTi_dpsi3    + dTe_dpsi3
      dT_dz	= dTi_dz       + dTe_dz
      dT_dz2	= dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
    		     zT,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    endif
  
    call FFprime(    xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd,zFFprime,dFFprime_dpsi,dFFprime_dz, &
                                                               dFFprime_dpsi2,dFFprime_dz2, dFFprime_dpsi_dz, .true.)

  
    zjz     = zFFprime      - R*R *      (dn_dpsi    * zT + zn * dT_dpsi)

    dj_dpsi = dFFprime_dpsi - R*R *      (dn_dpsi2   * zT + zn * dT_dpsi2  + 2.d0 * dn_dpsi * dT_dpsi)

    dj_dR   =               - 2.d0 * R * (dn_dpsi    * zT + zn * dT_dpsi)

    dj_dZ   = dFFprime_dz   - R*R *      (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)

    ! --- Kinetic RE drift-surface equilibrium: the RE current is part of the
    ! --- GS source and must appear in the current variable zj = Delta*psi as
    ! --- well (second derivatives of the piecewise-linear Nprof vanish a.e.
    ! --- and are left out of the 4th degree of freedom)
    if (re_kinetic_equilibrium) then
      call re_eq_source_derivs(psi, R, S_re, dS_re_dpsi, dS_re_dR)
      zjz     = zjz     + S_re
      dj_dpsi = dj_dpsi + dS_re_dpsi
      dj_dR   = dj_dR   + dS_re_dR
    endif
  
    dj_dR_dR = - 2.d0     * (dn_dpsi     * zT + zn * dT_dpsi)
  
    dj_dZ_dZ = dFFprime_dz2   - R*R * ( dn_dpsi_dz2 * zT   + dn_dpsi_dz * dT_dz  + dn_dz * dT_dpsi_dz  + dn_dz2 * dT_dpsi &
                                      +  dn_dpsi_dz  * dT_dz + dn_dpsi    * dT_dz2 + zn    * dT_dpsi_dz2 + dn_dz  * dT_dpsi_dz)
  
    dj_dpsi2 = dFFprime_dpsi2 - R*R * (dn_dpsi3 * zT + 3.d0 * dn_dpsi * dT_dpsi2 + 3.d0 * dn_dpsi2 * dT_dpsi + zn * dT_dpsi3 )
  
    dj_dR_dZ   = - 2.d0 * R * (dn_dpsi_dz * zT + dn_dpsi * dT_dz + zn * dT_dpsi_dz + dn_dz * dT_dpsi)
  
    dj_dR_dpsi = - 2.d0 * R * (dn_dpsi2   * zT + zn * dT_dpsi2   + 2.d0 * dn_dpsi * dT_dpsi)
  
    dj_dZ_dpsi = dFFprime_dpsi_dz - R*R * ( dn_dpsi2_dz * zT    + dn_dz * dT_dpsi2     + 2.d0 * dn_dpsi_dz * dT_dpsi  &
                                            + dn_dpsi2    * dT_dz + zn    * dT_dpsi2_dz  + 2.d0 * dn_dpsi    * dT_dpsi_dz)
  
  
    node_list%node(i)%values(1,1,3) = zjz
  
    node_list%node(i)%values(1,2,3) = dj_dpsi * node_list%node(i)%values(1,2,1) &
                                    + dj_dR   * node_list%node(i)%x(1,2,1)        &
                                    + dj_dZ   * node_list%node(i)%x(1,2,2)
  
    node_list%node(i)%values(1,3,3) = dj_dpsi * node_list%node(i)%values(1,3,1) &
                                    + dj_dR   * node_list%node(i)%x(1,3,1)        &
                                    + dj_dZ   * node_list%node(i)%x(1,3,2)
  
    node_list%node(i)%values(1,4,3) = dj_dpsi  * node_list%node(i)%values(1,4,1) &
                                    + dj_dR    * node_list%node(i)%x(1,4,1)        &
                                    + dj_dZ    * node_list%node(i)%x(1,4,2)        &
                                    + dj_dR_dR * node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,1)  &
                                    + dj_dZ_dZ * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)  &
                                    + dj_dpsi2 * node_list%node(i)%values(1,2,1) * node_list%node(i)%values(1,3,1)  &
                                    + dj_dR_dZ * ( node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,2)          &
                                                 + node_list%node(i)%x(1,3,1) * node_list%node(i)%x(1,2,2) )        &
                                    + dj_dR_dpsi*( node_list%node(i)%x(1,2,1) * node_list%node(i)%values(1,3,1)   &
                                                 + node_list%node(i)%x(1,3,1) * node_list%node(i)%values(1,2,1) ) &
                                    + dj_dZ_dpsi*( node_list%node(i)%x(1,2,2) * node_list%node(i)%values(1,3,1)   &
                                                 + node_list%node(i)%x(1,3,2) * node_list%node(i)%values(1,2,1) )

    ! --- Add contribution of current ropes
    if ((.not. restart) .and. (n_jropes .ne. 0)) then
      do nj=1,n_jropes
        rr = sqrt((R-R_jropes(nj))**2 + (Z-Z_jropes(nj))**2)
        drr_dR   = (R-R_jropes(nj)) / rr
        drr_dZ   = (Z-Z_jropes(nj)) / rr
        drr_dR2  = 1./rr - (R-R_jropes(nj)) / rr**2 * drr_dR
        drr_dZ2  = 1./rr - (Z-Z_jropes(nj)) / rr**2 * drr_dZ
        drr_dRdZ = - (R-R_jropes(nj)) / rr**2 * drr_dZ
        ww = w_jropes(nj)
        zjz        = 0.d0
        dj_dR      = 0.d0
        dj_dZ      = 0.d0
        dj_dR_dR   = 0.d0
        dj_dZ_dZ   = 0.d0
        dj_dR_dZ   = 0.d0
        if (rr .le. ww) then
          zjz        = current_jropes(nj) * (1.0 - (rr/ww)**2 )**2 * R
          dj_dR      = -4. * current_jropes(nj) * rr / ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR * R + zjz / R
          dj_dZ      = -4. * current_jropes(nj) * rr / ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ * R
          dj_dR_dR   = - zjz/R**2 + dj_dR/R - 4.  * current_jropes(nj) * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR    & 
                                            - 4.*R* current_jropes(nj)        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR**2 & 
                                            - 4.*R* current_jropes(nj) * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR2   & 
                                            + 8.*R* current_jropes(nj) * rr**2/ww**4                       * drr_dR**2 
          dj_dZ_dZ   = - 4.*R* current_jropes(nj)        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ**2 & 
                       - 4.*R* current_jropes(nj) * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dZ2   & 
                       + 8.*R* current_jropes(nj) * rr**2/ww**4                       * drr_dZ**2 
          dj_dR_dZ   = dj_dZ / R                                                                  &
                       - 4.*R* current_jropes(nj)        /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dR*drr_dZ & 
                       - 4.*R* current_jropes(nj) * rr   /ww**2 * (1.0 - (rr/ww)**2 ) * drr_dRdZ      & 
                       + 8.*R* current_jropes(nj) * rr**2/ww**4                       * drr_dR*drr_dZ 
        endif
       
        node_list%node(i)%values(1,1,3) = node_list%node(i)%values(1,1,3) + zjz
       
        node_list%node(i)%values(1,2,3) = node_list%node(i)%values(1,2,3)      &
                                        + dj_dR   * node_list%node(i)%x(1,2,1) &
                                        + dj_dZ   * node_list%node(i)%x(1,2,2)
       
        node_list%node(i)%values(1,3,3) = node_list%node(i)%values(1,3,3)      &
                                        + dj_dR   * node_list%node(i)%x(1,3,1) &
                                        + dj_dZ   * node_list%node(i)%x(1,3,2)
       
        node_list%node(i)%values(1,4,3) = node_list%node(i)%values(1,4,3)       &
                                        + dj_dR    * node_list%node(i)%x(1,4,1) &
                                        + dj_dZ    * node_list%node(i)%x(1,4,2) &
                                        + dj_dR_dR * node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,1)    &
                                        + dj_dZ_dZ * node_list%node(i)%x(1,2,2) * node_list%node(i)%x(1,3,2)    &
                                        + dj_dR_dZ * ( node_list%node(i)%x(1,2,1) * node_list%node(i)%x(1,3,2)  &
                                                     + node_list%node(i)%x(1,3,1) * node_list%node(i)%x(1,2,2) )
      enddo
    endif

  enddo
  
  ! --- Variable projection is better at higher order...
  ! --- (by the way, we could use this for n_order=3 and remove all the above as well, 
  ! --- and remove all derivatives from profiles functions, which are not really needed, 
  ! --- except dn_dpsi and dT_dpsi for current profile...)
  if (n_order .ge. 5) then
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    if (xpoint2) then
      ES%xpoint = xpoint2
      ES%Z_xpoint = Z_xpoint
    endif
    ES%psi_bnd  = psi_bnd
    ES%psi_axis = psi_axis
    ES%Z_xpoint = Z_xpoint
    ES%xpoint   = xpoint
    ES%xcase    = xcase
    call Poisson(my_id,0,node_list,element_list,bnd_node_list,bnd_elm_list, &
                 var_psi,var_zj,1, psi_axis,psi_bnd,xpoint2,xcase2,Z_xpoint,freeboundary_equil,refinement,1)
  endif
#endif
  ! --- END of filling data for current variable "zj" (R-MHD only)
  
  ! --- Find flux surfaces and plot them; determine the q-profile.  
  if (xpoint2 .and. (n_flux .gt. 1)) then
    
    call find_axis(my_id,node_list,element_list,psi_axis,R_axis,Z_axis,i_elm_axis,s_axis,t_axis,ifail)
    call find_xpoint(my_id,node_list,element_list,psi_xpoint,R_xpoint2,Z_xpoint2,i_elm_xpoint,s_xpoint,t_xpoint,xcase2,ifail)
    
    n_grids = 0
    sigmas  = 0.d0
    
    ! Build up some arrays to send as routine parameters to define_flux_values
    sigmas(1)  = SIG_closed(1); sigmas(2)  = SIG_theta
    sigmas(3)  = SIG_open    ; sigmas(4)  = SIG_outer   ; sigmas(5)  = SIG_inner
    sigmas(6)  = SIG_private ; sigmas(7)  = SIG_up_priv
    sigmas(8)  = SIG_leg_0   ; sigmas(9)  = SIG_leg_1
    sigmas(10) = SIG_up_leg_0; sigmas(11) = SIG_up_leg_1
    sigmas(12) = dPSI_open   ; sigmas(13) = dPSI_outer  ; sigmas(14) = dPSI_inner
    sigmas(15) = dPSI_private; sigmas(16) = dPSI_up_priv
    sigmas(17) = SIG_theta_up
    sigmas(18) = SIG_closed(2); sigmas(19) = SIG_closed(3)
    sigmas(20) = xr_closed(1) ; sigmas(21) = xr_closed(2)
    sigmas(22) = xr_closed(3)
  
    n_grids(1) = 2*n_flux   ; n_grids(2) = n_tht
    n_grids(3) = 2*n_open   ; n_grids(4) = 2*n_outer  ; n_grids(5) = 2*n_inner
    n_grids(6) = 2*n_private; n_grids(7) = 2*n_up_priv
    n_grids(8) = n_leg      ; n_grids(9) = n_up_leg
    if (xcase .eq. LOWER_XPOINT) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(7) = 0
      n_grids(9) = 0
    endif
    if (xcase .eq. UPPER_XPOINT) then
      n_grids(4) = 0
      n_grids(5) = 0
      n_grids(6) = 0
      n_grids(8) = 0
    endif
  
    ! Allocate surface_list structure (that's for plotting only)
    if (xcase2 .eq. LOWER_XPOINT) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_private
    if (xcase2 .eq. UPPER_XPOINT) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_up_priv
    if (xcase2 .eq. DOUBLE_NULL ) surface_list%n_psi = 2*n_flux + 2*n_open + 2*n_outer + 2*n_inner + 2*n_private + 2*n_up_priv
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    ! Allocate sep_list structure (that's for plotting only)  
    sep_list%n_psi =3
    if(xcase .eq. DOUBLE_NULL) sep_list%n_psi =6
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    
    ! Define the flux values to be plotted...
    psi_axis = psi_axis+0.01 !Just offset a little, because finding surfaces along the side of an element (on the xpoint grid) can be hard...
    call define_flux_values(node_list, element_list, surface_list, sep_list, xcase2, psi_xpoint, n_grids, sigmas)
    psi_axis = psi_axis-0.01 !Put it back, it's not used anyway, but just for principle!
    
  else
    surface_list%n_psi = 200  
    if (allocated(surface_list%psi_values)) call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
    call tr_allocate(surface_list%psi_values,1,surface_list%n_psi,"surface_list%psi_values",CAT_GRID)
    
    do i = 1, surface_list%n_psi
      surface_list%psi_values(i) = 1.25d0*(float(i)/float(surface_list%n_psi))**2 * (psi_bnd - psi_axis) + psi_axis
    enddo
    
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,surface_list)
  
    sep_list%n_psi =1
    if (allocated(sep_list%psi_values)) call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
    call tr_allocate(sep_list%psi_values,1,sep_list%n_psi,"sep_list%psi_values",CAT_GRID)
    sep_list%psi_values(1) = psi_bnd
  
    call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,sep_list)
  endif
  
  if (freeboundary_equil) then
    !call plot_coils(.true.)
    call plot_flux_surfaces(node_list,element_list,surface_list,.false.,4,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
  
    call plot_flux_surfaces(node_list,element_list,surface_list,.true.,4,xpoint2,xcase2)
    call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    !call plot_coils(.false.)
  else
    if (xpoint2 .and. (n_flux .gt. 1)) then
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,xpoint2,xcase2)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    else
      call plot_flux_surfaces(node_list,element_list,surface_list,.true.,1,.false.,0)
      call plot_flux_surfaces(node_list,element_list,sep_list,.false.,1,xpoint2,xcase2)
    endif
  endif
  
  if (nice_q) then
    if (xpoint2) then
      ES%xpoint = xpoint2
      ES%Z_xpoint = Z_xpoint
    endif
    ES%psi_bnd  = psi_bnd
    ES%psi_axis = psi_axis
    ES%Z_xpoint = Z_xpoint
    ES%xpoint   = xpoint
    ES%xcase    = xcase
    call q_profile(node_list,element_list,surface_list,psi_axis,psi_bnd,psi_xpoint,Z_xpoint)
  endif
  
  !================ Temperature and density profiles =f(psi_norm) similar to q(psi_norm) needed to calculate neoclassical coef===========
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  call tr_allocate(T_profile,1,surface_list%n_psi,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  call tr_allocate(density_profile,1,surface_list%n_psi,"density_profile",CAT_GRID)
  
  do i=2,surface_list%n_psi
    psi= surface_list%psi_values(i)
    
    if (with_TiTe) then
      call temperature_i(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
           Ti_prof,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)
  
      call temperature_e(xpoint2, xcase2, Z, Z_xpoint, psi,psi_axis,psi_bnd, &
           Te_prof,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)
      T_prof      = Ti_prof      + Te_prof
      dT_dpsi     = dTi_dpsi     + dTe_dpsi
      dT_dpsi2    = dTi_dpsi2    + dTe_dpsi2
      dT_dpsi3    = dTi_dpsi3    + dTe_dpsi3
      dT_dz       = dTi_dz       + dTe_dz
      dT_dz2      = dTi_dz2      + dTe_dz2
      dT_dpsi_dz  = dTi_dpsi_dz  + dTe_dpsi_dz
      dT_dpsi2_dz = dTi_dpsi2_dz + dTe_dpsi2_dz
      dT_dpsi_dz2 = dTi_dpsi_dz2 + dTe_dpsi_dz2 
    else
      call temperature(.false.,xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,T_prof,dT_dpsi,dT_dz,dT_dpsi2,dT_dz2,             &
        dT_dpsi_dz,dT_dpsi3,dT_dpsi_dz2, dT_dpsi2_dz)
    endif

    call density( .false., xcase2,0., Z_xpoint, psi,psi_axis,psi_bnd,density_prof,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,             &
          dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)
  
    T_profile(i)=T_prof
    density_profile(i)=density_prof
  end do
  
  write(*,*) '***************************************'
  write(*,*) 'output T and rho profiles (in JOREK units) for neoclassical profile calculation'
  ! --- Write out T and rho profiles to "T_rho_profiles.dat".
  open(432, file='T_rho_profiles.dat', action='write', status='replace')
  do i=2, surface_list%n_psi
     write(432,'(3ES13.5)') T_profile(i), density_profile(i)
  end do
  close(432)
  !========================= end modif ===========================================
  
  if (allocated(surface_list%psi_values))    call tr_deallocate(surface_list%psi_values,"surface_list%psi_values",CAT_GRID)
  if (allocated(surface_list%flux_surfaces)) deallocate(surface_list%flux_surfaces)
  if (allocated(sep_list%psi_values))        call tr_deallocate(sep_list%psi_values,"sep_list%psi_values",CAT_GRID)
  if (allocated(sep_list%flux_surfaces))     deallocate(sep_list%flux_surfaces)
  
  if (allocated(T_profile)) call tr_deallocate(T_profile,"T_profile",CAT_GRID)
  if (allocated(density_profile)) call tr_deallocate(density_profile,"density_profile",CAT_GRID)
  
end if ! my_id == 0

if (freeboundary_equil) then
  call broadcast_elements(my_id, element_list)
  call broadcast_nodes(my_id, node_list)  !--- This is required for boundary_check
  call broadcast_boundary(my_id, bnd_elm_list, bnd_node_list)
  call boundary_check(my_id)
  deallocate(response_m_eq)
endif

equil_initialized = .true.

return

contains

!-----------------------------------------------------------------------
!> One outer q-matching update for the kinetic RE equilibrium: evaluate
!> q(psihat_n) on the present psi by flux-surface integration (the standard
!> machinery), then transplant-update Nprof.
!>
!> Shared by the fixed-boundary phase and the free-boundary phase. Kept as
!> ONE copy deliberately: the two phases differ only in what has moved the
!> boundary beforehand, and a duplicated version would silently drift --
!> in particular the ph_top rule below, whose absence caused the edge
!> current bump on the 100 keV case and would be far harder to spot on top
!> of a moving free boundary.
!>
!> Rank 0 only. The GS assembly that consumes Nprof is rank-0 too
!> (poisson sets a_mat%comm = MPI_COMM_SELF and does the whole assembly and
!> solve inside my_id == 0), so the module state never has to leave task 0;
!> only the loop VERDICT has to be broadcast, or the ranks would disagree
!> about when to leave a loop containing the collective vacuum_equil call.
subroutine re_eq_q_transplant(n_inner)
  implicit none
  integer, intent(in) :: n_inner

  call update_equil_state(my_id,node_list, element_list, bnd_elm_list, xpoint, xcase)
  call re_eq_update_labels(my_id, node_list, element_list, bnd_node_list)

  ! --- LCFS geometry, measured once per OUTER iteration.
  !     Called directly rather than by ungating the call inside
  !     update_equil_state: that one is guarded by equil_initialized, which is
  !     only set at the very end of this routine, so the LCFS fields read zero
  !     throughout the solve (they did in every RE run so far). Ungating it
  !     there would trace a flux surface on every INNER Picard iteration --
  !     several hundred per run -- and would charge that to every JOREK
  !     equilibrium, RE or not. Here it costs one trace per outer iteration and
  !     touches no shared path.
  !     R_geo and a are the quantities in which "the same plasma at a different
  !     RE energy" is meaningful: unlike the magnetic axis they do not move
  !     with the internal Shafranov shift, which grows strongly with the drift
  !     parameter (measured: the drift axis moves 4.5 cm out from 100 keV to
  !     10 MeV while R_axis was held fixed at 2.59).
  call LCFS_shape_parameters(node_list, element_list)
  write(*,'(A,F9.5,A,F9.5,A,F9.5)') ' re_eq: LCFS  R_geo = ', ES%LCFS_Rgeo, &
    '   a = ', ES%LCFS_a, '   kappa = ', ES%LCFS_kappa
  ! The limiter contact is logged with it because the size control leans on the
  ! inboard edge being held by the wall: if this point wanders, the mechanism
  ! is not what we think it is and the response will drift. Measured on the JET
  ! case it moves ~1 mm in R over a whole solve, i.e. it is a fixed geometric
  ! feature -- but that is a property of the case, not a guarantee.
  write(*,'(A,F9.5,A,F9.5,A,F9.5)') '        limiter contact R = ', ES%R_lim, &
    '   Z = ', ES%Z_lim, '   inboard LCFS edge = ', ES%LCFS_Rgeo - ES%LCFS_a

  n_lev_q = re_eq_n_q_levels
  surface_list_q%n_psi = n_lev_q + 1     ! entry 1 (magnetic axis) is skipped by determine_q_profile
  if (allocated(surface_list_q%psi_values)) call tr_deallocate(surface_list_q%psi_values,"surface_list_q%psi_values",CAT_GRID)
  call tr_allocate(surface_list_q%psi_values,1,surface_list_q%n_psi,"surface_list_q%psi_values",CAT_GRID)
  if (.not. allocated(ph_lev)) then
    call tr_allocate(ph_lev, 1,n_lev_q,               "ph_lev", CAT_GRID)
    call tr_allocate(q_lev,  1,surface_list_q%n_psi,  "q_lev",  CAT_GRID)
    call tr_allocate(rad_lev,1,surface_list_q%n_psi,  "rad_lev",CAT_GRID)
  endif
  ! q evaluation levels in psihat. For a diverted (X-point) case the top
  ! level is pulled in from 0.985 to 0.95: q -> infinity at the separatrix
  ! and the flux-surface tracer should not be asked to follow surfaces
  ! hugging it. The controllable range is well below this anyway (the beam
  ! edge, and with l_beam<1 the vacuum annulus), so nothing matchable is
  ! lost -- q between the beam edge and the separatrix is an outcome.
  !
  ! Limiter case: the top level must REACH the outermost controllable label,
  ! otherwise re_eq_outer_update clamps every label beyond it to the same
  ! argument and that whole band receives one identical, psihat-unresolved
  ! push -- which, with the absorbing edge pinning the last label to zero,
  ! is the edge current bump seen on the 100 keV hollow-q case (labels ran
  ! to psihat_n = 0.9888 against a top level of 0.985). re_eq_ph_beam_max is
  ! the previous iteration's value and is 0 before the first outer update,
  ! hence the 0.985 floor; the 0.995 cap keeps the flux-surface tracer off
  ! the boundary.
  ph_top = merge(0.95d0, min(max(0.985d0, re_eq_ph_beam_max), 0.995d0), xpoint2)
  do i_lev = 1, n_lev_q
    ph_lev(i_lev) = 0.02d0 + (ph_top - 0.02d0) &
                            * dble(i_lev-1) / dble(n_lev_q-1)
    surface_list_q%psi_values(i_lev+1) = ES%psi_axis + ph_lev(i_lev) * (ES%psi_bnd - ES%psi_axis)
  enddo
  surface_list_q%psi_values(1) = ES%psi_axis + 0.01d0 * (ES%psi_bnd - ES%psi_axis)

  call find_flux_surfaces(my_id,xpoint2,xcase2,node_list,element_list,surface_list_q)
  call determine_q_profile(node_list,element_list,surface_list_q,ES%psi_axis,ES%psi_xpoint,ES%Z_xpoint, &
                           q_lev,rad_lev)
  if (allocated(surface_list_q%flux_surfaces)) deallocate(surface_list_q%flux_surfaces)

  call re_eq_outer_update(my_id, node_list, element_list, n_lev_q, ph_lev, q_lev(2:n_lev_q+1), &
                          n_inner, re_eq_converged)

end subroutine re_eq_q_transplant

end subroutine equilibrium
