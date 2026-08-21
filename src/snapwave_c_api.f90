module snapwave_c_api
   !************************************************************************
   ! Coarse C ABI facade over the existing SnapWave Fortran model.
   !
    ! Exposes the bind(C) entry point snapwave_run_c, which mirrors the
    ! lifecycle of the stand-alone program in src/snapwave.f90:
   !    1. load the Rust-resolved configuration
   !    2. initialize domain
   !    3. read observation points
   !    4. read boundary conditions (incl. wind if specified)
   !    5. initialize NetCDF output
   !    6. run timesteps
   !    7. finalize output
   !
    ! The solver internals and the module global state are unchanged; the
    ! Rust wrapper provides main() and calls this facade. Since plan.md
    ! Phase 4 the wrapper resolves the whole configuration in Rust
    ! (defaults, validation, post-processing) and passes it down as
    ! canonical key=value text — Fortran no longer reads SnapWave.inp or
    ! decides defaults on this route. All sibling input and output files
    ! are still resolved by the Fortran readers relative to the current
    ! working directory, so the caller is expected to chdir to the input
    ! file's parent directory before calling.
   !************************************************************************
   use iso_c_binding
   implicit none
   !
   ! Width of the nameobs global (character*32), see snapwave_data.
   integer, parameter :: WIDTH_NAMEOBS = 32
   !
   ! The bind(C) mirror of the Rust #[repr(C)] struct SnapWaveStateC
   ! (src_rust/state.rs). Field order and widths must stay identical on
   ! both sides (AGENTS.md FFI rule). Absent data is zero extents plus
   ! empty (zero-length) buffers, never null pointers.
   !
   type, bind(C) :: snapwave_state_t
      ! ---- mesh
      integer(c_int) :: no_nodes
      integer(c_int) :: no_faces
      integer(c_int) :: max_nodes
      integer(c_int) :: sferic
      type(c_ptr)    :: x
      type(c_ptr)    :: y
      type(c_ptr)    :: zb
      type(c_ptr)    :: msk
      type(c_ptr)    :: face_nodes
      ! ---- boundary enclosure polyline
      integer(c_int) :: n_bndenc
      type(c_ptr)    :: x_bndenc
      type(c_ptr)    :: y_bndenc
      ! ---- neumann polyline
      integer(c_int) :: n_neu
      type(c_ptr)    :: x_neu
      type(c_ptr)    :: y_neu
      ! ---- observation points
      integer(c_int) :: nobs
      type(c_ptr)    :: xobs
      type(c_ptr)    :: yobs
      type(c_ptr)    :: names
      ! ---- boundary forcing series
      integer(c_int) :: boundary_mode
      integer(c_int) :: nwbnd
      integer(c_int) :: ntwbnd
      type(c_ptr)    :: x_bwv
      type(c_ptr)    :: y_bwv
      type(c_ptr)    :: t_bwv
      type(c_ptr)    :: hs_bwv
      type(c_ptr)    :: tp_bwv
      type(c_ptr)    :: wd_bwv
      type(c_ptr)    :: ds_bwv
      type(c_ptr)    :: zs_bwv
   end type snapwave_state_t
contains
   function snapwave_run_c(config, config_len) bind(C, name="snapwave_run_c") result(status)
      use snapwave_data
      use snapwave_input
      use snapwave_domain
      use snapwave_boundaries
      use snapwave_solver
      use snapwave_ncoutput
      use snapwave_obspoints
      use snapwave_results
      use omp_lib
      implicit none

      type(c_ptr), value :: config
      integer(c_int), value :: config_len
      integer(c_int) :: status

      character(len=:), allocatable :: ftext
      character(kind=c_char), dimension(:), pointer :: cchars
      integer :: i
      integer :: du
      integer :: ios

      ! plan.md Phase 4: the Rust wrapper resolves the entire configuration
      ! (defaults, validation, post-processing) and passes it here as
      ! canonical key=value text. Fortran no longer reads SnapWave.inp or
      ! decides defaults on this route — it only stores what Rust resolved.
      allocate (character(len=config_len) :: ftext)
      call c_f_pointer(config, cchars, [config_len])
      do i = 1, config_len
         ftext(i:i) = achar(iachar(cchars(i)))
      end do

      open (newunit=du, status='scratch', action='readwrite', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_c_api: cannot open scratch file for config'
         status = 1_c_int
         return
      end if
      call write_config_lines(du, ftext)
      call read_resolved_input(du)
      close (du)
      !
      call initialize_snapwave_domain()     ! Read mesh, finds upwind neighbors, etc.
      !
      call read_obs_points()
      !
      ! Read boundary conditions
      !
      call read_boundary_data()
      !
      ! read wind data if specified
      !
      call read_wind_data()
      !
      ! Initialize NetCDF output
      !
      call ncoutput_init()
      !
      ! Start time loop
      !
      call run_time_loop()
      !
      call ncoutput_finalize()
      !
      status = 0_c_int
      !
    end function snapwave_run_c

   !************************************************************************
   ! plan.md Phase 7: run the model without NetCDF output, streaming the
   ! output-time state to a capture file that the Rust wrapper replays
   ! through its own NetCDF writers (src_rust/output.rs). Shares the
   ! timestep loop with snapwave_run_c via run_time_loop; only the output
   ! sink differs (capture stream vs NetCDF).
   !************************************************************************
   function snapwave_run_capture_c(config, config_len, capture_path, capture_path_len) &
         bind(C, name="snapwave_run_capture_c") result(status)
      use snapwave_data
      use snapwave_input
      use snapwave_domain
      use snapwave_boundaries
      use snapwave_ncoutput
      use snapwave_obspoints
      implicit none

      type(c_ptr), value :: config
      integer(c_int), value :: config_len
      type(c_ptr), value :: capture_path
      integer(c_int), value :: capture_path_len
      integer(c_int) :: status

      character(len=:), allocatable :: ftext
      character(len=1024) :: cpath
      character(kind=c_char), dimension(:), pointer :: cchars
      integer :: i, du, ios, cu

      status = 1_c_int

      allocate (character(len=config_len) :: ftext)
      call c_f_pointer(config, cchars, [config_len])
      do i = 1, config_len
         ftext(i:i) = achar(iachar(cchars(i)))
      end do
      call c_f_pointer(capture_path, cchars, [capture_path_len])
      cpath = ' '
      do i = 1, capture_path_len
         cpath(i:i) = achar(iachar(cchars(i)))
      end do

      open (newunit=du, status='scratch', action='readwrite', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_run_capture_c: cannot open scratch file for config'
         return
      end if
      call write_config_lines(du, ftext)
      call read_resolved_input(du)
      close (du)

      call initialize_snapwave_domain()
      call read_obs_points()
      call read_boundary_data()
      call read_wind_data()

      open (newunit=cu, file=trim(cpath), form='unformatted', access='stream', &
            action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_run_capture_c: cannot open capture file: ', trim(cpath)
         return
      end if

      call run_capture_tail(cu)

      close (cu)
      status = 0_c_int
      !
   end function snapwave_run_capture_c

   !************************************************************************
   ! plan.md Phase 8: coarse Fortran entry point that consumes
   ! Rust-prepared state. Mirrors snapwave_run_capture_c, except the
   ! mesh, the enclosure/neumann polylines, the observation points and
   ! the boundary series are not read from disk: they were parsed and
   ! validated in Rust (plan.md Phases 6-7) and are associated here
   ! with Rust-owned memory through associate_state_from_rust, so the
   ! allocation ownership of those non-solver arrays is Rust's
   ! (Phase 8, steps 4-5). Derived quantities (surrounding points,
   ! upwind neighbours, make_map_fm / find_boundary_indices weights)
   ! and the remaining file-backed inputs (wind, fw/fwig, vegetation)
   ! stay Fortran until the Phase 9 interpolation migration; the legacy
   ! make-built oracle keeps reading every file itself.
   !************************************************************************
   function snapwave_run_capture_state_c(config, config_len, capture_path, capture_path_len, state) &
         bind(C, name="snapwave_run_capture_state_c") result(status)
      use snapwave_data
      use snapwave_input
      use snapwave_domain
      use snapwave_boundaries
      use snapwave_ncoutput
      use snapwave_obspoints
      implicit none

      type(c_ptr), value :: config
      integer(c_int), value :: config_len
      type(c_ptr), value :: capture_path
      integer(c_int), value :: capture_path_len
      type(c_ptr), value :: state
      integer(c_int) :: status

      character(len=:), allocatable :: ftext
      character(len=1024) :: cpath
      character(kind=c_char), dimension(:), pointer :: cchars
      type(snapwave_state_t), pointer :: st
      integer :: i, du, ios, cu, istat

      status = 1_c_int

      allocate (character(len=config_len) :: ftext)
      call c_f_pointer(config, cchars, [config_len])
      do i = 1, config_len
         ftext(i:i) = achar(iachar(cchars(i)))
      end do
      call c_f_pointer(capture_path, cchars, [capture_path_len])
      cpath = ' '
      do i = 1, capture_path_len
         cpath(i:i) = achar(iachar(cchars(i)))
      end do

      open (newunit=du, status='scratch', action='readwrite', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_run_capture_state_c: cannot open scratch file for config'
         return
      end if
      call write_config_lines(du, ftext)
      call read_resolved_input(du)
      close (du)

      ! Associate the snapwave_data globals with the Rust-owned buffers.
      call c_f_pointer(state, st)
      call associate_state_from_rust(st, istat)
      if (istat /= 0) then
         write (*, *) 'ERROR: snapwave_run_capture_state_c: rejecting Rust state'
         return
      end if

      call initialize_snapwave_domain(mesh_from_rust=.true.)
      call init_obs_points_from_state()
      call init_boundary_from_state()
      call read_wind_data()

      open (newunit=cu, file=trim(cpath), form='unformatted', access='stream', &
            action='write', status='replace', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_run_capture_state_c: cannot open capture file: ', trim(cpath)
         return
      end if

      call run_capture_tail(cu)

      close (cu)
      status = 0_c_int
      !
   end function snapwave_run_capture_state_c

!************************************************************************
    ! plan.md Phase 10: Rust-owned time loop and output scheduling.
    !
    ! Instead of the Fortran run_time_loop driving the model, Rust now
    ! owns the time loop and calls these coarse entry points:
    !   - snapwave_init_capture_c       (Rust-state route: init + capture)
    !   - snapwave_init_legacy_capture_c (legacy route: init + capture)
    !   - snapwave_timestep_c           (one solver step)
    !   - snapwave_capture_map_c        (capture map output at current t)
    !   - snapwave_capture_his_c        (capture his output at current t)
    !   - snapwave_finalize_capture_c   (close capture stream)
    !
    ! The existing snapwave_run_capture_c / snapwave_run_capture_state_c
    ! stay available for the --legacy-mesh parity hook and the comparison
    ! hooks; they are unchanged.
    !************************************************************************

    !************************************************************************
    ! Phase 10: initialize the model with Rust-owned state, open the
    ! capture stream and write the static output data, but do NOT run the
    ! time loop. The caller (Rust) will drive the loop through
    ! snapwave_timestep_c / snapwave_capture_{map,his}_c and finalize
    ! with snapwave_finalize_capture_c.
    !
    ! Mirrors snapwave_run_capture_state_c up to (and including)
    ! ncoutput_init, then returns — the time loop is Rust-owned.
    !************************************************************************
    function snapwave_init_capture_c(config, config_len, capture_path, capture_path_len, state) &
          bind(C, name="snapwave_init_capture_c") result(status)
       use snapwave_data
       use snapwave_input
       use snapwave_domain
       use snapwave_boundaries
       use snapwave_ncoutput
       use snapwave_obspoints
       implicit none

       type(c_ptr), value :: config
       integer(c_int), value :: config_len
       type(c_ptr), value :: capture_path
       integer(c_int), value :: capture_path_len
       type(c_ptr), value :: state
       integer(c_int) :: status

       character(len=:), allocatable :: ftext
       character(len=1024) :: cpath
       character(kind=c_char), dimension(:), pointer :: cchars
       type(snapwave_state_t), pointer :: st
       integer :: i, du, ios, cu, istat

       status = 1_c_int

       allocate (character(len=config_len) :: ftext)
       call c_f_pointer(config, cchars, [config_len])
       do i = 1, config_len
          ftext(i:i) = achar(iachar(cchars(i)))
       end do
       call c_f_pointer(capture_path, cchars, [capture_path_len])
       cpath = ' '
       do i = 1, capture_path_len
          cpath(i:i) = achar(iachar(cchars(i)))
       end do

       open (newunit=du, status='scratch', action='readwrite', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_init_capture_c: cannot open scratch file for config'
          return
       end if
       call write_config_lines(du, ftext)
       call read_resolved_input(du)
       close (du)

       ! Associate the snapwave_data globals with the Rust-owned buffers.
       call c_f_pointer(state, st)
       call associate_state_from_rust(st, istat)
       if (istat /= 0) then
          write (*, *) 'ERROR: snapwave_init_capture_c: rejecting Rust state'
          return
       end if

       call initialize_snapwave_domain(mesh_from_rust=.true.)
       call init_obs_points_from_state()
       call init_boundary_from_state()
       call read_wind_data()

       open (newunit=cu, file=trim(cpath), form='unformatted', access='stream', &
             action='write', status='replace', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_init_capture_c: cannot open capture file: ', trim(cpath)
          return
       end if

       call run_init_tail(cu)

       status = 0_c_int
       !
    end function snapwave_init_capture_c

    !************************************************************************
    ! Phase 10: initialize the model by reading files (legacy route),
    ! open the capture stream and write the static output data, but do
    ! NOT run the time loop. Mirrors snapwave_run_capture_c up to (and
    ! including) ncoutput_init.
    !************************************************************************
    function snapwave_init_legacy_capture_c(config, config_len, capture_path, capture_path_len) &
          bind(C, name="snapwave_init_legacy_capture_c") result(status)
       use snapwave_data
       use snapwave_input
       use snapwave_domain
       use snapwave_boundaries
       use snapwave_ncoutput
       use snapwave_obspoints
       implicit none

       type(c_ptr), value :: config
       integer(c_int), value :: config_len
       type(c_ptr), value :: capture_path
       integer(c_int), value :: capture_path_len
       integer(c_int) :: status

       character(len=:), allocatable :: ftext
       character(len=1024) :: cpath
       character(kind=c_char), dimension(:), pointer :: cchars
       integer :: i, du, ios, cu

       status = 1_c_int

       allocate (character(len=config_len) :: ftext)
       call c_f_pointer(config, cchars, [config_len])
       do i = 1, config_len
          ftext(i:i) = achar(iachar(cchars(i)))
       end do
       call c_f_pointer(capture_path, cchars, [capture_path_len])
       cpath = ' '
       do i = 1, capture_path_len
          cpath(i:i) = achar(iachar(cchars(i)))
       end do

       open (newunit=du, status='scratch', action='readwrite', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_init_legacy_capture_c: cannot open scratch file for config'
          return
       end if
       call write_config_lines(du, ftext)
       call read_resolved_input(du)
       close (du)

       call initialize_snapwave_domain()
       call read_obs_points()
       call read_boundary_data()
       call read_wind_data()

       open (newunit=cu, file=trim(cpath), form='unformatted', access='stream', &
             action='write', status='replace', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_init_legacy_capture_c: cannot open capture file: ', trim(cpath)
          return
       end if

       call run_init_tail(cu)

       status = 0_c_int
       !
    end function snapwave_init_legacy_capture_c

    !************************************************************************
    ! Phase 10: one solver timestep. Calls update_boundary_conditions(t)
    ! and compute_wave_field(t) — exactly what run_time_loop does per
    ! iteration. The caller (Rust) owns the time loop and output
    ! scheduling.
    !************************************************************************
    function snapwave_timestep_c(t, it) bind(C, name="snapwave_timestep_c") result(status)
       use snapwave_data
       use snapwave_boundaries
       use snapwave_solver
       implicit none

       real(c_double), value :: t
       integer(c_int), value :: it
       integer(c_int) :: status

       call update_boundary_conditions(t)
       call compute_wave_field(t)

       status = 0_c_int
       !
    end function snapwave_timestep_c

    !************************************************************************
    ! Phase 10: capture map output at the current time. Calls
    ! ncoutput_update_map(t, ntmapout), which writes to the capture
    ! stream when capture mode is active (set by run_init_tail).
    !************************************************************************
    function snapwave_capture_map_c(t, ntmapout) bind(C, name="snapwave_capture_map_c") result(status)
       use snapwave_data
       use snapwave_ncoutput
       implicit none

       real(c_double), value :: t
       integer(c_int), value :: ntmapout
       integer(c_int) :: status

       call ncoutput_update_map(t, ntmapout)

       status = 0_c_int
       !
    end function snapwave_capture_map_c

    !************************************************************************
    ! Phase 10: capture history output at the current time. Calls
    ! update_obs_points() then ncoutput_update_his(t, nthisout), which
    ! writes to the capture stream when capture mode is active.
    !************************************************************************
    function snapwave_capture_his_c(t, nthisout) bind(C, name="snapwave_capture_his_c") result(status)
       use snapwave_data
       use snapwave_ncoutput
       use snapwave_obspoints
       implicit none

       real(c_double), value :: t
       integer(c_int), value :: nthisout
       integer(c_int) :: status

       call update_obs_points()
       call ncoutput_update_his(t, nthisout)

       status = 0_c_int
       !
    end function snapwave_capture_his_c

    !************************************************************************
    ! Phase 10: finalize the capture run. Calls ncoutput_finalize(),
    ! closes the capture stream file (which was opened by the init
    ! function — the unit is saved in snapwave_ncoutput's capture_unit),
    ! and resets capture mode.
    !************************************************************************
    function snapwave_finalize_capture_c() bind(C, name="snapwave_finalize_capture_c") result(status)
       use snapwave_ncoutput, only: ncoutput_finalize, ncoutput_capture_end, capture_unit
       implicit none

       integer(c_int) :: status
       integer :: ios

       call ncoutput_finalize()
       close (capture_unit, iostat=ios)
       call ncoutput_capture_end()

       status = 0_c_int
       !
    end function snapwave_finalize_capture_c

    !************************************************************************
    ! Shared Phase 10 init tail: open the capture stream and write the
    ! static output data (ncoutput_capture_begin + ncoutput_init). Called
    ! by both snapwave_init_capture_c and snapwave_init_legacy_capture_c
    ! after the domain/boundary/wind initialisation is done.
    !************************************************************************
    subroutine run_init_tail(cu)
       !
       use snapwave_ncoutput
       implicit none
       !
       integer, intent(in) :: cu
       !
       call ncoutput_capture_begin(cu)
       call ncoutput_init()
       !
    end subroutine run_init_tail

    !************************************************************************
    ! Associate the Phase 8 snapwave_data globals with the Rust-owned
    ! buffers of one snapwave_state_t (plan.md Phase 8, step 4: the
    ! allocation ownership of these non-solver arrays is Rust's; the
    ! Fortran side only points at the memory). nameobs is the one
    ! exception: character(len=32) arrays do not associate portably, so
    ! the names are copied from the packed 32-byte records.
    !
    ! Buffer layouts (pinned by tests in src_rust/ffi_layout.rs):
    !   face_nodes : (4, no_faces)      column-major
    !   *_bwv      : (nwbnd, ntwbnd)    column-major
    !   msk        : integer*1 (c_int8_t)
    !************************************************************************
    subroutine associate_state_from_rust(st, status)
      use snapwave_data
      implicit none

      type(snapwave_state_t), intent(in) :: st
      integer, intent(out) :: status

      character(kind=c_char), dimension(:), pointer :: cname
      integer :: i, j
      integer(c_size_t) :: sz_nodes, sz_faces, sz_enc, sz_neu, sz_obs, sz_bnd, sz_t

      status = 1

      if (st%no_nodes < 0 .or. st%no_faces < 0 .or. st%max_nodes < 0 .or. st%nobs < 0 &
            .or. st%n_bndenc < 0 .or. st%n_neu < 0 .or. st%nwbnd < 0 .or. st%ntwbnd < 0) then
         write (*, *) 'ERROR: associate_state_from_rust: negative array extent in Rust state'
         return
      end if

      ! ---- mesh
      no_nodes = st%no_nodes
      no_faces = st%no_faces
      max_nodes = st%max_nodes
      sferic    = st%sferic
      sz_nodes = int(no_nodes, c_size_t)
      sz_faces = int(no_faces, c_size_t)
      call c_f_pointer(st%x, x, [sz_nodes])
      call c_f_pointer(st%y, y, [sz_nodes])
      call c_f_pointer(st%zb, zb, [sz_nodes])
      call c_f_pointer(st%msk, msk, [sz_nodes])
      call c_f_pointer(st%face_nodes, face_nodes, [4_c_size_t, sz_faces])
      ! xs/ys are legacy-only (allocated, never filled); the state route
      ! leaves them unallocated, matching their unused status.

      ! ---- polylines
      n_bndenc = st%n_bndenc
      sz_enc = int(n_bndenc, c_size_t)
      call c_f_pointer(st%x_bndenc, x_bndenc, [sz_enc])
      call c_f_pointer(st%y_bndenc, y_bndenc, [sz_enc])
      n_neu = st%n_neu
      sz_neu = int(n_neu, c_size_t)
      call c_f_pointer(st%x_neu, x_neu, [sz_neu])
      call c_f_pointer(st%y_neu, y_neu, [sz_neu])

      ! ---- observation points
      nobs = st%nobs
      sz_obs = int(nobs, c_size_t)
      call c_f_pointer(st%xobs, xobs, [sz_obs])
      call c_f_pointer(st%yobs, yobs, [sz_obs])
      if (nobs > 0) then
         call c_f_pointer(st%names, cname, [int(WIDTH_NAMEOBS, c_size_t) * sz_obs])
         allocate (nameobs(nobs))
         do i = 1, nobs
            do j = 1, WIDTH_NAMEOBS
               nameobs(i) (j:j) = achar(iachar(cname(WIDTH_NAMEOBS*(i - 1) + j)))
            end do
         end do
      end if

      ! ---- boundary forcing series (wd/ds already converted to radians)
      nwbnd  = st%nwbnd
      ntwbnd = st%ntwbnd
      sz_bnd = int(nwbnd, c_size_t)
      sz_t = int(ntwbnd, c_size_t)
      call c_f_pointer(st%x_bwv, x_bwv, [sz_bnd])
      call c_f_pointer(st%y_bwv, y_bwv, [sz_bnd])
      call c_f_pointer(st%t_bwv, t_bwv, [sz_t])
      call c_f_pointer(st%hs_bwv, hs_bwv, [sz_bnd, sz_t])
      call c_f_pointer(st%tp_bwv, tp_bwv, [sz_bnd, sz_t])
      call c_f_pointer(st%wd_bwv, wd_bwv, [sz_bnd, sz_t])
      call c_f_pointer(st%ds_bwv, ds_bwv, [sz_bnd, sz_t])
      call c_f_pointer(st%zs_bwv, zs_bwv, [sz_bnd, sz_t])

      status = 0
      !
   end subroutine associate_state_from_rust

   subroutine run_capture_tail(cu)
      !
      ! Shared Phase 7 capture tail of snapwave_run_capture_c and
      ! snapwave_run_capture_state_c: initialize the output in capture
      ! mode, run the timestep loop, finalize.
      !
      use snapwave_ncoutput
      implicit none
      !
      integer, intent(in) :: cu
      !
      call ncoutput_capture_begin(cu)
      call ncoutput_init()
      call run_time_loop()
      call ncoutput_finalize()
      call ncoutput_capture_end()
      !
   end subroutine run_capture_tail

   subroutine run_time_loop()
      !
      ! Timestep loop shared by snapwave_run_c (NetCDF output),
      ! snapwave_run_capture_c and snapwave_run_capture_state_c (Phase 7/8
      ! capture output). The model stepping and the output scheduling are
      ! identical; the facades differ only in how ncoutput_* handles each
      ! output time and in where the input data comes from (files vs the
      ! Rust-owned state).
      !
      use snapwave_data
      use snapwave_boundaries
      use snapwave_solver
      use snapwave_ncoutput
      use snapwave_obspoints
      use omp_lib
      implicit none
      !
      real*8  :: t
      real*8  :: output_tol
      integer :: it
      !
      !$omp parallel
      !$omp single
      write(*,'(A,I2,A)') 'Running with ', omp_get_num_threads(), ' OMP threads.'
      !$omp end single
      !$omp end parallel
      it = 0
      t  = tstart
      map_output_count = 0
      his_output_count = 0
      next_map_output = tstart
      next_his_output = tstart
      output_tol = max(1.0d-6, abs(dble(timestep))*1.0d-6)
      !
      write(*,*)'Start time loop'
      do while (t<=tstop)
         !
         ! New time step
         !
         it = it + 1
         !
         call update_boundary_conditions(t) ! includes theta_grid creation
         !
         call compute_wave_field(t)
         !
         if (his_filename /= '' .and. nobs > 0 .and. t >= next_his_output - output_tol) then
            his_output_count = his_output_count + 1
            call update_obs_points()
            call ncoutput_update_his(t, his_output_count)
            do while (next_his_output <= t + output_tol)
               next_his_output = next_his_output + dble(his_interval)
            end do
         end if
         !
         if (ja_save_each_iter == 0 .and. map_filename /= '' .and. t >= next_map_output - output_tol) then
            map_output_count = map_output_count + 1
            call ncoutput_update_map(t, map_output_count)
            do while (next_map_output <= t + output_tol)
               next_map_output = next_map_output + dble(map_interval)
            end do
         end if
         !
         t = t + timestep
         !
      enddo
      !
   end subroutine run_time_loop

   !************************************************************************
   ! plan.md Phase 7 step 2 hook: read the mesh NetCDF with the unchanged
   ! nc_read_net reader (plus the two initialize_snapwave_domain
   ! post-processing steps) and dump the resulting globals, so the Rust port
   ! (src_rust/mesh.rs) can be pinned against the numerical oracle.
   !************************************************************************
   function snapwave_mesh_dump_c(config, config_len, dump_path, dump_path_len) &
         bind(C, name="snapwave_mesh_dump_c") result(status)
      use snapwave_data
      use snapwave_input
      use snapwave_ncoutput, only: nc_read_net
      implicit none

      type(c_ptr), value :: config
      integer(c_int), value :: config_len
      type(c_ptr), value :: dump_path
      integer(c_int), value :: dump_path_len
      integer(c_int) :: status

      character(len=:), allocatable :: ftext
      character(len=1024) :: dpath
      character(kind=c_char), dimension(:), pointer :: cchars
      integer :: i, ios, dunit, du, k

      status = 1_c_int

      allocate (character(len=config_len) :: ftext)
      call c_f_pointer(config, cchars, [config_len])
      do i = 1, config_len
         ftext(i:i) = achar(iachar(cchars(i)))
      end do
      call c_f_pointer(dump_path, cchars, [dump_path_len])
      dpath = ' '
      do i = 1, dump_path_len
         dpath(i:i) = achar(iachar(cchars(i)))
      end do

      open (newunit=du, status='scratch', action='readwrite', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_mesh_dump_c: cannot open scratch file for config'
         return
      end if
      call write_config_lines(du, ftext)
      call read_resolved_input(du)
      close (du)

      ! Read the mesh, then reproduce the two initialize_snapwave_domain
      ! post-processing steps (zb sign flip, fourth-node sentinel).
      call nc_read_net()
      zb = -posdwn*zb
      do k = 1, no_faces
         if (face_nodes(4, k) == 0) face_nodes(4, k) = -999
      end do

      open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_mesh_dump_c: cannot open dump file: ', trim(dpath)
         return
      end if

      call dump_mesh_globals(dunit)
      close (dunit)
      !
      status = 0_c_int
      !
   end function snapwave_mesh_dump_c

   subroutine dump_mesh_globals(u)
      !
      ! Dump the mesh globals produced by nc_read_net + post-processing, in
      ! the flat format parsed by src_rust/mesh.rs. face_nodes is dumped for
      ! rows 1..max_nodes only (the deterministic part; for a pure-triangle
      ! mesh the fourth row is never written).
      !
      use snapwave_data
      implicit none
      integer, intent(in) :: u
      integer :: k, j
      !
      write (u, '(A,1X,I0)') 'no_nodes', no_nodes
      write (u, '(A,1X,I0)') 'no_faces', no_faces
      write (u, '(A,1X,I0)') 'max_nodes', max_nodes
      write (u, '(A,1X,I0)') 'sferic', sferic
      !
      write (u, '(A,1X,I0)') 'x', no_nodes
      do k = 1, no_nodes
         call dump_r8_line(u, x(k))
      end do
      write (u, '(A,1X,I0)') 'y', no_nodes
      do k = 1, no_nodes
         call dump_r8_line(u, y(k))
      end do
      write (u, '(A,1X,I0)') 'zb', no_nodes
      do k = 1, no_nodes
         call dump_r4_line(u, zb(k))
      end do
      write (u, '(A,1X,I0)') 'msk', no_nodes
      do k = 1, no_nodes
         write (u, '(I0)') msk(k)
      end do
      write (u, '(A,1X,I0)') 'face_nodes', max_nodes*no_faces
      do k = 1, no_faces
         do j = 1, max_nodes
            write (u, '(I0)') face_nodes(j, k)
         end do
      end do
      !
   end subroutine dump_mesh_globals

   !************************************************************************
   ! Comparison hooks (plan.md Phase 3 step 5 and Phase 4).
   !
   ! Both dumps render the resolved snapwave_data globals as canonical
   ! key=value text so the Rust wrapper can compare its own parse result
   ! against the Fortran side:
   !   - snapwave_read_input_dump_c  parses the input file with the legacy
   !     Fortran reader (the numerical oracle) and dumps the globals;
   !   - snapwave_load_config_dump_c  loads the Rust-resolved configuration
   !     text (read_resolved_input) and dumps the same globals, pinning the
   !     Phase 4 config handoff.
   ! The model is not run by either hook.
   !
   ! Dump conventions (mirrored by src_rust/input_compare.rs):
   !   integer   -> decimal (I0)
   !   real*4/8  -> IEEE-754 bit pattern in zero-padded hex (Z8.8 / Z16.16),
   !                so values compare without decimal formatting drift
   !   logical   -> 1 / 0
   !   character -> trimmed value after the first '=' (keys never contain
   !                '=', values may contain blanks)
   !************************************************************************
   function snapwave_read_input_dump_c(input_path, input_path_len, dump_path, dump_path_len) &
         bind(C, name="snapwave_read_input_dump_c") result(status)
      use snapwave_data
      use snapwave_input
      implicit none

      type(c_ptr), value :: input_path
      integer(c_int), value :: input_path_len
      type(c_ptr), value :: dump_path
      integer(c_int), value :: dump_path_len
      integer(c_int) :: status

      character(len=1024) :: fpath, dpath
      character(kind=c_char), dimension(:), pointer :: cchars
      integer :: i, ios, dunit, istat
      logical :: exists

      status = 1_c_int

      call c_f_pointer(input_path, cchars, [input_path_len])
      fpath = ' '
      do i = 1, input_path_len
         fpath(i:i) = achar(iachar(cchars(i)))
      end do
      call c_f_pointer(dump_path, cchars, [dump_path_len])
      dpath = ' '
      do i = 1, dump_path_len
         dpath(i:i) = achar(iachar(cchars(i)))
      end do

      ! The Rust wrapper validates existence before calling; fail cleanly
      ! (status return, not a runtime error) if the contract is broken.
      inquire (file=trim(fpath), exist=exists)
      if (.not. exists) then
         write (*, *) 'ERROR: snapwave_read_input_dump_c: input file not found: ', trim(fpath)
         return
      end if

      call read_snapwave_input(input_file=trim(fpath), status=istat)
      if (istat /= 0) then
         write (*, *) 'ERROR: snapwave_read_input_dump_c: input reader failed'
         return
      end if

      open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_read_input_dump_c: cannot open dump file: ', trim(dpath)
         return
      end if

      call dump_globals(dunit)

      close (dunit)
      !
      status = 0_c_int
      !
   end function snapwave_read_input_dump_c

   !************************************************************************
   ! plan.md Phase 4 hook: load the Rust-resolved configuration text and
   ! dump the resulting globals, so `--compare-input` can pin that the
   ! config handoff (Rust -> text -> Fortran globals) round-trips exactly.
   !************************************************************************
   function snapwave_load_config_dump_c(config, config_len, dump_path, dump_path_len) &
         bind(C, name="snapwave_load_config_dump_c") result(status)
      use snapwave_data
      use snapwave_input
      implicit none

      type(c_ptr), value :: config
      integer(c_int), value :: config_len
      type(c_ptr), value :: dump_path
      integer(c_int), value :: dump_path_len
      integer(c_int) :: status

      character(len=:), allocatable :: ftext
      character(len=1024) :: dpath
      character(kind=c_char), dimension(:), pointer :: cchars
      integer :: i, ios, dunit, du

      status = 1_c_int

      allocate (character(len=config_len) :: ftext)
      call c_f_pointer(config, cchars, [config_len])
      do i = 1, config_len
         ftext(i:i) = achar(iachar(cchars(i)))
      end do
      call c_f_pointer(dump_path, cchars, [dump_path_len])
      dpath = ' '
      do i = 1, dump_path_len
         dpath(i:i) = achar(iachar(cchars(i)))
      end do

      open (newunit=du, status='scratch', action='readwrite', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_load_config_dump_c: cannot open scratch file for config'
         return
      end if
      call write_config_lines(du, ftext)
      call read_resolved_input(du)
      close (du)

      open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: snapwave_load_config_dump_c: cannot open dump file: ', trim(dpath)
         return
      end if

      call dump_globals(dunit)

      close (dunit)
      !
      status = 0_c_int
      !
   end function snapwave_load_config_dump_c

   subroutine dump_globals(u)
      !
      ! Dump every snapwave_data global that the input readers set, in the
      ! canonical key=value format (see the module docs). Shared by both
      ! comparison hooks so the two dumps stay in lock-step.
      !
      use snapwave_data
      implicit none
      integer, intent(in) :: u
      !
      ! ---- time / control
      call dump_char (u, 'trefstr', trefstr)
      call dump_char (u, 'tstartstr', tstartstr)
      call dump_char (u, 'tstopstr', tstopstr)
      call dump_real8(u, 'tstart', tstart)
      call dump_real8(u, 'tstop', tstop)
      call dump_real4(u, 'timestep', timestep)
      call dump_real4(u, 'dt', dt)
      call dump_int  (u, 'niter', niter)
      call dump_real4(u, 'crit', crit)
      call dump_log  (u, 'restart', restart)

      ! ---- grid / domain (mmax/nmax already include the +2 dummy rows)
      call dump_int  (u, 'mmax', mmax)
      call dump_int  (u, 'nmax', nmax)
      call dump_real4(u, 'dx', dx)
      call dump_real4(u, 'dy', dy)
      call dump_real4(u, 'x0', x0)
      call dump_real4(u, 'y0', y0)
      call dump_real4(u, 'rotation', rotation)
      call dump_real4(u, 'posdwn', posdwn)
      call dump_int  (u, 'sferic', sferic)
      call dump_real4(u, 'dtheta', dtheta)
      call dump_real4(u, 'sector', sector)
      call dump_char (u, 'gridfile', gridfile)
      call dump_char (u, 'depfile', depfile)
      call dump_char (u, 'mskfile', mskfile)
      call dump_char (u, 'indfile', indfile)
      call dump_char (u, 'upwfile', upwfile)

      ! ---- boundary forcing
      call dump_char (u, 'jonswapfile', jonswapfile)
      call dump_char (u, 'bndfile', bndfile)
      call dump_char (u, 'encfile', encfile)
      call dump_char (u, 'neumannfile', neumannfile)
      call dump_char (u, 'bhsfile', bhsfile)
      call dump_char (u, 'btpfile', btpfile)
      call dump_char (u, 'bwdfile', bwdfile)
      call dump_char (u, 'bdsfile', bdsfile)
      call dump_char (u, 'bzsfile', bzsfile)
      call dump_char (u, 'obsfile', obsfile)
      call dump_real4(u, 'tol', tol)

      ! ---- wind
      call dump_char (u, 'u10str', u10str)
      call dump_char (u, 'u10dirstr', u10dirstr)
      call dump_char (u, 'windlistfile', windlistfile)
      call dump_int  (u, 'mwind', mwind)
      call dump_int  (u, 'wind', wind)

      ! ---- output
      call dump_char (u, 'map_filename', map_filename)
      call dump_char (u, 'his_filename', his_filename)
      call dump_real4(u, 'map_interval', map_interval)
      call dump_real4(u, 'his_interval', his_interval)
      call dump_int  (u, 'map_dep', map_dep)
      call dump_int  (u, 'map_Hm0', map_Hm0)
      call dump_int  (u, 'map_Hig', map_Hig)
      call dump_int  (u, 'map_Tp', map_Tp)
      call dump_int  (u, 'map_dir', map_dir)
      call dump_int  (u, 'map_dirspr', map_dirspr)
      call dump_int  (u, 'map_cg', map_cg)
      call dump_int  (u, 'map_Dw', map_Dw)
      call dump_int  (u, 'map_Df', map_Df)
      call dump_int  (u, 'map_SwE', map_SwE)
      call dump_int  (u, 'map_SwA', map_SwA)
      call dump_int  (u, 'map_sig', map_sig)
      call dump_int  (u, 'map_u10', map_u10)
      call dump_int  (u, 'map_Dveg', map_Dveg)
      call dump_int  (u, 'map_ee', map_ee)
      call dump_int  (u, 'map_ctheta', map_ctheta)
      call dump_int  (u, 'ja_save_each_iter', ja_save_each_iter)

      ! ---- diagnostics
      call dump_log  (u, 'writetestfiles', writetestfiles)

      ! ---- vegetation
      call dump_int  (u, 'ja_vegetation', ja_vegetation)
      call dump_char (u, 'vegmapfile', vegmapfile)

      ! ---- solver physics knobs
      call dump_real4(u, 'gamma', gamma)
      call dump_real4(u, 'alpha', alpha)
      call dump_real4(u, 'gammax', gammax)
      call dump_real4(u, 'hmin', hmin)
      call dump_real4(u, 'fwcutoff', fwcutoff)
      call dump_char (u, 'fwstr', fwstr)
      call dump_char (u, 'fw_igstr', fw_igstr)
      call dump_real4(u, 'Tpini', Tpini)
      call dump_real4(u, 'zsini', zsini)
      call dump_real4(u, 'sigmin', sigmin)
      call dump_real4(u, 'sigmax', sigmax)
      call dump_int  (u, 'jadcgdx', jadcgdx)
      call dump_real4(u, 'c_dispT', c_dispT)
      call dump_int  (u, 'ig', ig)
      call dump_int  (u, 'upwindref', upwindref)
      !
   end subroutine dump_globals

   subroutine dump_int(u, key, val)
      !
      integer, intent(in) :: u, val
      character(len=*), intent(in) :: key
      write (u, '(A,"=",I0)') trim(key), val
   end subroutine dump_int

   subroutine dump_real4(u, key, val)
      !
      integer, intent(in) :: u
      real*4, intent(in) :: val
      character(len=*), intent(in) :: key
      integer(c_int32_t) :: bits
      bits = transfer(val, bits)
      write (u, '(A,"=",Z8.8)') trim(key), bits
   end subroutine dump_real4

   subroutine dump_real8(u, key, val)
      !
      integer, intent(in) :: u
      real*8, intent(in) :: val
      character(len=*), intent(in) :: key
      integer(c_int64_t) :: bits
      bits = transfer(val, bits)
      write (u, '(A,"=",Z16.16)') trim(key), bits
   end subroutine dump_real8

   subroutine dump_log(u, key, val)
      !
      integer, intent(in) :: u
      logical, intent(in) :: val
      character(len=*), intent(in) :: key
      write (u, '(A,"=",I1)') trim(key), merge(1, 0, val)
   end subroutine dump_log

    subroutine dump_char(u, key, val)
       !
       integer, intent(in) :: u
       character(len=*), intent(in) :: key
       character(len=*), intent(in) :: val
       write (u, '(A,"=",A)') trim(key), trim(val)
    end subroutine dump_char

    !************************************************************************
    ! plan.md Phase 6 hook: read the auxiliary *text* inputs (observation
    ! points, boundary conditions, wind, enclosure/neumann polylines) with
    ! the unchanged Fortran readers and dump the resulting globals, so the
    ! Rust parsers (src_rust/text_input.rs) can be pinned against the
    ! numerical oracle. Mirrors snapwave_run_c up to read_wind_data but
    ! stops before the NetCDF output and the time loop.
    !
    ! Dump format (parsed by src_rust/text_compare.rs):
    !   "section <name>" then "key value" lines. Array keys (x, y, name, t,
    !   hs, tp, wd, ds, zs, u10, u10dir) are followed by a count and then
    !   one value per line: reals as IEEE-754 bit patterns in zero-padded
    !   hex (real*8 Z16.16, real*4 Z8.8), names as trimmed text. Scalar keys
    !   (n, mode, nwbnd, ntwbnd, ntu10bnd) carry a single value.
    !************************************************************************
    function snapwave_text_dump_c(config, config_len, dump_path, dump_path_len) &
          bind(C, name="snapwave_text_dump_c") result(status)
       use snapwave_data
       use snapwave_input
       use snapwave_domain
       use snapwave_boundaries
       use snapwave_obspoints
       implicit none

       type(c_ptr), value :: config
       integer(c_int), value :: config_len
       type(c_ptr), value :: dump_path
       integer(c_int), value :: dump_path_len
       integer(c_int) :: status

       character(len=:), allocatable :: ftext
       character(len=1024) :: dpath
       character(kind=c_char), dimension(:), pointer :: cchars
       integer :: i, ios, dunit, du

       status = 1_c_int

       allocate (character(len=config_len) :: ftext)
       call c_f_pointer(config, cchars, [config_len])
       do i = 1, config_len
          ftext(i:i) = achar(iachar(cchars(i)))
       end do
       call c_f_pointer(dump_path, cchars, [dump_path_len])
       dpath = ' '
       do i = 1, dump_path_len
          dpath(i:i) = achar(iachar(cchars(i)))
       end do

       open (newunit=du, status='scratch', action='readwrite', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_text_dump_c: cannot open scratch file for config'
          return
       end if
       call write_config_lines(du, ftext)
       call read_resolved_input(du)
       close (du)

       ! Read mesh (and, inside it, the enclosure + Neumann polylines),
       ! observation points, boundary conditions and wind — everything the
       ! Phase 6 parsers mirror — but do not run the model.
       call initialize_snapwave_domain()
       call read_obs_points()
       call read_boundary_data()
       call read_wind_data()

       open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_text_dump_c: cannot open dump file: ', trim(dpath)
          return
       end if

       call dump_text_globals(dunit)

       close (dunit)
       !
       status = 0_c_int
       !
    end function snapwave_text_dump_c

    subroutine dump_text_globals(u)
       !
       ! Dump the snapwave_data globals produced by the auxiliary text
       ! readers (obs points, boundary, wind, enclosure, neumann) in the
       ! format described above. Keep the field order in lock-step with
       ! src_rust/text_compare.rs.
       !
       use snapwave_data
       implicit none
       integer, intent(in) :: u
       integer :: i, ib, itb

       ! ---- observation points
       write (u, '(A)') 'section obs'
       write (u, '(A,1X,I0)') 'n', nobs
       if (nobs > 0) then
          write (u, '(A,1X,I0)') 'x', nobs
          do i = 1, nobs
             call dump_r8_line(u, xobs(i))
          end do
          write (u, '(A,1X,I0)') 'y', nobs
          do i = 1, nobs
             call dump_r8_line(u, yobs(i))
          end do
          write (u, '(A,1X,I0)') 'name', nobs
          do i = 1, nobs
             write (u, '(A)') trim(nameobs(i))
          end do
       end if

       ! ---- boundary conditions
       write (u, '(A)') 'section boundary'
       if (len_trim(jonswapfile) > 0) then
          ! Single-point JONSWAP file (read_boundary_data_singlepoint).
          write (u, '(A)') 'mode single'
          write (u, '(A,1X,I0)') 'nwbnd', nwbnd
          write (u, '(A,1X,I0)') 'ntwbnd', ntwbnd
          write (u, '(A,1X,I0)') 't', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, t_bwv(itb))
          end do
          write (u, '(A,1X,I0)') 'hs', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, hs_bwv(1, itb))
          end do
          write (u, '(A,1X,I0)') 'tp', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, tp_bwv(1, itb))
          end do
          write (u, '(A,1X,I0)') 'wd', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, wd_bwv(1, itb))
          end do
          write (u, '(A,1X,I0)') 'ds', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, ds_bwv(1, itb))
          end do
          write (u, '(A,1X,I0)') 'zs', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, zs_bwv(1, itb))
          end do
       else if (nwbnd > 0) then
          ! Space- and time-varying files (read_boundary_data_timeseries).
          write (u, '(A)') 'mode timeseries'
          write (u, '(A,1X,I0)') 'nwbnd', nwbnd
          write (u, '(A,1X,I0)') 'ntwbnd', ntwbnd
          write (u, '(A,1X,I0)') 'x', nwbnd
          do ib = 1, nwbnd
             call dump_r8_line(u, x_bwv(ib))
          end do
          write (u, '(A,1X,I0)') 'y', nwbnd
          do ib = 1, nwbnd
             call dump_r8_line(u, y_bwv(ib))
          end do
          write (u, '(A,1X,I0)') 't', ntwbnd
          do itb = 1, ntwbnd
             call dump_r4_line(u, t_bwv(itb))
          end do
          write (u, '(A,1X,I0)') 'hs', nwbnd * ntwbnd
          do itb = 1, ntwbnd
             do ib = 1, nwbnd
                call dump_r4_line(u, hs_bwv(ib, itb))
             end do
          end do
          write (u, '(A,1X,I0)') 'tp', nwbnd * ntwbnd
          do itb = 1, ntwbnd
             do ib = 1, nwbnd
                call dump_r4_line(u, tp_bwv(ib, itb))
             end do
          end do
          write (u, '(A,1X,I0)') 'wd', nwbnd * ntwbnd
          do itb = 1, ntwbnd
             do ib = 1, nwbnd
                call dump_r4_line(u, wd_bwv(ib, itb))
             end do
          end do
          write (u, '(A,1X,I0)') 'ds', nwbnd * ntwbnd
          do itb = 1, ntwbnd
             do ib = 1, nwbnd
                call dump_r4_line(u, ds_bwv(ib, itb))
             end do
          end do
          write (u, '(A,1X,I0)') 'zs', nwbnd * ntwbnd
          do itb = 1, ntwbnd
             do ib = 1, nwbnd
                call dump_r4_line(u, zs_bwv(ib, itb))
             end do
          end do
       else
          write (u, '(A)') 'mode none'
       end if

       ! ---- wind
       write (u, '(A)') 'section wind'
       write (u, '(A,1X,I0)') 'ntu10bnd', ntu10bnd
       if (len_trim(windlistfile) > 0) then
          write (u, '(A)') 'mode list'
          write (u, '(A,1X,I0)') 't', ntu10bnd
          do i = 1, ntu10bnd
             call dump_r4_line(u, t_u10_bwv(i))
          end do
       else
          write (u, '(A)') 'mode uniform'
          write (u, '(A,1X,I0)') 'u10', 1
          call dump_r4_line(u, u10_bwv(1, 1))
          write (u, '(A,1X,I0)') 'u10dir', 1
          call dump_r4_line(u, u10dir_bwv(1, 1))
       end if

       ! ---- boundary enclosure polyline
       write (u, '(A)') 'section enc'
       write (u, '(A,1X,I0)') 'n', n_bndenc
       if (n_bndenc > 0) then
          write (u, '(A,1X,I0)') 'x', n_bndenc
          do i = 1, n_bndenc
             call dump_r8_line(u, x_bndenc(i))
          end do
          write (u, '(A,1X,I0)') 'y', n_bndenc
          do i = 1, n_bndenc
             call dump_r8_line(u, y_bndenc(i))
          end do
       end if

       ! ---- neumann polyline
       write (u, '(A)') 'section neu'
       write (u, '(A,1X,I0)') 'n', n_neu
       if (n_neu > 0) then
          write (u, '(A,1X,I0)') 'x', n_neu
          do i = 1, n_neu
             call dump_r8_line(u, x_neu(i))
          end do
          write (u, '(A,1X,I0)') 'y', n_neu
          do i = 1, n_neu
             call dump_r8_line(u, y_neu(i))
          end do
       end if

       write (u, '(A)') 'section end'
       !
    end subroutine dump_text_globals

    subroutine dump_r8_line(u, v)
       !
       integer, intent(in) :: u
       real*8, intent(in) :: v
       integer(c_int64_t) :: bits
       bits = transfer(v, bits)
       write (u, '(Z16.16)') bits
    end subroutine dump_r8_line

    subroutine dump_r4_line(u, v)
       !
       integer, intent(in) :: u
       real*4, intent(in) :: v
       integer(c_int32_t) :: bits
       bits = transfer(v, bits)
       write (u, '(Z8.8)') bits
    end subroutine dump_r4_line

    !************************************************************************
    ! plan.md Phase 9 hook: compute the *derived geometry* — surrounding
    ! points and upwind neighbours (initialize_snapwave_domain), observation
    ! interpolation weights (read_obs_points -> make_map_fm) and boundary
    ! support-point mapping (read_boundary_data -> find_boundary_indices) —
    ! with the unchanged Fortran routines and dump the resulting globals, so
    ! the Rust ports (src_rust/geometry.rs, src_rust/interp.rs) can be pinned
    ! against the numerical oracle. Mirrors snapwave_text_dump_c but stops
    ! before read_wind_data (wind is not geometry) and dumps the geometry
    ! globals instead of the parsed-text globals.
    !
    ! Dump format (parsed by src_rust/geometry_compare.rs): `section <name>`
    ! blocks of `key value` lines. Array keys (kp, dhdx, dhdy, w360, prev360,
    ! ds360, msk, neumannconnected, nmindbnd, neubnd, wobs, irefobs, nrefobs,
    ! ind1_bwv_cst, ind2_bwv_cst, fac_bwv_cst) carry a count and then one
    ! value per line: reals as IEEE-754 bit patterns (real*8 Z16.16, real*4
    ! Z8.8), integers decimal. Scalar keys (no_nodes, np, ntheta360, nb,
    ! nnmb, nobs, nwbnd) carry a single value. Array order is Fortran
    ! column-major over the natural do-loop order; see dump_geometry_globals.
    !************************************************************************
    function snapwave_geometry_dump_c(config, config_len, dump_path, dump_path_len) &
          bind(C, name="snapwave_geometry_dump_c") result(status)
       use snapwave_data
       use snapwave_input
       use snapwave_domain
       use snapwave_boundaries
       use snapwave_obspoints
       implicit none

       type(c_ptr), value :: config
       integer(c_int), value :: config_len
       type(c_ptr), value :: dump_path
       integer(c_int), value :: dump_path_len
       integer(c_int) :: status

       character(len=:), allocatable :: ftext
       character(len=1024) :: dpath
       character(kind=c_char), dimension(:), pointer :: cchars
       integer :: i, ios, dunit, du

       status = 1_c_int

       allocate (character(len=config_len) :: ftext)
       call c_f_pointer(config, cchars, [config_len])
       do i = 1, config_len
          ftext(i:i) = achar(iachar(cchars(i)))
       end do
       call c_f_pointer(dump_path, cchars, [dump_path_len])
       dpath = ' '
       do i = 1, dump_path_len
          dpath(i:i) = achar(iachar(cchars(i)))
       end do

       open (newunit=du, status='scratch', action='readwrite', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_geometry_dump_c: cannot open scratch file for config'
          return
       end if
       call write_config_lines(du, ftext)
       call read_resolved_input(du)
       close (du)

       ! Compute the geometry: mesh + enclosure/neumann (surrounding points,
       ! upwind neighbours, mask refinement), observation interpolation
       ! weights, and boundary support-point mapping.
       call initialize_snapwave_domain()
       call read_obs_points()
       call read_boundary_data()

       open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
       if (ios /= 0) then
          write (*, *) 'ERROR: snapwave_geometry_dump_c: cannot open dump file: ', trim(dpath)
          return
       end if

       call dump_geometry_globals(dunit)

       close (dunit)
       !
       status = 0_c_int
       !
    end function snapwave_geometry_dump_c

    subroutine dump_geometry_globals(u)
       !
       ! Dump the derived-geometry globals produced by
       ! initialize_snapwave_domain + read_obs_points + read_boundary_data,
       ! in the format described above. Keep the field/array order in
       ! lock-step with src_rust/geometry_compare.rs.
       !
       use snapwave_data
       implicit none
       integer, intent(in) :: u
       integer :: k, ip, itheta

       write (u, '(A)') 'section domain'
       write (u, '(A,1X,I0)') 'no_nodes', no_nodes
       write (u, '(A,1X,I0)') 'np', np
       write (u, '(A,1X,I0)') 'ntheta360', ntheta360
       write (u, '(A,1X,I0)') 'nb', nb
       write (u, '(A,1X,I0)') 'nnmb', nnmb
       write (u, '(A,1X,I0)') 'kp', np*no_nodes
       do k = 1, no_nodes
          do ip = 1, np
             write (u, '(I0)') kp(ip, k)
          end do
       end do
       write (u, '(A,1X,I0)') 'dhdx', no_nodes
       do k = 1, no_nodes
          call dump_r4_line(u, dhdx(k))
       end do
       write (u, '(A,1X,I0)') 'dhdy', no_nodes
       do k = 1, no_nodes
          call dump_r4_line(u, dhdy(k))
       end do
       write (u, '(A,1X,I0)') 'w360', 2*ntheta360*no_nodes
       do k = 1, no_nodes
          do itheta = 1, ntheta360
             do ip = 1, 2
                call dump_r4_line(u, w360(ip, itheta, k))
             end do
          end do
       end do
       write (u, '(A,1X,I0)') 'prev360', 2*ntheta360*no_nodes
       do k = 1, no_nodes
          do itheta = 1, ntheta360
             do ip = 1, 2
                write (u, '(I0)') prev360(ip, itheta, k)
             end do
          end do
       end do
       write (u, '(A,1X,I0)') 'ds360', ntheta360*no_nodes
       do k = 1, no_nodes
          do itheta = 1, ntheta360
             call dump_r4_line(u, ds360(itheta, k))
          end do
       end do
       write (u, '(A,1X,I0)') 'msk', no_nodes
       do k = 1, no_nodes
          write (u, '(I0)') msk(k)
       end do
       write (u, '(A,1X,I0)') 'neumannconnected', no_nodes
       do k = 1, no_nodes
          write (u, '(I0)') neumannconnected(k)
       end do
       write (u, '(A,1X,I0)') 'nmindbnd', nb
       do k = 1, nb
          write (u, '(I0)') nmindbnd(k)
       end do
       if (nnmb > 0) then
          write (u, '(A,1X,I0)') 'neubnd', nnmb
          do k = 1, nnmb
             write (u, '(I0)') neubnd(k)
          end do
       end if

       write (u, '(A)') 'section obs'
       write (u, '(A,1X,I0)') 'nobs', nobs
       if (nobs > 0) then
          write (u, '(A,1X,I0)') 'wobs', 4*nobs
          do k = 1, nobs
             do ip = 1, 4
                call dump_r8_line(u, wobs(ip, k))
             end do
          end do
          write (u, '(A,1X,I0)') 'irefobs', 4*nobs
          do k = 1, nobs
             do ip = 1, 4
                write (u, '(I0)') irefobs(ip, k)
             end do
          end do
          write (u, '(A,1X,I0)') 'nrefobs', no_nodes
          do k = 1, no_nodes
             write (u, '(I0)') nrefobs(k)
          end do
       end if

       write (u, '(A)') 'section boundary'
       write (u, '(A,1X,I0)') 'nwbnd', nwbnd
       write (u, '(A,1X,I0)') 'nb', nb
       if (nwbnd > 0 .and. nb > 0) then
          write (u, '(A,1X,I0)') 'ind1_bwv_cst', nb
          do k = 1, nb
             write (u, '(I0)') ind1_bwv_cst(k)
          end do
          write (u, '(A,1X,I0)') 'ind2_bwv_cst', nb
          do k = 1, nb
             write (u, '(I0)') ind2_bwv_cst(k)
          end do
          write (u, '(A,1X,I0)') 'fac_bwv_cst', nb
          do k = 1, nb
             call dump_r4_line(u, fac_bwv_cst(k))
          end do
       end if

write (u, '(A)') 'section end'
        !
     end subroutine dump_geometry_globals

     !************************************************************************
     ! plan.md Phase 11 hook: run the unchanged Fortran solver for one
     ! timestep and dump the resulting solver-state globals, so the Rust
     ! ports (src_rust/solver.rs) can be pinned against the numerical oracle.
     !
     ! Dump format (parsed by src_rust/solver_compare.rs): `section <name>`
     ! blocks of `key value` lines. Array keys carry a count and then one
     ! value per line: reals as IEEE-754 bit patterns (real*4 Z8.8),
     ! integers decimal. Scalar keys carry a single value. Array order is
     ! Fortran column-major over the natural do-loop order.
     !************************************************************************
     function snapwave_solver_dump_c(config, config_len, dump_path, dump_path_len) &
           bind(C, name="snapwave_solver_dump_c") result(status)
        use snapwave_data
        use snapwave_input
        use snapwave_domain
        use snapwave_boundaries
        use snapwave_obspoints
        use snapwave_solver
        use snapwave_ncoutput
        implicit none

        type(c_ptr), value :: config
        integer(c_int), value :: config_len
        type(c_ptr), value :: dump_path
        integer(c_int), value :: dump_path_len
        integer(c_int) :: status

        character(len=:), allocatable :: ftext
        character(len=1024) :: dpath
        character(kind=c_char), dimension(:), pointer :: cchars
        integer :: i, ios, dunit, du

        status = 1_c_int

        allocate (character(len=config_len) :: ftext)
        call c_f_pointer(config, cchars, [config_len])
        do i = 1, config_len
           ftext(i:i) = achar(iachar(cchars(i)))
        end do
        call c_f_pointer(dump_path, cchars, [dump_path_len])
        dpath = ' '
        do i = 1, dump_path_len
           dpath(i:i) = achar(iachar(cchars(i)))
        end do

        open (newunit=du, status='scratch', action='readwrite', iostat=ios)
        if (ios /= 0) then
           write (*, *) 'ERROR: snapwave_solver_dump_c: cannot open scratch file for config'
           return
        end if
        call write_config_lines(du, ftext)
        call read_resolved_input(du)
        close (du)

        ! Initialize domain, read obs/boundary/wind, then run one timestep.
        call initialize_snapwave_domain()
        call read_obs_points()
        call read_boundary_data()
        call read_wind_data()

        ! Initialize output (needed for some globals like map_filename)
        call ncoutput_init()

        ! Update boundary conditions for t=tstart and run one solver step
        call update_boundary_conditions(tstart)
        call compute_wave_field(tstart)

        open (newunit=dunit, file=trim(dpath), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
           write (*, *) 'ERROR: snapwave_solver_dump_c: cannot open dump file: ', trim(dpath)
           return
        end if

        call dump_solver_globals(dunit)

        close (dunit)
        !
        status = 0_c_int
        !
     end function snapwave_solver_dump_c

     subroutine dump_solver_globals(u)
        !
        ! Dump the solver-state globals produced by one call to
        ! compute_wave_field, in the format described above. Keep the
        ! field/array order in lock-step with src_rust/solver_compare.rs.
        !
        use snapwave_data
        implicit none
        integer, intent(in) :: u
        integer :: k, itheta

        write (u, '(A)') 'section solver'
        write (u, '(A,1X,I0)') 'no_nodes', no_nodes
        write (u, '(A,1X,I0)') 'ntheta', ntheta
        write (u, '(A,1X,I0)') 'ig', ig
        write (u, '(A,1X,I0)') 'wind', wind

        write (u, '(A,1X,I0)') 'H', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, H(k))
        end do
        write (u, '(A,1X,I0)') 'Dw', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, Dw(k))
        end do
        write (u, '(A,1X,I0)') 'Df', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, Df(k))
        end do
        write (u, '(A,1X,I0)') 'F', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, F(k))
        end do
        write (u, '(A,1X,I0)') 'thetam', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, thetam(k))
        end do
        write (u, '(A,1X,I0)') 'Tp', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, Tp(k))
        end do
        write (u, '(A,1X,I0)') 'sig', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, sig(k))
        end do
        write (u, '(A,1X,I0)') 'kwav', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, kwav(k))
        end do
        write (u, '(A,1X,I0)') 'cg', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, Cg(k))
        end do
        write (u, '(A,1X,I0)') 'sinhkh', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, sinhkh(k))
        end do
        write (u, '(A,1X,I0)') 'Hmx', no_nodes
        do k = 1, no_nodes
           call dump_r4_line(u, Hmx(k))
        end do

        ! Directional energy density (first 5 nodes to keep dump manageable)
        write (u, '(A,1X,I0)') 'ee_nodes', min(5, no_nodes)
        do k = 1, min(5, no_nodes)
           write (u, '(A,1X,I0)') 'ee', ntheta
           do itheta = 1, ntheta
              call dump_r4_line(u, ee(itheta, k))
           end do
        end do

        if (ig == 1) then
           write (u, '(A,1X,I0)') 'H_ig', no_nodes
           do k = 1, no_nodes
              call dump_r4_line(u, H_ig(k))
           end do
        end if

        if (wind == 1) then
           write (u, '(A,1X,I0)') 'SwE', no_nodes
           do k = 1, no_nodes
              call dump_r4_line(u, SwE(k))
           end do
           write (u, '(A,1X,I0)') 'SwA', no_nodes
           do k = 1, no_nodes
              call dump_r4_line(u, SwA(k))
           end do
        end if

        write (u, '(A)') 'section end'
        !
     end subroutine dump_solver_globals
 end module snapwave_c_api
