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
      real*8  :: t
      real*8  :: output_tol
      integer :: it
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
       call ncoutput_finalize()
       !
       status = 0_c_int
       !
    end function snapwave_run_c

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
end module snapwave_c_api
