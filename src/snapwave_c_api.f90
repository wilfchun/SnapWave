module snapwave_c_api
   !************************************************************************
   ! Coarse C ABI facade over the existing SnapWave Fortran model.
   !
   ! Exposes a single bind(C) entry point, snapwave_run_c, that mirrors the
   ! lifecycle of the stand-alone program in src/snapwave.f90:
   !    1. read input
   !    2. initialize domain
   !    3. read observation points
   !    4. read boundary conditions (incl. wind if specified)
   !    5. initialize NetCDF output
   !    6. run timesteps
   !    7. finalize output
   !
   ! The solver internals and the module global state are unchanged; the
   ! Rust wrapper provides main() and calls this facade. The input file is
   ! expected to be found in the current working directory (read_snapwave_input
   ! probes snapwave.inp/SnapWave.inp/... itself), so the caller is expected
   ! to chdir to the input file's parent directory and pass only the file
   ! name. The path is validated here as well before the model runs.
   !************************************************************************
   use iso_c_binding
   implicit none
contains
   function snapwave_run_c(input_path, input_path_len) bind(C, name="snapwave_run_c") result(status)
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

      type(c_ptr), value :: input_path
      integer(c_int), value :: input_path_len
      integer(c_int) :: status

      character(len=1024) :: fpath
      character(kind=c_char), dimension(:), pointer :: cchars
      real*8  :: t
      real*8  :: output_tol
      integer :: it
      integer :: i
      logical :: exists

      ! Convert C string/path to Fortran character (not NUL terminated;
      ! the length is passed explicitly).
      call c_f_pointer(input_path, cchars, [input_path_len])
      fpath = ' '
      do i = 1, input_path_len
         fpath(i:i) = achar(iachar(cchars(i)))
      end do

      inquire (file=trim(fpath), exist=exists)
      if (.not. exists) then
         write (*, *) 'ERROR: snapwave_c_api: input file not found: ', trim(fpath)
         status = 1_c_int
         return
      end if

      call read_snapwave_input()            ! Reads snapwave.inp
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
end module snapwave_c_api
