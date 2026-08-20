module snapwave_input
   implicit none
contains

   subroutine read_snapwave_input(input_file, status)
      !
      ! Reads snapwave.inp
      !
      ! The optional input_file argument (plan.md Phase 3, step 6) lets the
      ! C ABI facade pass the Rust-selected input file explicitly. Without
      ! it, the legacy stand-alone behaviour is preserved: probe the usual
      ! file names in the current working directory.
      !
      ! The optional status argument (plan.md Phase 4, step 5) converts the
      ! hard `stop 1` error paths into status returns when the caller is the
      ! C ABI facade, so a bad configuration becomes a clean Rust-side
      ! failure instead of killing the whole process. The legacy stand-alone
      ! program omits status and keeps the original `stop` behaviour.
      !
      use snapwave_data
      use snapwave_date
      !
      implicit none
      !
      character(len=*), optional, intent(in) :: input_file
      integer, optional, intent(out) :: status
      !
      integer :: dtsec
      integer :: irestart
      integer :: iwritetestfiles
      !
      character(len=1024) :: filename
      integer :: ios, ii
      logical :: exists

! List of possible reasonable filenames
      character(len=*), parameter :: possible_names(*) = [ &
                                     'snapwave.inp', 'SnapWave.inp', 'SNAPWAVE.INP', &
                                     'snapwave.INP', 'Snapwave.INP', 'SNAPWAVE.inp']

      if (present(status)) status = 0

      write (*, *) 'Reading input file ...'

      if (present(input_file)) then
         !
         ! Rust-selected input file: no probing; the wrapper has already
         ! validated and selected this exact file (plan.md Phase 3).
         !
         filename = trim(input_file)
         inquire (file=trim(filename), exist=exists, iostat=ios)
         if (.not. (exists .and. ios == 0)) then
            write (*, *) 'ERROR: input file not found: ', trim(filename)
            if (present(status)) then
               status = 1
               return
            end if
            stop 1
         end if
      else
         !
         ! Legacy stand-alone behaviour: probe the usual names in the CWD.
         !
         do ii = 1, size(possible_names)
            filename = trim(possible_names(ii))
            inquire (file=filename, exist=exists, iostat=ios)
            if (exists .and. ios == 0) exit
         end do

         if (.not. exists) then
            write (*, *) 'ERROR: none of the expected input files were found:'
            do ii = 1, size(possible_names)
               write (*, *) '   - ', trim(possible_names(ii))
            end do
            if (present(status)) then
               status = 1
               return
            end if
            stop 1
         end if
      end if

      open (unit=500, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'ERROR: found file but failed to open: ', trim(filename)
         if (present(status)) then
            status = 1
            return
         end if
         stop 1
      end if

      write (*, *) 'Successfully opened: ', trim(filename)
      !
      ! Input section
      !
      call read_int_input(500, 'nmax', nmax, 0)
      call read_int_input(500, 'mmax', mmax, 0)
      call read_real_input(500, 'dx', dx, 0.0)
      call read_real_input(500, 'dy', dy, 0.0)
      call read_real_input(500, 'x0', x0, 0.0)
      call read_real_input(500, 'y0', y0, 0.0)
      call read_real_input(500, 'rotation', rotation, 0.0)
      call read_real_input(500, 'posdwn', posdwn, -1.0)
      call read_char_input(500, 'tref', trefstr, '20000101 000000')
      call read_char_input(500, 'tstart', tstartstr, '20000101 000000')
      call read_char_input(500, 'tstop', tstopstr, '20000101 000000')
      call read_real_input(500, 'timestep', timestep, 3600.0)
      call read_int_input(500, 'niter', niter, 10)
      call read_real_input(500, 'crit', crit, 0.00001)
      call read_real_input(500, 'dt', dt, 36000.0)
      call read_real_input(500, 'gamma', gamma, 0.7)
      call read_real_input(500, 'alpha', alpha, 1.0)
      call read_real_input(500, 'hmin', hmin, 0.1)
      call read_real_input(500, 'gammax', gammax, 0.6)
      call read_char_input(500, 'gridfile', gridfile, '.txt')
      call read_int_input(500, 'sferic', sferic, 0)
      call read_char_input(500, 'fw', fwstr, '0.01')
      call read_char_input(500, 'fwig', fw_igstr, '0.015')
      call read_real_input(500, 'fwcutoff', fwcutoff, 200.0)
      call read_real_input(500, 'tol', tol, 10.0)
      call read_real_input(500, 'dtheta', dtheta, 10.0)
      call read_real_input(500, 'sector', sector, 180.0)
      call read_char_input(500, 'jonswapfile', jonswapfile, '')
      call read_char_input(500, 'windlistfile', windlistfile, '')
      call read_char_input(500, 'bndfile', bndfile, 'none')
      call read_char_input(500, 'encfile', encfile, 'none')
      call read_char_input(500, 'neumannfile', neumannfile, 'none')
      call read_char_input(500, 'bhsfile', bhsfile, '')
      call read_char_input(500, 'btpfile', btpfile, '')
      call read_char_input(500, 'bwdfile', bwdfile, '')
      call read_char_input(500, 'bdsfile', bdsfile, '')
      call read_char_input(500, 'bzsfile', bzsfile, '')
      call read_char_input(500, 'upwfile', upwfile, '')
      call read_char_input(500, 'mskfile', mskfile, '')
      call read_char_input(500, 'indfile', indfile, '')
      call read_char_input(500, 'depfile', depfile, '')
      call read_char_input(500, 'obsfile', obsfile, 'none')
      call read_char_input(500, 'map_file', map_filename, '')
      call read_char_input(500, 'his_file', his_filename, '')
      call read_real_input(500, 'map_interval', map_interval, timestep)
      call read_real_input(500, 'his_interval', his_interval, timestep)
      if (map_filename /= '' .and. map_interval <= 0.0) then
         write (*, *) 'ERROR: map_interval must be positive.'
         if (present(status)) then
            status = 1
            return
         end if
         stop 1
      end if
      if (his_filename /= '' .and. his_interval <= 0.0) then
         write (*, *) 'ERROR: his_interval must be positive.'
         if (present(status)) then
            status = 1
            return
         end if
         stop 1
      end if
      call read_int_input(500, 'map_depth', map_dep, 1)
      call read_int_input(500, 'map_Hm0', map_Hm0, 1)
      call read_int_input(500, 'map_Hig', map_Hig, 0)
      call read_int_input(500, 'map_Tp', map_Tp, 1)
      call read_int_input(500, 'map_dir', map_dir, 1)
      call read_int_input(500, 'map_dirspr', map_dirspr, 0)
      call read_int_input(500, 'map_Cg', map_Cg, 0)
      call read_int_input(500, 'map_Dw', map_Dw, 0)
      call read_int_input(500, 'map_Df', map_Df, 0)
      call read_int_input(500, 'map_SwE', map_SwE, 0)
      call read_int_input(500, 'map_SwA', map_SwA, 0)
      call read_int_input(500, 'map_sig', map_sig, 0)
      call read_int_input(500, 'map_u10', map_u10, 0)
      call read_int_input(500, 'map_Dveg', map_Dveg, 0)
      call read_int_input(500, 'writetestfiles', iwritetestfiles, 0)
      call read_int_input(500, 'ja_save_each_iter', ja_save_each_iter, 0)

      call read_int_input(500, 'map_ee', map_ee, 0)
      call read_int_input(500, 'map_ctheta', map_ctheta, 0)
      call read_int_input(500, 'restart', irestart, 0)
      !
      call read_char_input(500, 'u10', u10str, '0.0')
      call read_char_input(500, 'u10dir', u10dirstr, '270.0')
      call read_real_input(500, 'Tpini', Tpini, 1.0)
      call read_int_input(500, 'mwind', mwind, 2)
      call read_real_input(500, 'sigmin', sigmin, 8.0 * atan(1.0) / 25.0)
      call read_real_input(500, 'sigmax', sigmax, 8.0 * atan(1.0) / 1.0)
      call read_int_input(500, 'jadcgdx', jadcgdx, 1)
      call read_real_input(500, 'c_dispT', c_dispT, 1.0)
      call read_real_input(500, 'zsini', zsini, 0.0)
      call read_int_input(500, 'ig', ig, 0)
      call read_int_input(500, 'upwindref', upwindref, 0)
      !
      ! Vegetation input
      !
      call read_int_input(500, 'ja_vegetation', ja_vegetation, 0)
      call read_char_input(500, 'vegmapfile', vegmapfile, '.txt')
      !
      wind = 0
      if ((u10str == '0.0') .and. (windlistfile == '')) then
         !
         wind = 0
         !
         write (*, *) '   Uniform wave period in entire domain.'
         !
      else
         !
         write (*, *) '   Wind growth turned on.'
         wind = 1
      end if
      !
      close (500)
      !
      call time_difference(trefstr, tstartstr, dtsec) ! time difference in seconds between tstart and tref
      tstart = dtsec * 1.0 ! time difference in seconds between tstop and tstart
      call time_difference(trefstr, tstopstr, dtsec)
      tstop = dtsec * 1.0 ! time difference in seconds between tstop and tstart
      !
      mmax = mmax + 2 ! Original mmax and nmax are for number of cells in bathy grid. Add two dummy rows.
      nmax = nmax + 2
      !
      restart = .true.
      writetestfiles = .true.
      if (irestart == 0) restart = .false.
      if (iwritetestfiles == 0) writetestfiles = .false.
      !
   end subroutine

   subroutine read_real_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      real*4, intent(out) :: value
      real*4, intent(in) :: default
      integer j, stat
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) value
            exit
         end if
      end do
      !
   end subroutine

   subroutine read_real_array_input(fileid, keyword, value, default, nr)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      integer, intent(in) :: nr
      real*4, dimension(:), intent(out), allocatable :: value
      real*4, intent(in) :: default
      integer j, stat, m
      !
      allocate (value(nr))
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) (value(m), m=1, nr)
            exit
         end if
      end do
      !
   end subroutine

   subroutine read_int_input(fileid, keyword, value, default)
      !
      character(*), intent(in) :: keyword
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer, intent(in) :: fileid
      integer, intent(out) :: value
      integer, intent(in) :: default
      integer j, stat
      !
      value = default
      rewind (fileid)
      do while (.true.)
         read (fileid, '(a)', iostat=stat) line
         if (stat < 0) exit
         j = index(line, '=')
         keystr = trim(line(1:j - 1))
         if (trim(keystr) == trim(keyword)) then
            valstr = trim(line(j + 1:256))
            read (valstr, *) value
            exit
         end if
      end do
      !
   end subroutine

    subroutine read_char_input(fileid, keyword, value, default)
       !
       character(*), intent(in) :: keyword
       character(len=256) :: keystr
       character(len=256) :: valstr
       character(len=256) :: line
       integer, intent(in) :: fileid
       character(*), intent(in) :: default
       character(*), intent(out) :: value
       integer j, stat
       !
       value = default
       rewind (fileid)
       do while (.true.)
          read (fileid, '(a)', iostat=stat) line
          if (stat < 0) exit
          j = index(line, '=')
          keystr = trim(line(1:j - 1))
          if (trim(keystr) == trim(keyword)) then
             valstr = adjustl(trim(line(j + 1:256)))
             value = valstr
             exit
          end if
       end do
       !
    end subroutine

    subroutine read_real8_input(fileid, keyword, value, default)
       !
       ! Same keyword-scanning as read_real_input, but for real*8 globals
       ! (tstart/tstop). Needed by read_resolved_input (plan.md Phase 4),
       ! which receives those two already computed by the Rust parser.
       !
       character(*), intent(in) :: keyword
       character(len=256) :: keystr
       character(len=256) :: valstr
       character(len=256) :: line
       integer, intent(in) :: fileid
       real*8, intent(out) :: value
       real*8, intent(in) :: default
       integer j, stat
       !
       value = default
       rewind (fileid)
       do while (.true.)
          read (fileid, '(a)', iostat=stat) line
          if (stat < 0) exit
          j = index(line, '=')
          keystr = trim(line(1:j - 1))
          if (trim(keystr) == trim(keyword)) then
             valstr = trim(line(j + 1:256))
             read (valstr, *) value
             exit
          end if
       end do
       !
    end subroutine

    subroutine write_config_lines(unit, text)
       !
       ! Write the Rust-resolved configuration text (one `key=value` per
       ! line, newline-separated) as records of a scratch file, so the
       ! existing read_*_input keyword scanners can consume it
       ! (plan.md Phase 4). len_trim drops only trailing blanks, which the
       ! resolved text never carries (values are already trimmed).
       !
       integer, intent(in) :: unit
       character(len=*), intent(in) :: text
       integer :: i, n, start
       !
       n = len_trim(text)
       start = 1
       do i = 1, n
          if (text(i:i) == achar(10)) then
             if (i > start) write (unit, '(a)') text(start:i - 1)
             start = i + 1
          end if
       end do
       if (start <= n) write (unit, '(a)') text(start:n)
       !
    end subroutine

    subroutine read_resolved_input(unit)
       !
       ! plan.md Phase 4: consume the fully-resolved configuration produced
       ! by the Rust parser (every key already defaulted, validated and
       ! post-processed) and store it into the snapwave_data globals. On the
       ! Rust route Fortran no longer reads SnapWave.inp nor decides
       ! defaults; it only receives the resolved values.
       !
       ! No defaults, no post-processing (mmax/nmax are already the +2
       ! model-facing values, tstart/tstop are already seconds, restart /
       ! wind / writetestfiles are already booleans) and no `stop` paths:
       ! Rust is the authority here. All keys are always present, so the
       ! default arguments of the read_*_input helpers never fire.
       !
       use snapwave_data
       implicit none
       !
       integer, intent(in) :: unit
       integer :: itemp
       !
       ! ---- time / control
       call read_char_input(unit, 'trefstr', trefstr, '20000101 000000')
       call read_char_input(unit, 'tstartstr', tstartstr, '20000101 000000')
       call read_char_input(unit, 'tstopstr', tstopstr, '20000101 000000')
       call read_real8_input(unit, 'tstart', tstart, 0.0d0)
       call read_real8_input(unit, 'tstop', tstop, 0.0d0)
       call read_real_input(unit, 'timestep', timestep, 3600.0)
       call read_real_input(unit, 'dt', dt, 36000.0)
       call read_int_input(unit, 'niter', niter, 10)
       call read_real_input(unit, 'crit', crit, 0.00001)
       call read_int_input(unit, 'restart', itemp, 0)
       restart = (itemp /= 0)
       !
       ! ---- grid / domain (mmax/nmax already include the +2 dummy rows)
       call read_int_input(unit, 'mmax', mmax, 0)
       call read_int_input(unit, 'nmax', nmax, 0)
       call read_real_input(unit, 'dx', dx, 0.0)
       call read_real_input(unit, 'dy', dy, 0.0)
       call read_real_input(unit, 'x0', x0, 0.0)
       call read_real_input(unit, 'y0', y0, 0.0)
       call read_real_input(unit, 'rotation', rotation, 0.0)
       call read_real_input(unit, 'posdwn', posdwn, -1.0)
       call read_int_input(unit, 'sferic', sferic, 0)
       call read_real_input(unit, 'dtheta', dtheta, 10.0)
       call read_real_input(unit, 'sector', sector, 180.0)
       call read_char_input(unit, 'gridfile', gridfile, '.txt')
       call read_char_input(unit, 'depfile', depfile, '')
       call read_char_input(unit, 'mskfile', mskfile, '')
       call read_char_input(unit, 'indfile', indfile, '')
       call read_char_input(unit, 'upwfile', upwfile, '')
       !
       ! ---- boundary forcing
       call read_char_input(unit, 'jonswapfile', jonswapfile, '')
       call read_char_input(unit, 'bndfile', bndfile, 'none')
       call read_char_input(unit, 'encfile', encfile, 'none')
       call read_char_input(unit, 'neumannfile', neumannfile, 'none')
       call read_char_input(unit, 'bhsfile', bhsfile, '')
       call read_char_input(unit, 'btpfile', btpfile, '')
       call read_char_input(unit, 'bwdfile', bwdfile, '')
       call read_char_input(unit, 'bdsfile', bdsfile, '')
       call read_char_input(unit, 'bzsfile', bzsfile, '')
       call read_char_input(unit, 'obsfile', obsfile, 'none')
       call read_real_input(unit, 'tol', tol, 10.0)
       !
       ! ---- wind
       call read_char_input(unit, 'u10str', u10str, '0.0')
       call read_char_input(unit, 'u10dirstr', u10dirstr, '270.0')
       call read_char_input(unit, 'windlistfile', windlistfile, '')
       call read_int_input(unit, 'mwind', mwind, 2)
       call read_int_input(unit, 'wind', wind, 0)
       !
       ! ---- output
       call read_char_input(unit, 'map_filename', map_filename, '')
       call read_char_input(unit, 'his_filename', his_filename, '')
       call read_real_input(unit, 'map_interval', map_interval, 3600.0)
       call read_real_input(unit, 'his_interval', his_interval, 3600.0)
       call read_int_input(unit, 'map_dep', map_dep, 1)
       call read_int_input(unit, 'map_Hm0', map_Hm0, 1)
       call read_int_input(unit, 'map_Hig', map_Hig, 0)
       call read_int_input(unit, 'map_Tp', map_Tp, 1)
       call read_int_input(unit, 'map_dir', map_dir, 1)
       call read_int_input(unit, 'map_dirspr', map_dirspr, 0)
       call read_int_input(unit, 'map_cg', map_cg, 0)
       call read_int_input(unit, 'map_Dw', map_Dw, 0)
       call read_int_input(unit, 'map_Df', map_Df, 0)
       call read_int_input(unit, 'map_SwE', map_SwE, 0)
       call read_int_input(unit, 'map_SwA', map_SwA, 0)
       call read_int_input(unit, 'map_sig', map_sig, 0)
       call read_int_input(unit, 'map_u10', map_u10, 0)
       call read_int_input(unit, 'map_Dveg', map_Dveg, 0)
       call read_int_input(unit, 'map_ee', map_ee, 0)
       call read_int_input(unit, 'map_ctheta', map_ctheta, 0)
       call read_int_input(unit, 'ja_save_each_iter', ja_save_each_iter, 0)
       !
       ! ---- diagnostics
       call read_int_input(unit, 'writetestfiles', itemp, 0)
       writetestfiles = (itemp /= 0)
       !
       ! ---- vegetation
       call read_int_input(unit, 'ja_vegetation', ja_vegetation, 0)
       call read_char_input(unit, 'vegmapfile', vegmapfile, '.txt')
       !
       ! ---- solver physics knobs
       call read_real_input(unit, 'gamma', gamma, 0.7)
       call read_real_input(unit, 'alpha', alpha, 1.0)
       call read_real_input(unit, 'gammax', gammax, 0.6)
       call read_real_input(unit, 'hmin', hmin, 0.1)
       call read_real_input(unit, 'fwcutoff', fwcutoff, 200.0)
       call read_char_input(unit, 'fwstr', fwstr, '0.01')
       call read_char_input(unit, 'fw_igstr', fw_igstr, '0.015')
       call read_real_input(unit, 'Tpini', Tpini, 1.0)
       call read_real_input(unit, 'zsini', zsini, 0.0)
       call read_real_input(unit, 'sigmin', sigmin, 0.0)
       call read_real_input(unit, 'sigmax', sigmax, 0.0)
       call read_int_input(unit, 'jadcgdx', jadcgdx, 1)
       call read_real_input(unit, 'c_dispT', c_dispT, 1.0)
       call read_int_input(unit, 'ig', ig, 0)
       call read_int_input(unit, 'upwindref', upwindref, 0)
       !
    end subroutine

 end module snapwave_input
