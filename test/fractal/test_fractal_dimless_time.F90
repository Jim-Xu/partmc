! Copyright (C) 2009-2012, 2017 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The test_fractal_dimless_time program.

!> Read NetCDF output files and write out the time evolution of aerosol
!> dimensionless number concentrations in text format.

program test_fractal_dimless_time

  use pmc_aero_state
  use pmc_aero_particle
  use pmc_output
  use pmc_mpi
  use pmc_util
  use pmc_constants
  use getopt_m

  integer, parameter :: REGIME_FREE = 1
  integer, parameter :: REGIME_CONT = 2

  character(len=PMC_MAX_FILENAME_LEN) :: in_prefix, out_filename
  character(len=PMC_MAX_FILENAME_LEN), allocatable :: filename_list(:)
  character(len=1000) :: tmp_str
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: index, i_repeat, i_spec, out_unit
  integer :: i_file, n_file, regime
  real(kind=dp) :: time, del_t, density, n_init
  type(aero_particle_t), pointer :: aero_particle
  character(len=PMC_UUID_LEN) :: uuid, run_uuid
  real(kind=dp), allocatable :: particle_num_concs(:), particle_masses(:)
  real(kind=dp), allocatable :: particle_volumes(:)
  real(kind=dp), allocatable :: times(:), time_num_concs(:)
  real(kind=dp), allocatable :: dimless_times(:), dimless_time_num_concs(:)
  type(option_s) :: opts(3)

  call pmc_mpi_init()

  opts(1) = option_s("free", .false., 'f')
  opts(2) = option_s("cont", .false., 'c')
  opts(3) = option_s("output", .true., 'o')

  regime = 0
  n_init = 1d14
  out_filename = ""

  do
     select case(getopt("fco:", opts))
     case(char(0))
        exit
     case('f')
        regime = REGIME_FREE
     case('c')
        regime = REGIME_CONT
     case('o')
        out_filename = optarg
     end select
  end do

  if (regime == 0) then
     call die_msg(738176431, 'missing aerosol regime: ' &
          // 'please specify --free or --cont')
  end if

  if (optind /= command_argument_count()) then
     call die_msg(467447497, &
          'expected exactly one non-option prefix argument')
  end if

  call get_command_argument(optind, in_prefix)

  if (out_filename == "") then
     out_filename = trim(in_prefix) // "_dimless_time.txt"
  end if

  allocate(filename_list(0))
  call input_filename_list(in_prefix, filename_list)
  n_file = size(filename_list)
  call assert_msg(582892256, n_file > 0, &
       "no NetCDF files found with prefix: " // trim(in_prefix))

  call input_state(filename_list(1), index, time, del_t, &
       i_repeat, uuid, aero_data=aero_data, &
       aero_state=aero_state, env_state=env_state)
  run_uuid = uuid

  allocate(times(n_file))
  allocate(dimless_times(n_file))
  allocate(time_num_concs(n_file))
  allocate(dimless_time_num_concs(n_file))

  allocate(particle_num_concs(0))
  allocate(particle_masses(0))
  allocate(particle_volumes(0))

  do i_file = 1,n_file
     call input_state(filename_list(i_file), index, time, &
          del_t, i_repeat, uuid, aero_data=aero_data, &
          aero_state=aero_state, env_state=env_state)

     call assert_msg(545569193, uuid == run_uuid, &
          "UUID mismatch between " // trim(filename_list(1)) // " and " &
          // trim(filename_list(i_file)))

     times(i_file) = time
     particle_masses = aero_state_masses(aero_state, aero_data)
     call ensure_real_array_size(particle_volumes, aero_state%apa%n_part)
     particle_volumes = aero_particle_volume( &
          aero_state%apa%particle(1:aero_state%apa%n_part))
     density = sum(particle_masses) / sum(particle_volumes)
     if (regime == REGIME_FREE) then
        dimless_times(i_file) = (6d0 * const%boltzmann * env_state%temp &
             * aero_data%fractal%prime_radius / density)**(1d0 / 2d0)   &
             * n_init * times(i_file)
     elseif (regime == REGIME_CONT) then
        dimless_times(i_file) = (2d0 * const%boltzmann * env_state%temp &
             / 3d0 / const%air_dyn_visc) * n_init * times(i_file)
     else
        call die(612347164)
     end if
     particle_num_concs = aero_state_num_concs(aero_state, aero_data)
     time_num_concs(i_file) = sum(particle_num_concs)
     dimless_time_num_concs(i_file) = time_num_concs(i_file) / n_init
  end do

  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column  1: dimensionless time"
  write(*,'(a)') "    column  2: dimensionless aerosol number concentration"

  call open_file_write(out_filename, out_unit)
  do i_file = 1,n_file
     write(out_unit,'(e30.15e3)',advance='no') dimless_times(i_file)
     write(out_unit,'(e30.15e3)',advance='no') dimless_time_num_concs(i_file)
     write(out_unit,'(a)') ''
  end do
  call close_file(out_unit)

  deallocate(times)
  deallocate(dimless_times)
  deallocate(time_num_concs)
  deallocate(dimless_time_num_concs)
  deallocate(particle_num_concs)
  deallocate(particle_masses)
  deallocate(particle_volumes)
  deallocate(filename_list)

  call pmc_mpi_finalize()

end program test_fractal_dimless_time
