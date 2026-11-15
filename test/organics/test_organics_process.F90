! Copyright (C) 2021 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The organics test process program.

!> Read a NetCDF output file and write the O:C, OM:OC, organic mass
!> and organic carbon mass of the particles of the three monodisperse
!> populations.
program test_organics_process

  use pmc_output

  !> Number of populations.
  integer, parameter :: N_POP = 3
  !> Upper diameter limits of the populations (m).
  real(kind=dp), parameter :: POP_DIAM_MAX(N_POP) = [1.5d-7, 3d-7, 6d-7]

  character(len=PMC_MAX_FILENAME_LEN) :: prefix, in_filename, out_filename
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: index, repeat, unit, i_pop, i_part, n_pop_part
  real(kind=dp) :: time, del_t, diam_min
  character(len=PMC_UUID_LEN) :: uuid
  real(kind=dp), allocatable :: diameters(:), o_cs(:), om_ocs(:), &
       org_masses(:), org_carbon_masses(:)
  logical, allocatable :: in_pop(:)

  call pmc_mpi_init()

  call assert_msg(350372166, command_argument_count() == 1, &
       "usage: test_organics_process <netcdf_prefix>")
  call get_command_argument(1, prefix)

  call make_filename(in_filename, prefix, ".nc", 1, 1)
  call input_state(in_filename, index, time, del_t, repeat, &
       uuid, aero_data=aero_data, aero_state=aero_state, &
       env_state=env_state)

  diameters = aero_state_diameters(aero_state, aero_data)
  o_cs = aero_state_o_c(aero_state, aero_data)
  om_ocs = aero_state_om_oc(aero_state, aero_data)
  org_masses = aero_state_organic_masses(aero_state, aero_data)
  org_carbon_masses = aero_state_organic_carbon_masses(aero_state, &
       aero_data)

  ! one line per population with the smallest and the largest value of
  ! each quantity (all particles of a population are identical)
  out_filename = trim(prefix) // "_o_c.txt"
  call open_file_write(out_filename, unit)
  allocate(in_pop(aero_state_n_part(aero_state)))
  diam_min = 0d0
  do i_pop = 1,N_POP
     in_pop = (diameters > diam_min) .and. (diameters <= POP_DIAM_MAX(i_pop))
     n_pop_part = count(in_pop)
     call assert_msg(703158942, n_pop_part > 0, "empty population")
     write(*,'(a,i1,a,i6,a,g14.6,a,g14.6)') 'population ', i_pop, ': ', &
          n_pop_part, ' particles, O:C = ', maxval(o_cs, mask=in_pop), &
          ', OM:OC = ', maxval(om_ocs, mask=in_pop)
     write(unit,'(8e30.20)') &
          minval(o_cs, mask=in_pop), maxval(o_cs, mask=in_pop), &
          minval(om_ocs, mask=in_pop), maxval(om_ocs, mask=in_pop), &
          minval(org_masses, mask=in_pop), maxval(org_masses, mask=in_pop), &
          minval(org_carbon_masses, mask=in_pop), &
          maxval(org_carbon_masses, mask=in_pop)
     diam_min = POP_DIAM_MAX(i_pop)
  end do
  call close_file(unit)

  call pmc_mpi_finalize()

end program test_organics_process
