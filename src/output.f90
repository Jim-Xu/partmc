! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_output module.

!> Write data in NetCDF format.
module pmc_output

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_netcdf
  use pmc_env_state
  use pmc_util
  use pmc_gas_data
  use pmc_mpi
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current state.
  subroutine output_state_netcdf(prefix, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop, record_removals)

    !> Prefix of state file.
    character(len=*), intent(in) :: prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current timestep (s).
    real*8, intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop
    !> Whether to output particle removal info.
    logical, intent(in) :: record_removals
    
    character(len=len(prefix)+100) :: filename
    character(len=500) :: history
    type(spec_file_t) :: file
    type(env_state_t) :: env_write
    type(gas_state_t) :: gas_state_write
    type(aero_state_t) :: aero_state_write
    integer :: ierr, status, buffer_size, i_proc, position
    character, allocatable :: buffer(:)
    integer :: ncid

    ! only root node actually writes to the file
    if (pmc_mpi_rank() == 0) then
       write(filename, '(a,a,i4.4,a,i8.8,a)') trim(prefix), &
            '_', i_loop, '_', index, '.nc'
       call pmc_nc_check_msg(nf90_create(filename, NF90_CLOBBER, ncid), &
            "opening " // trim(filename))
       
       call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
            "PartMC output file"))
       call iso8601_date_and_time(history)
       history((len_trim(history)+1):) = " created by PartMC"
       call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))
       
       call pmc_nc_check(nf90_enddef(ncid))

       call pmc_nc_write_real(ncid, time, "time", "s")
       call pmc_nc_write_real(ncid, del_t, "timestep", "s")
       call pmc_nc_write_integer(ncid, i_loop, "loop", "1")
       call pmc_nc_write_integer(ncid, index, "timestep_index", "1")

       call env_state_output_netcdf(env_state, ncid)
       call gas_data_output_netcdf(gas_data, ncid)
       call gas_state_output_netcdf(gas_state, ncid, gas_data)
       call aero_data_output_netcdf(aero_data, ncid)
       call aero_state_output_netcdf(aero_state, ncid, bin_grid, &
            aero_data, record_removals)
    end if

#ifdef PMC_USE_MPI
    ! FIXME: old stuff from text output

    ! write everyone else's state
    do i_proc = 1,(pmc_mpi_size() - 1)
       call pmc_mpi_barrier()
       ! compute and send buffer_size from remote node
       if (pmc_mpi_rank() == i_proc) then
          buffer_size = 0
          buffer_size = buffer_size + pmc_mpi_pack_size_env_state(env_state)
          buffer_size = buffer_size + pmc_mpi_pack_size_gas_state(gas_state)
          buffer_size = buffer_size + pmc_mpi_pack_size_aero_state(aero_state)
          call mpi_send(buffer_size, 1, MPI_INTEGER, 0, 57, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
       end if
       ! get buffer_size at root node
       if (pmc_mpi_rank() == 0) then
          call mpi_recv(buffer_size, 1, MPI_INTEGER, i_proc, 57, &
               MPI_COMM_WORLD, status, ierr)
          call pmc_mpi_check_ierr(ierr)
       end if
       ! send buffer from remote node
       if (pmc_mpi_rank() == i_proc) then
          allocate(buffer(buffer_size))
          position = 0
          call pmc_mpi_pack_env_state(buffer, position, env_state)
          call pmc_mpi_pack_gas_state(buffer, position, gas_state)
          call pmc_mpi_pack_aero_state(buffer, position, aero_state)
          call assert(542772170, position == buffer_size)
          call mpi_send(buffer, buffer_size, MPI_CHARACTER, 0, 58, &
               MPI_COMM_WORLD, ierr)
          call pmc_mpi_check_ierr(ierr)
          deallocate(buffer)
       end if
       ! get buffer at root node
       if (pmc_mpi_rank() == 0) then
          allocate(buffer(buffer_size))
          call mpi_recv(buffer, buffer_size, MPI_CHARACTER, i_proc, 58, &
               MPI_COMM_WORLD, status, ierr)
          position = 0
          call pmc_mpi_unpack_env_state(buffer, position, env_write)
          call pmc_mpi_unpack_gas_state(buffer, position, gas_state_write)
          call pmc_mpi_unpack_aero_state(buffer, position, aero_state_write)
          call assert(518174881, position == buffer_size)
          deallocate(buffer)
       end if
       ! process state at root node
       if (pmc_mpi_rank() == 0) then
          !call write_integer(file, 'processor', i_proc)
          !call write_env_state(file, env_write)
          !call write_gas_state(file, gas_state_write)
          !call write_aero_state(file, aero_state_write)
          
          call env_state_deallocate(env_write)
          call gas_state_deallocate(gas_state_write)
          call aero_state_deallocate(aero_state_write)
       end if
    end do
#endif
       
    if (pmc_mpi_rank() == 0) then
       call pmc_nc_check(nf90_close(ncid))
    end if
    
  end subroutine output_state_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the current state.
  subroutine input_state_netcdf(filename, bin_grid, aero_data, &
       aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop)

    !> Prefix of state file.
    character(len=*), intent(in) :: filename
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(out) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(out) :: gas_state
    !> Environment state.
    type(env_state_t), intent(out) :: env_state
    !> Filename index.
    integer, intent(out) :: index
    !> Current time (s).
    real*8, intent(out) :: time
    !> Current timestep (s).
    real*8, intent(out) :: del_t
    !> Current loop number.
    integer, intent(out) :: i_loop
    
    type(spec_file_t) :: file
    type(env_state_t) :: env_write
    type(gas_state_t) :: gas_state_write
    type(aero_state_t) :: aero_state_write
    integer :: ierr, status, buffer_size, i_proc, position
    character, allocatable :: buffer(:)
    integer :: ncid
    character(len=1000) :: unit

    ! only root node actually reads from the file
    if (pmc_mpi_rank() == 0) then
       call pmc_nc_open_read(filename, ncid)

       call pmc_nc_read_real(ncid, time, "time", unit)
       call pmc_nc_read_real(ncid, del_t, "timestep", unit)
       call pmc_nc_read_integer(ncid, i_loop, "loop", unit)
       call pmc_nc_read_integer(ncid, index, "timestep_index", unit)

       call env_state_input_netcdf(env_state, ncid)
       call gas_data_input_netcdf(gas_data, ncid)
       call gas_state_input_netcdf(gas_state, ncid, gas_data)
       call aero_data_input_netcdf(aero_data, ncid)
       call aero_state_input_netcdf(aero_state, ncid, bin_grid, &
            aero_data)

       call pmc_nc_close(ncid)
    end if
    
  end subroutine input_state_netcdf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current sectional data.
  subroutine output_sectional(prefix, bin_grid, aero_data, &
       aero_binned, gas_data, gas_state, env_state, index, time, del_t)
    
    !> Prefix of filename to write
    character(len=*), intent(in) :: prefix
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Binned aerosol data.
    type(aero_binned_t), intent(in) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current output time-step (s).
    real*8, intent(in) :: del_t

    integer :: ncid
    character(len=len(prefix)+100) :: filename
    character(len=500) :: history

    write(filename, '(a,a,i8.8,a)') trim(prefix), &
         '_', index, '.nc'
    call pmc_nc_check_msg(nf90_create(filename, NF90_CLOBBER, ncid), &
         "opening " // trim(filename))

    ! write header attributes
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC sectional output file"))
    call iso8601_date_and_time(history)
    history((len_trim(history)+1):) = " created by PartMC"
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))

    call pmc_nc_check(nf90_enddef(ncid))

    ! write data
    call pmc_nc_write_real(ncid, time, "time", "s")
    call pmc_nc_write_real(ncid, del_t, "timestep", "s")
    call pmc_nc_write_integer(ncid, index, "timestep_index", "1")
    
    call env_state_output_netcdf(env_state, ncid)
    call gas_data_output_netcdf(gas_data, ncid)
    call gas_state_output_netcdf(gas_state, ncid, gas_data)
    call aero_data_output_netcdf(aero_data, ncid)
    call aero_binned_output_netcdf(aero_binned, ncid, bin_grid, &
         aero_data)

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_sectional
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_output
