! Copyright (C) 2009-2013, 2016, 2017, 2021 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program for the 8_organics scenario.

!> Read NetCDF output files and compute the organic aerosol diagnostics:
!> bulk mass concentrations of groups of organic species, bulk O:C and
!> OM:OC, and per-particle distributions of O:C and of the secondary
!> organic mass fraction.
program process

  use pmc_output
  use pmc_stats

  character(len=PMC_MAX_FILENAME_LEN), parameter :: prefix &
       = "out/organics"

  !> Number of groups of organic species.
  integer, parameter :: N_GROUP = 6
  !> Maximum number of species per group.
  integer, parameter :: N_GROUP_SPEC = 14
  !> Names of the groups (used in the names of the output variables).
  character(len=20), parameter :: GROUP_NAME(N_GROUP) = [ &
       character(len=20) :: "poa", "soa_roc", "soa_monoterpene", &
       "soa_isoprene", "soa_glyoxal", "soa_cloud"]
  !> Whether the groups are secondary organic matter.
  logical, parameter :: GROUP_SECONDARY(N_GROUP) = [.false., .true., &
       .true., .true., .true., .true.]

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, o_c_grid, frac_grid
  type(aero_data_t) :: aero_data
  type(aero_state_t) :: aero_state
  type(env_state_t) :: env_state
  integer :: ncid, index, repeat, i_index, i_repeat, n_index, n_repeat
  integer :: i_group, i_name, i_spec, i_part
  real(kind=dp) :: time, del_t, tot_num_conc, tot_dry_mass_conc
  real(kind=dp) :: carbon_conc, oxygen_conc, sec_carbon_conc, &
       sec_oxygen_conc, sec_mass_conc, org_mass_conc
  character(len=PMC_UUID_LEN) :: uuid
  character(len=AERO_NAME_LEN) :: group_spec(N_GROUP_SPEC, N_GROUP)
  logical, allocatable :: is_secondary(:)
  real(kind=dp), allocatable :: times(:), dry_diameters(:), num_concs(:), &
       dry_masses(:), sec_masses(:), sec_fracs(:), org_carbon_masses(:), &
       o_cs(:), spec_mass_concs(:), group_mass_concs(:), num_dist(:), &
       diam_o_c_dist(:,:), diam_sec_frac_dist(:,:)
  type(stats_1d_t) :: stats_num_dist, stats_tot_num_conc, &
       stats_tot_dry_mass_conc, stats_tot_org_mass_conc, &
       stats_tot_sec_org_mass_conc, stats_tot_o_c, stats_tot_om_oc, &
       stats_sec_o_c, stats_group_mass_conc(N_GROUP)
  type(stats_2d_t) :: stats_diam_o_c_dist, stats_diam_sec_frac_dist

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(o_c_grid, BIN_GRID_TYPE_LINEAR, 70, 0d0, 1.4d0)
  call bin_grid_make(frac_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)

  allocate(times(n_index))

  ! groups of the organic aerosol species of CRACMM3, as in the species
  ! definition file of CMAQ (SpecDef_Conc_cracmm3.txt)
  group_spec = ""
  ! primary organic matter: the alkane-like reactive organic carbon
  ! (APOM of CMAQ) and the nonvolatile OC of the initial and background
  ! particles
  group_spec(1:7,1) = [character(len=AERO_NAME_LEN) :: &
       "AROCN2ALK", "AROCN1ALK", "AROCP0ALK", "AROCP1ALK", "AROCP2ALK", &
       "AROCP3ALK", "OC"]
  ! oxygenated reactive organic carbon, aromatic peroxides and the SOA of
  ! phenols and other species (AORGA of CMAQ): SOA of aromatics, alkanes,
  ! sesquiterpenes and organic gases of intermediate volatility; the
  ! oxygenated primary organic matter of cooking is emitted as these
  ! species and counts as secondary here, as in CMAQ
  group_spec(1:14,2) = [character(len=AERO_NAME_LEN) :: &
       "AROCN2OXY2", "AROCN2OXY4", "AROCN2OXY8", "AROCN1OXY1", &
       "AROCN1OXY3", "AROCN1OXY6", "AROCP0OXY2", "AROCP0OXY4", &
       "AROCP1OXY1", "AROCP1OXY3", "AROCP2OXY2", "AROCP3OXY2", "AOP3", &
       "ASOAT"]
  ! SOA from monoterpenes: highly oxygenated molecules, organic nitrates
  ! and the photolysis product of the HOM (part of AORGB of CMAQ)
  group_spec(1:5,3) = [character(len=AERO_NAME_LEN) :: &
       "AHOM", "AELHOM", "ATRPN", "AHONIT", "AMTN1"]
  ! SOA from isoprene: uptake of the epoxydiols, of the tetrafunctional
  ! compounds and of the nitrates (part of AORGB of CMAQ)
  group_spec(1:4,4) = [character(len=AERO_NAME_LEN) :: &
       "AISO3NOS", "AISO3OS", "AISO4", "AISO5"]
  ! SOA from the uptake of methylglyoxal
  group_spec(1:1,5) = [character(len=AERO_NAME_LEN) :: "AGLYOLIG"]
  ! SOA from cloud processing (not formed in PartMC-MOSAIC)
  group_spec(1:1,6) = [character(len=AERO_NAME_LEN) :: "AORGC"]

  o_cs = [ real(kind=dp) :: ] ! silence compiler warnings
  sec_fracs = [ real(kind=dp) :: ]

  do i_index = 1,n_index
     do i_repeat = 1,n_repeat
        call make_filename(in_filename, prefix, ".nc", i_index, i_repeat)
        write(*,*) "Processing " // trim(in_filename)
        call input_state(in_filename, index, time, del_t, repeat, &
             uuid, aero_data=aero_data, aero_state=aero_state, &
             env_state=env_state)

        times(i_index) = time

        dry_diameters = aero_state_dry_diameters(aero_state, aero_data)
        num_concs = aero_state_num_concs(aero_state, aero_data)
        num_dist = bin_grid_histogram_1d(diam_grid, dry_diameters, num_concs)
        call stats_1d_add(stats_num_dist, num_dist)

        tot_num_conc = sum(num_concs)
        call stats_1d_add_entry(stats_tot_num_conc, tot_num_conc, i_index)

        dry_masses = aero_state_masses(aero_state, aero_data, &
             exclude=(/"H2O"/))
        tot_dry_mass_conc = sum(dry_masses * num_concs)
        call stats_1d_add_entry(stats_tot_dry_mass_conc, tot_dry_mass_conc, &
             i_index)

        ! bulk mass concentration of every species
        if (allocated(spec_mass_concs)) deallocate(spec_mass_concs)
        allocate(spec_mass_concs(aero_data_n_spec(aero_data)))
        spec_mass_concs = 0d0
        do i_part = 1,aero_state_n_part(aero_state)
           spec_mass_concs = spec_mass_concs + num_concs(i_part) &
                * aero_particle_species_masses( &
                aero_state%apa%particle(i_part), aero_data)
        end do

        ! bulk mass concentrations of the groups, and the species that
        ! are secondary organic matter
        if (allocated(is_secondary)) deallocate(is_secondary)
        allocate(is_secondary(aero_data_n_spec(aero_data)))
        is_secondary = .false.
        if (allocated(group_mass_concs)) deallocate(group_mass_concs)
        allocate(group_mass_concs(N_GROUP))
        group_mass_concs = 0d0
        do i_group = 1,N_GROUP
           do i_name = 1,N_GROUP_SPEC
              if (len_trim(group_spec(i_name, i_group)) == 0) cycle
              i_spec = aero_data_spec_by_name(aero_data, &
                   group_spec(i_name, i_group))
              if (i_spec == 0) cycle
              group_mass_concs(i_group) = group_mass_concs(i_group) &
                   + spec_mass_concs(i_spec)
              is_secondary(i_spec) = GROUP_SECONDARY(i_group)
           end do
           call stats_1d_add_entry(stats_group_mass_conc(i_group), &
                group_mass_concs(i_group), i_index)
        end do

        ! bulk O:C and OM:OC of all organic matter and bulk O:C of the
        ! secondary organic matter, all weighted by carbon; oxygen_conc
        ! is the carbon concentration times O:C
        org_mass_conc = sum(spec_mass_concs, mask=(aero_data%om_oc > 0d0))
        carbon_conc = sum(spec_mass_concs / aero_data%om_oc, &
             mask=(aero_data%om_oc > 0d0))
        oxygen_conc = sum(spec_mass_concs / aero_data%om_oc &
             * aero_data%o_c, mask=(aero_data%om_oc > 0d0))
        sec_mass_conc = sum(spec_mass_concs, mask=is_secondary)
        sec_carbon_conc = sum(spec_mass_concs / aero_data%om_oc, &
             mask=(is_secondary .and. (aero_data%om_oc > 0d0)))
        sec_oxygen_conc = sum(spec_mass_concs / aero_data%om_oc &
             * aero_data%o_c, &
             mask=(is_secondary .and. (aero_data%om_oc > 0d0)))
        call stats_1d_add_entry(stats_tot_org_mass_conc, org_mass_conc, &
             i_index)
        call stats_1d_add_entry(stats_tot_sec_org_mass_conc, sec_mass_conc, &
             i_index)
        if (carbon_conc > 0d0) then
           call stats_1d_add_entry(stats_tot_o_c, &
                oxygen_conc / carbon_conc, i_index)
           call stats_1d_add_entry(stats_tot_om_oc, &
                org_mass_conc / carbon_conc, i_index)
        else
           call stats_1d_add_entry(stats_tot_o_c, 0d0, i_index)
           call stats_1d_add_entry(stats_tot_om_oc, 0d0, i_index)
        end if
        if (sec_carbon_conc > 0d0) then
           call stats_1d_add_entry(stats_sec_o_c, &
                sec_oxygen_conc / sec_carbon_conc, i_index)
        else
           call stats_1d_add_entry(stats_sec_o_c, 0d0, i_index)
        end if

        ! distribution of the per-particle O:C, over the particles that
        ! contain organic carbon
        org_carbon_masses = aero_state_organic_carbon_masses(aero_state, &
             aero_data)
        o_cs = aero_state_o_c(aero_state, aero_data)
        diam_o_c_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             o_c_grid, o_cs, merge(num_concs, 0d0, org_carbon_masses > 0d0))
        call stats_2d_add(stats_diam_o_c_dist, diam_o_c_dist)

        ! distribution of the secondary organic fraction of the dry mass
        if (allocated(sec_masses)) deallocate(sec_masses)
        allocate(sec_masses(aero_state_n_part(aero_state)))
        do i_part = 1,aero_state_n_part(aero_state)
           sec_masses(i_part) = sum(aero_particle_species_masses( &
                aero_state%apa%particle(i_part), aero_data), &
                mask=is_secondary)
        end do
        sec_fracs = sec_masses / dry_masses
        diam_sec_frac_dist = bin_grid_histogram_2d(diam_grid, &
             dry_diameters, frac_grid, sec_fracs, num_concs)
        call stats_2d_add(stats_diam_sec_frac_dist, diam_sec_frac_dist)

     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "8_organics process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(o_c_grid, ncid, "o_c", unit="1")
     call bin_grid_output_netcdf(frac_grid, ncid, "sec_frac", unit="1")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)

     call stats_2d_output_netcdf(stats_diam_o_c_dist, ncid, "diam_o_c_dist", &
          dim_name_1="diam", dim_name_2="o_c", unit="m^{-3}")
     call stats_2d_clear(stats_diam_o_c_dist)

     call stats_2d_output_netcdf(stats_diam_sec_frac_dist, ncid, &
          "diam_sec_frac_dist", dim_name_1="diam", dim_name_2="sec_frac", &
          unit="m^{-3}")
     call stats_2d_clear(stats_diam_sec_frac_dist)

     call pmc_nc_close(ncid)
  end do

  call make_filename(out_filename, prefix, "_process.nc")
  write(*,*) "Writing " // trim(out_filename)
  call pmc_nc_open_write(out_filename, ncid)
  call pmc_nc_write_info(ncid, uuid, "8_organics process")
  call pmc_nc_write_real_1d(ncid, times, "time", dim_name="time", unit="s")
  call stats_1d_output_netcdf(stats_tot_num_conc, ncid, "tot_num_conc", &
       dim_name="time", unit="m^{-3}")
  call stats_1d_output_netcdf(stats_tot_dry_mass_conc, ncid, &
       "tot_dry_mass_conc", dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_org_mass_conc, ncid, &
       "tot_org_mass_conc", dim_name="time", unit="kg m^{-3}")
  call stats_1d_output_netcdf(stats_tot_sec_org_mass_conc, ncid, &
       "tot_sec_org_mass_conc", dim_name="time", unit="kg m^{-3}")
  do i_group = 1,N_GROUP
     call stats_1d_output_netcdf(stats_group_mass_conc(i_group), ncid, &
          "mass_conc_" // trim(GROUP_NAME(i_group)), dim_name="time", &
          unit="kg m^{-3}")
  end do
  call stats_1d_output_netcdf(stats_tot_o_c, ncid, "tot_o_c", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_tot_om_oc, ncid, "tot_om_oc", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_sec_o_c, ncid, "sec_o_c", &
       dim_name="time", unit="1")
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

end program process
