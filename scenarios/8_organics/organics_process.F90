! Copyright (C) 2009-2013, 2016, 2017, 2021 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The process program for the 8_organics scenario.

!> Read NetCDF output files and compute the organic aerosol diagnostics:
!> bulk mass concentrations of groups of organic species, bulk O:C and
!> OM:OC, and per-particle distributions of O:C and of the secondary
!> organic mass fraction; the liquid-liquid phase separation of the
!> particles that contain organic matter and inorganic salt, from the
!> separation relative humidity of their OM:OC; and the critical
!> supersaturation of every particle (kappa-Koehler theory with the
!> surface tension of water and with the effective surface tension of
!> the organic film) with the CCN concentrations at fixed
!> supersaturations; between these two limits, the CCN concentrations
!> with the film only on the particles whose organic matter has an O:C
!> up to a limit, for three such limits from laboratory studies of the
!> phase separation near saturation. The phase separation and the CCN
!> concentrations are also computed for two internal mixtures of the
!> same particles, to show what the per-particle composition changes:
!> every particle with the average dry composition of all particles, and
!> every particle with the average dry composition of the particles in
!> its size section.
!> The number and mass concentrations, the O:C, the phase separation and
!> the CCN concentrations are also computed for the particles of each
!> class by origin: the urban primary particles, the urban and regional
!> particles joined by coagulation, the regional particles of the
!> initial air, and the regional particles that the exchange with the
!> background air brought in.
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
  !> Number of supersaturations of the CCN concentrations.
  integer, parameter :: N_SS = 5
  !> Supersaturations of the CCN concentrations (%).
  real(kind=dp), parameter :: CCN_SS(N_SS) = [0.1d0, 0.2d0, 0.3d0, &
       0.5d0, 1d0]
  !> Names of the supersaturations in the output variables.
  character(len=4), parameter :: CCN_SS_NAME(N_SS) = [ &
       character(len=4) :: "s010", "s020", "s030", "s050", "s100"]
  !> Number of limits of the O:C for the organic film at activation.
  integer, parameter :: N_FILM = 3
  !> Limits of the O:C of the organic matter of a particle (1): up to
  !> the limit, the organic matter is taken to stay in a phase of its
  !> own until the particle activates, so that the particle keeps its
  !> organic film (see llps_ccn below). These are the O:C up to which
  !> two liquid phases were observed at relative humidities of 95 to
  !> 100% or at supersaturation: 0.44 in secondary organic material
  !> without salts and in particles of one organic compound
  !> (Renbaum-Wolff et al., 2016; Song et al., 2017, 2018, 2020; Ham et
  !> al., 2019), 0.67 as the average of particles of two organic
  !> compounds (Song et al., 2020), and 0.75 as the average of the
  !> organic matter mixed with ammonium sulfate (Malek et al., 2023).
  real(kind=dp), parameter :: FILM_O_C(N_FILM) = [0.44d0, 0.67d0, 0.75d0]
  !> Names of the limits in the output variables.
  character(len=5), parameter :: FILM_NAME(N_FILM) = [ &
       character(len=5) :: "oc044", "oc067", "oc075"]
  !> Inorganic species that form the salt for the phase separation
  !> (those missing from the aerosol data are skipped).
  character(len=AERO_NAME_LEN), parameter :: SALT_SPEC(5) = [ &
       character(len=AERO_NAME_LEN) :: "SO4", "NO3", "NH4", "Cl", "Na"]
  !> Number of internal mixtures.
  integer, parameter :: N_MIX = 2
  !> Names of the internal mixtures in the output variables: all
  !> particles mixed, and the particles of each size section mixed.
  character(len=7), parameter :: MIX_NAME(N_MIX) = [ &
       character(len=7) :: "intmix", "sectmix"]
  !> Whether the internal mixtures are mixed by size section.
  logical, parameter :: MIX_SECTIONAL(N_MIX) = [.false., .true.]
  !> Number of size sections.
  integer, parameter :: N_SECT = 8
  !> Lower edge of the first and upper edge of the last size section (dry
  !> diameter, m); the sections are equally spaced in the logarithm of
  !> the diameter, the upper edge of each twice its lower edge, and the
  !> first and the last take the particles below and above. These are
  !> the eight sections of MOSAIC in WRF-Chem (version 3.9.1,
  !> chem/module_mosaic_driver.F, init_data_mosaic_asect).
  real(kind=dp), parameter :: SECT_DIAM_MIN = 39.0625d-9
  real(kind=dp), parameter :: SECT_DIAM_MAX = 10d-6
  !> Number of classes of particles by origin.
  integer, parameter :: N_ORIGIN = 4
  !> Names of the classes by origin (used in the names of the output
  !> variables; see particle_origins below).
  character(len=13), parameter :: ORIGIN_NAME(N_ORIGIN) = [ &
       character(len=13) :: "city", "city_regional", "regional", &
       "background"]
  !> Source of the regional particles: the name of their mode in the
  !> initial and in the background distribution (aero_init_dist.dat and
  !> aero_back_dist.dat, written by make_case.py). The other sources are
  !> the urban modes of the initial air.
  character(len=AERO_SOURCE_NAME_LEN), parameter :: REGIONAL_SOURCE &
       = "regional"

  character(len=PMC_MAX_FILENAME_LEN) :: in_filename, out_filename
  type(bin_grid_t) :: diam_grid, o_c_grid, frac_grid, sc_grid, kappa_grid, &
       srh_grid, entry_grid
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
       diam_o_c_dist(:,:), diam_sec_frac_dist(:,:), org_masses(:), &
       srhs(:), scs(:), scs_est(:), kappas(:), &
       diam_srh_dist(:,:), diam_sc_dist(:,:), diam_sc_est_dist(:,:), &
       diam_kappa_dist(:,:), entry_times(:), entry_srh_dist(:,:), &
       scs_film(:,:)
  logical, allocatable :: is_mixed(:), is_separated(:), has_film(:,:)
  real(kind=dp) :: tot_org_part_mass_conc
  integer :: i_ss, i_film
  type(aero_state_t) :: mixed_aero_state
  integer :: i_mix
  real(kind=dp) :: srh_intmix, tot_org_mix_mass_conc
  real(kind=dp), allocatable :: mix_org_masses(:), mix_srhs(:), &
       mix_kappas(:), mix_scs(:), mix_scs_est(:), mix_scs_film(:,:)
  logical, allocatable :: mix_is_mixed(:), mix_is_separated(:), &
       mix_has_film(:,:)
  integer :: i_origin
  integer, allocatable :: origins(:)
  logical, allocatable :: in_origin(:), is_dry(:), is_salt(:)
  real(kind=dp) :: origin_org_mass_conc, origin_carbon_conc, &
       origin_oxygen_conc
  real(kind=dp), allocatable :: origin_spec_mass_concs(:,:)
  type(stats_1d_t) :: stats_origin_num_dist(N_ORIGIN), &
       stats_origin_num_conc(N_ORIGIN), &
       stats_origin_mixed_num_conc(N_ORIGIN), &
       stats_origin_llps_num_conc(N_ORIGIN), &
       stats_origin_llps_org_mass_conc(N_ORIGIN), &
       stats_origin_ccn_conc(N_SS, N_ORIGIN), &
       stats_origin_ccn_est_conc(N_SS, N_ORIGIN), &
       stats_origin_ccn_film_conc(N_SS, N_FILM, N_ORIGIN), &
       stats_origin_dry_mass_conc(N_ORIGIN), &
       stats_origin_salt_mass_conc(N_ORIGIN), &
       stats_origin_org_mass_conc(N_ORIGIN), &
       stats_origin_group_mass_conc(N_GROUP, N_ORIGIN), &
       stats_origin_o_c(N_ORIGIN), stats_origin_om_oc(N_ORIGIN)
  type(stats_1d_t) :: stats_num_dist, stats_tot_num_conc, &
       stats_tot_dry_mass_conc, stats_tot_org_mass_conc, &
       stats_tot_sec_org_mass_conc, stats_tot_o_c, stats_tot_om_oc, &
       stats_sec_o_c, stats_group_mass_conc(N_GROUP), stats_rel_humid, &
       stats_mixed_num_frac, stats_llps_num_frac, stats_llps_org_frac, &
       stats_ccn_conc(N_SS), stats_ccn_est_conc(N_SS), stats_srh_intmix, &
       stats_mix_llps_num_frac(N_MIX), stats_mix_llps_org_frac(N_MIX), &
       stats_mix_ccn_conc(N_SS, N_MIX), stats_mix_ccn_est_conc(N_SS, N_MIX), &
       stats_film_num_frac(N_FILM), stats_ccn_film_conc(N_SS, N_FILM), &
       stats_mix_ccn_film_conc(N_SS, N_FILM, N_MIX)
  type(stats_2d_t) :: stats_diam_o_c_dist, stats_diam_sec_frac_dist, &
       stats_diam_srh_dist, stats_diam_sc_dist, stats_diam_sc_est_dist, &
       stats_diam_kappa_dist, stats_entry_srh_dist

  call pmc_mpi_init()

  call input_n_files(prefix, n_repeat, n_index)

  call bin_grid_make(diam_grid, BIN_GRID_TYPE_LOG, 180, 1d-9, 1d-3)
  call bin_grid_make(o_c_grid, BIN_GRID_TYPE_LINEAR, 70, 0d0, 1.4d0)
  call bin_grid_make(frac_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  call bin_grid_make(sc_grid, BIN_GRID_TYPE_LOG, 50, 1d-4, 1d0)
  call bin_grid_make(kappa_grid, BIN_GRID_TYPE_LINEAR, 80, 0d0, 0.8d0)
  call bin_grid_make(srh_grid, BIN_GRID_TYPE_LINEAR, 50, 0d0, 1d0)
  ! entry time of the particles in the box (s): one bin per hour of the
  ! two days of the case, centered on the whole hours, so that the
  ! particles of the initial air (entry time 0) are in the first bin
  call bin_grid_make(entry_grid, BIN_GRID_TYPE_LINEAR, 49, -1800d0, &
       174600d0)

  allocate(times(n_index))

  ! groups of the organic aerosol species of CRACMM3, as in the species
  ! definition file of CMAQ (SpecDef_Conc_cracmm3.txt)
  group_spec = ""
  ! primary organic matter: the alkane-like reactive organic carbon
  ! (APOM of CMAQ) and the nonvolatile OC (which the present case of the
  ! scenario does not use)
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
  ! SOA from cloud processing, which PartMC-MOSAIC does not form; in the
  ! present case of the scenario it is the aged regional organic matter
  ! of the initial and background particles
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

        ! liquid-liquid phase separation and critical supersaturation of
        ! every particle (see llps_ccn below)
        call llps_ccn(aero_state, aero_data, env_state, org_masses, &
             is_mixed, srhs, is_separated, kappas, scs, scs_est, has_film, &
             scs_film)
        call stats_1d_add_entry(stats_rel_humid, env_state%rel_humid, &
             i_index)
        if (tot_num_conc > 0d0) then
           call stats_1d_add_entry(stats_mixed_num_frac, &
                sum(num_concs, mask=is_mixed) / tot_num_conc, i_index)
           call stats_1d_add_entry(stats_llps_num_frac, &
                sum(num_concs, mask=is_separated) / tot_num_conc, i_index)
        else
           call stats_1d_add_entry(stats_mixed_num_frac, 0d0, i_index)
           call stats_1d_add_entry(stats_llps_num_frac, 0d0, i_index)
        end if
        tot_org_part_mass_conc = sum(org_masses * num_concs)
        if (tot_org_part_mass_conc > 0d0) then
           call stats_1d_add_entry(stats_llps_org_frac, &
                sum(org_masses * num_concs, mask=is_separated) &
                / tot_org_part_mass_conc, i_index)
        else
           call stats_1d_add_entry(stats_llps_org_frac, 0d0, i_index)
        end if
        diam_srh_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             srh_grid, srhs, merge(num_concs, 0d0, is_mixed))
        call stats_2d_add(stats_diam_srh_dist, diam_srh_dist)
        ! the same against the time at which the oldest constituent of the
        ! particle entered the box: 0 for the particles of the initial air,
        ! the time of the exchange for those of the background air (see
        ! particle_origins below)
        if (allocated(entry_times)) deallocate(entry_times)
        allocate(entry_times(aero_state_n_part(aero_state)))
        do i_part = 1,aero_state_n_part(aero_state)
           entry_times(i_part) &
                = aero_state%apa%particle(i_part)%least_create_time
        end do
        entry_srh_dist = bin_grid_histogram_2d(entry_grid, entry_times, &
             srh_grid, srhs, merge(num_concs, 0d0, is_mixed))
        call stats_2d_add(stats_entry_srh_dist, entry_srh_dist)

        do i_ss = 1,N_SS
           call stats_1d_add_entry(stats_ccn_conc(i_ss), &
                sum(num_concs, mask=(scs <= CCN_SS(i_ss) / 100d0)), i_index)
           call stats_1d_add_entry(stats_ccn_est_conc(i_ss), &
                sum(num_concs, mask=(scs_est <= CCN_SS(i_ss) / 100d0)), &
                i_index)
        end do
        ! with the organic film only on the particles whose organic matter
        ! has an O:C up to each limit, and the number fraction of these
        ! particles
        do i_film = 1,N_FILM
           if (tot_num_conc > 0d0) then
              call stats_1d_add_entry(stats_film_num_frac(i_film), &
                   sum(num_concs, mask=has_film(:, i_film)) / tot_num_conc, &
                   i_index)
           else
              call stats_1d_add_entry(stats_film_num_frac(i_film), 0d0, &
                   i_index)
           end if
           do i_ss = 1,N_SS
              call stats_1d_add_entry(stats_ccn_film_conc(i_ss, i_film), &
                   sum(num_concs, mask=(scs_film(:, i_film) &
                   <= CCN_SS(i_ss) / 100d0)), i_index)
           end do
        end do
        diam_sc_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs, num_concs)
        call stats_2d_add(stats_diam_sc_dist, diam_sc_dist)
        diam_sc_est_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             sc_grid, scs_est, num_concs)
        call stats_2d_add(stats_diam_sc_est_dist, diam_sc_est_dist)
        diam_kappa_dist = bin_grid_histogram_2d(diam_grid, dry_diameters, &
             kappa_grid, kappas, num_concs)
        call stats_2d_add(stats_diam_kappa_dist, diam_kappa_dist)

        ! the same for the two internal mixtures of these particles (see
        ! make_internal_mixture below), with the number concentrations of
        ! the particles, which the mixing does not change. srh_intmix is
        ! the separation relative humidity of the mixture of all
        ! particles: that of the OM/OC of all organic matter
        if (carbon_conc > 0d0) then
           srh_intmix = 1d0 / (1d0 + exp(7.7d0 * org_mass_conc &
                / carbon_conc - 15.8d0))
        else
           srh_intmix = 0d0
        end if
        call stats_1d_add_entry(stats_srh_intmix, srh_intmix, i_index)
        do i_mix = 1,N_MIX
           call make_internal_mixture(aero_state, aero_data, &
                MIX_SECTIONAL(i_mix), mixed_aero_state)
           call llps_ccn(mixed_aero_state, aero_data, env_state, &
                mix_org_masses, mix_is_mixed, mix_srhs, mix_is_separated, &
                mix_kappas, mix_scs, mix_scs_est, mix_has_film, mix_scs_film)
           if (.not. MIX_SECTIONAL(i_mix)) then
              call assert_msg(537120864, all((.not. mix_is_mixed) &
                   .or. (abs(mix_srhs - srh_intmix) < 1d-9)), &
                   "particles of the internal mixture differ in SRH")
           end if
           if (tot_num_conc > 0d0) then
              call stats_1d_add_entry(stats_mix_llps_num_frac(i_mix), &
                   sum(num_concs, mask=mix_is_separated) / tot_num_conc, &
                   i_index)
           else
              call stats_1d_add_entry(stats_mix_llps_num_frac(i_mix), 0d0, &
                   i_index)
           end if
           tot_org_mix_mass_conc = sum(mix_org_masses * num_concs)
           if (tot_org_mix_mass_conc > 0d0) then
              call stats_1d_add_entry(stats_mix_llps_org_frac(i_mix), &
                   sum(mix_org_masses * num_concs, mask=mix_is_separated) &
                   / tot_org_mix_mass_conc, i_index)
           else
              call stats_1d_add_entry(stats_mix_llps_org_frac(i_mix), 0d0, &
                   i_index)
           end if
           do i_ss = 1,N_SS
              call stats_1d_add_entry(stats_mix_ccn_conc(i_ss, i_mix), &
                   sum(num_concs, mask=(mix_scs <= CCN_SS(i_ss) / 100d0)), &
                   i_index)
              call stats_1d_add_entry(stats_mix_ccn_est_conc(i_ss, i_mix), &
                   sum(num_concs, &
                   mask=(mix_scs_est <= CCN_SS(i_ss) / 100d0)), i_index)
              do i_film = 1,N_FILM
                 call stats_1d_add_entry( &
                      stats_mix_ccn_film_conc(i_ss, i_film, i_mix), &
                      sum(num_concs, mask=(mix_scs_film(:, i_film) &
                      <= CCN_SS(i_ss) / 100d0)), i_index)
              end do
           end do
        end do

        ! the same diagnostics for the particles of each class by origin
        ! (see particle_origins below): the number distribution and
        ! concentration; the number concentrations of the particles with
        ! organic matter and salt, of the separated particles and of the
        ! CCN; the organic matter in the separated particles; and the
        ! mass concentrations of the species in the particles of the
        ! class, from which follow its dry, salt, organic and group mass
        ! concentrations and the O:C and OM:OC of its organic matter
        ! (weighted by carbon, as for all particles above)
        origins = particle_origins(aero_state, aero_data)
        if (allocated(origin_spec_mass_concs)) &
             deallocate(origin_spec_mass_concs)
        allocate(origin_spec_mass_concs(aero_data_n_spec(aero_data), &
             N_ORIGIN))
        origin_spec_mass_concs = 0d0
        do i_part = 1,aero_state_n_part(aero_state)
           origin_spec_mass_concs(:, origins(i_part)) &
                = origin_spec_mass_concs(:, origins(i_part)) &
                + num_concs(i_part) * aero_particle_species_masses( &
                aero_state%apa%particle(i_part), aero_data)
        end do
        call assert_msg(620390712, all(abs(sum(origin_spec_mass_concs, 2) &
             - spec_mass_concs) <= 1d-10 * spec_mass_concs), &
             "the classes by origin do not add up to all particles")
        if (allocated(is_dry)) deallocate(is_dry)
        if (allocated(is_salt)) deallocate(is_salt)
        allocate(is_dry(aero_data_n_spec(aero_data)), &
             is_salt(aero_data_n_spec(aero_data)))
        is_dry = .true.
        if (aero_data%i_water > 0) is_dry(aero_data%i_water) = .false.
        is_salt = .false.
        do i_name = 1,size(SALT_SPEC)
           i_spec = aero_data_spec_by_name(aero_data, SALT_SPEC(i_name))
           if (i_spec > 0) is_salt(i_spec) = .true.
        end do
        do i_origin = 1,N_ORIGIN
           in_origin = (origins == i_origin)
           call stats_1d_add(stats_origin_num_dist(i_origin), &
                bin_grid_histogram_1d(diam_grid, dry_diameters, &
                merge(num_concs, 0d0, in_origin)))
           call stats_1d_add_entry(stats_origin_num_conc(i_origin), &
                sum(num_concs, mask=in_origin), i_index)
           call stats_1d_add_entry(stats_origin_mixed_num_conc(i_origin), &
                sum(num_concs, mask=(in_origin .and. is_mixed)), i_index)
           call stats_1d_add_entry(stats_origin_llps_num_conc(i_origin), &
                sum(num_concs, mask=(in_origin .and. is_separated)), i_index)
           call stats_1d_add_entry(stats_origin_llps_org_mass_conc(i_origin), &
                sum(org_masses * num_concs, &
                mask=(in_origin .and. is_separated)), i_index)
           do i_ss = 1,N_SS
              call stats_1d_add_entry(stats_origin_ccn_conc(i_ss, i_origin), &
                   sum(num_concs, mask=(in_origin &
                   .and. (scs <= CCN_SS(i_ss) / 100d0))), i_index)
              call stats_1d_add_entry( &
                   stats_origin_ccn_est_conc(i_ss, i_origin), &
                   sum(num_concs, mask=(in_origin &
                   .and. (scs_est <= CCN_SS(i_ss) / 100d0))), i_index)
              do i_film = 1,N_FILM
                 call stats_1d_add_entry( &
                      stats_origin_ccn_film_conc(i_ss, i_film, i_origin), &
                      sum(num_concs, mask=(in_origin &
                      .and. (scs_film(:, i_film) <= CCN_SS(i_ss) / 100d0))), &
                      i_index)
              end do
           end do
           call stats_1d_add_entry(stats_origin_dry_mass_conc(i_origin), &
                sum(origin_spec_mass_concs(:, i_origin), mask=is_dry), &
                i_index)
           call stats_1d_add_entry(stats_origin_salt_mass_conc(i_origin), &
                sum(origin_spec_mass_concs(:, i_origin), mask=is_salt), &
                i_index)
           origin_org_mass_conc = sum(origin_spec_mass_concs(:, i_origin), &
                mask=(aero_data%om_oc > 0d0))
           origin_carbon_conc = sum(origin_spec_mass_concs(:, i_origin) &
                / aero_data%om_oc, mask=(aero_data%om_oc > 0d0))
           origin_oxygen_conc = sum(origin_spec_mass_concs(:, i_origin) &
                / aero_data%om_oc * aero_data%o_c, &
                mask=(aero_data%om_oc > 0d0))
           call stats_1d_add_entry(stats_origin_org_mass_conc(i_origin), &
                origin_org_mass_conc, i_index)
           if (origin_carbon_conc > 0d0) then
              call stats_1d_add_entry(stats_origin_o_c(i_origin), &
                   origin_oxygen_conc / origin_carbon_conc, i_index)
              call stats_1d_add_entry(stats_origin_om_oc(i_origin), &
                   origin_org_mass_conc / origin_carbon_conc, i_index)
           else
              call stats_1d_add_entry(stats_origin_o_c(i_origin), 0d0, &
                   i_index)
              call stats_1d_add_entry(stats_origin_om_oc(i_origin), 0d0, &
                   i_index)
           end if
           do i_group = 1,N_GROUP
              call stats_1d_add_entry( &
                   stats_origin_group_mass_conc(i_group, i_origin), &
                   group_mass_conc(aero_data, group_spec(:, i_group), &
                   origin_spec_mass_concs(:, i_origin)), i_index)
           end do
        end do

     end do

     call make_filename(out_filename, prefix, "_process.nc", index)
     write(*,*) "Writing " // trim(out_filename)
     call pmc_nc_open_write(out_filename, ncid)
     call pmc_nc_write_info(ncid, uuid, "8_organics process")
     call bin_grid_output_netcdf(diam_grid, ncid, "diam", unit="m")
     call bin_grid_output_netcdf(o_c_grid, ncid, "o_c", unit="1")
     call bin_grid_output_netcdf(frac_grid, ncid, "sec_frac", unit="1")
     call bin_grid_output_netcdf(sc_grid, ncid, "sc", unit="1")
     call bin_grid_output_netcdf(kappa_grid, ncid, "kappa", unit="1")
     call bin_grid_output_netcdf(srh_grid, ncid, "srh", unit="1")
     call bin_grid_output_netcdf(entry_grid, ncid, "entry_time", unit="s")

     call stats_1d_output_netcdf(stats_num_dist, ncid, "num_dist", &
          dim_name="diam", unit="m^{-3}")
     call stats_1d_clear(stats_num_dist)
     do i_origin = 1,N_ORIGIN
        call stats_1d_output_netcdf(stats_origin_num_dist(i_origin), ncid, &
             "num_dist_" // trim(ORIGIN_NAME(i_origin)), dim_name="diam", &
             unit="m^{-3}")
        call stats_1d_clear(stats_origin_num_dist(i_origin))
     end do

     call stats_2d_output_netcdf(stats_diam_o_c_dist, ncid, "diam_o_c_dist", &
          dim_name_1="diam", dim_name_2="o_c", unit="m^{-3}")
     call stats_2d_clear(stats_diam_o_c_dist)

     call stats_2d_output_netcdf(stats_diam_sec_frac_dist, ncid, &
          "diam_sec_frac_dist", dim_name_1="diam", dim_name_2="sec_frac", &
          unit="m^{-3}")
     call stats_2d_clear(stats_diam_sec_frac_dist)

     call stats_2d_output_netcdf(stats_diam_srh_dist, ncid, "diam_srh_dist", &
          dim_name_1="diam", dim_name_2="srh", unit="m^{-3}")
     call stats_2d_clear(stats_diam_srh_dist)

     call stats_2d_output_netcdf(stats_entry_srh_dist, ncid, &
          "entry_time_srh_dist", dim_name_1="entry_time", dim_name_2="srh", &
          unit="m^{-3}")
     call stats_2d_clear(stats_entry_srh_dist)

     call stats_2d_output_netcdf(stats_diam_sc_dist, ncid, "diam_sc_dist", &
          dim_name_1="diam", dim_name_2="sc", unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_dist)

     call stats_2d_output_netcdf(stats_diam_sc_est_dist, ncid, &
          "diam_sc_est_dist", dim_name_1="diam", dim_name_2="sc", &
          unit="m^{-3}")
     call stats_2d_clear(stats_diam_sc_est_dist)

     call stats_2d_output_netcdf(stats_diam_kappa_dist, ncid, &
          "diam_kappa_dist", dim_name_1="diam", dim_name_2="kappa", &
          unit="m^{-3}")
     call stats_2d_clear(stats_diam_kappa_dist)

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
  call stats_1d_output_netcdf(stats_rel_humid, ncid, "rel_humid", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_mixed_num_frac, ncid, &
       "mixed_num_frac", dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_llps_num_frac, ncid, "llps_num_frac", &
       dim_name="time", unit="1")
  call stats_1d_output_netcdf(stats_llps_org_frac, ncid, &
       "llps_org_mass_frac", dim_name="time", unit="1")
  do i_ss = 1,N_SS
     call stats_1d_output_netcdf(stats_ccn_conc(i_ss), ncid, &
          "ccn_conc_" // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
     call stats_1d_output_netcdf(stats_ccn_est_conc(i_ss), ncid, &
          "ccn_est_conc_" // CCN_SS_NAME(i_ss), dim_name="time", &
          unit="m^{-3}")
  end do
  do i_film = 1,N_FILM
     call stats_1d_output_netcdf(stats_film_num_frac(i_film), ncid, &
          "film_num_frac_" // FILM_NAME(i_film), dim_name="time", unit="1")
     do i_ss = 1,N_SS
        call stats_1d_output_netcdf(stats_ccn_film_conc(i_ss, i_film), &
             ncid, "ccn_film_conc_" // FILM_NAME(i_film) // "_" &
             // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
     end do
  end do
  call stats_1d_output_netcdf(stats_srh_intmix, ncid, "srh_intmix", &
       dim_name="time", unit="1")
  do i_mix = 1,N_MIX
     call stats_1d_output_netcdf(stats_mix_llps_num_frac(i_mix), ncid, &
          "llps_num_frac_" // trim(MIX_NAME(i_mix)), dim_name="time", &
          unit="1")
     call stats_1d_output_netcdf(stats_mix_llps_org_frac(i_mix), ncid, &
          "llps_org_mass_frac_" // trim(MIX_NAME(i_mix)), &
          dim_name="time", unit="1")
     do i_ss = 1,N_SS
        call stats_1d_output_netcdf(stats_mix_ccn_conc(i_ss, i_mix), ncid, &
             "ccn_conc_" // trim(MIX_NAME(i_mix)) // "_" &
             // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
        call stats_1d_output_netcdf(stats_mix_ccn_est_conc(i_ss, i_mix), &
             ncid, "ccn_est_conc_" // trim(MIX_NAME(i_mix)) // "_" &
             // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
        do i_film = 1,N_FILM
           call stats_1d_output_netcdf( &
                stats_mix_ccn_film_conc(i_ss, i_film, i_mix), ncid, &
                "ccn_film_conc_" // trim(MIX_NAME(i_mix)) // "_" &
                // FILM_NAME(i_film) // "_" // CCN_SS_NAME(i_ss), &
                dim_name="time", unit="m^{-3}")
        end do
     end do
  end do
  do i_origin = 1,N_ORIGIN
     call stats_1d_output_netcdf(stats_origin_num_conc(i_origin), ncid, &
          "num_conc_" // trim(ORIGIN_NAME(i_origin)), dim_name="time", &
          unit="m^{-3}")
     call stats_1d_output_netcdf(stats_origin_mixed_num_conc(i_origin), &
          ncid, "mixed_num_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="m^{-3}")
     call stats_1d_output_netcdf(stats_origin_llps_num_conc(i_origin), &
          ncid, "llps_num_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="m^{-3}")
     call stats_1d_output_netcdf(stats_origin_llps_org_mass_conc(i_origin), &
          ncid, "llps_org_mass_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="kg m^{-3}")
     do i_ss = 1,N_SS
        call stats_1d_output_netcdf(stats_origin_ccn_conc(i_ss, i_origin), &
             ncid, "ccn_conc_" // trim(ORIGIN_NAME(i_origin)) // "_" &
             // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
        call stats_1d_output_netcdf( &
             stats_origin_ccn_est_conc(i_ss, i_origin), ncid, &
             "ccn_est_conc_" // trim(ORIGIN_NAME(i_origin)) // "_" &
             // CCN_SS_NAME(i_ss), dim_name="time", unit="m^{-3}")
        do i_film = 1,N_FILM
           call stats_1d_output_netcdf( &
                stats_origin_ccn_film_conc(i_ss, i_film, i_origin), ncid, &
                "ccn_film_conc_" // trim(ORIGIN_NAME(i_origin)) // "_" &
                // FILM_NAME(i_film) // "_" // CCN_SS_NAME(i_ss), &
                dim_name="time", unit="m^{-3}")
        end do
     end do
     call stats_1d_output_netcdf(stats_origin_dry_mass_conc(i_origin), &
          ncid, "dry_mass_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="kg m^{-3}")
     call stats_1d_output_netcdf(stats_origin_salt_mass_conc(i_origin), &
          ncid, "salt_mass_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="kg m^{-3}")
     call stats_1d_output_netcdf(stats_origin_org_mass_conc(i_origin), &
          ncid, "org_mass_conc_" // trim(ORIGIN_NAME(i_origin)), &
          dim_name="time", unit="kg m^{-3}")
     do i_group = 1,N_GROUP
        call stats_1d_output_netcdf( &
             stats_origin_group_mass_conc(i_group, i_origin), ncid, &
             "mass_conc_" // trim(GROUP_NAME(i_group)) // "_" &
             // trim(ORIGIN_NAME(i_origin)), dim_name="time", &
             unit="kg m^{-3}")
     end do
     call stats_1d_output_netcdf(stats_origin_o_c(i_origin), ncid, &
          "o_c_" // trim(ORIGIN_NAME(i_origin)), dim_name="time", unit="1")
     call stats_1d_output_netcdf(stats_origin_om_oc(i_origin), ncid, &
          "om_oc_" // trim(ORIGIN_NAME(i_origin)), dim_name="time", &
          unit="1")
  end do
  call pmc_nc_close(ncid)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Liquid-liquid phase separation and critical supersaturation of
  !> every particle of an aerosol state.
  !!
  !! The phase separation is diagnosed as in CMAQ (Pye et al., 2017,
  !! Eq. 7): the separation relative humidity of a particle follows from
  !! the OM/OC of its organic matter, SRH = 1 / (1 + exp(7.7 OM/OC -
  !! 15.8)), the relation that Pye et al. derived from the experiments
  !! with ammonium sulfate of You et al. (2013); a particle that contains
  !! organic matter and inorganic salt separates into an organic-rich and
  !! an inorganic-rich phase when the relative humidity is below its SRH.
  !!
  !! The critical supersaturation follows from kappa-Koehler theory, with
  !! the solute kappa of the dry particle (the volume average over the
  !! species; for the organic species 0.11 OM:OC - 0.10 as in CMAQ), once
  !! with the surface tension of water and once with the effective
  !! surface tension of an organic film
  !! (aero_particle_crit_rel_humid_est). For particles without soluble
  !! material the effective surface tension is not defined (that function
  !! then returns a critical relative humidity of 1), so they keep the
  !! value with the surface tension of water.
  !!
  !! These two critical supersaturations are limits: no particle and
  !! every particle carries the organic film until it activates. The
  !! film is that of the models with complete phase separation of Davies
  !! et al. (2019, after Ovadnevaite et al., 2017), in which the organic
  !! matter stays apart from the aqueous phase whatever its composition;
  !! in their equilibrium model the film is gone once the organic matter
  !! mixes with the aqueous phase. The separation relative humidity
  !! above does not tell which particles keep the film, because it is
  !! below 1 for every OM/OC and the particles activate above 1. Here
  !! the particles keep the film whose organic matter has an O:C (the
  !! average of the particle, weighted by carbon) of at most a limit;
  !! the others take the value with the surface tension of water. This
  !! is done for each limit of FILM_O_C: O:C up to which two liquid
  !! phases were observed near saturation or at supersaturation, in
  !! secondary organic material without salts (Renbaum-Wolff et al.,
  !! 2016; Song et al., 2017; Ham et al., 2019), in particles of one or
  !! two organic compounds (Song et al., 2018, 2020) and in mixtures
  !! with ammonium sulfate (Malek et al., 2023).
  subroutine llps_ccn(aero_state, aero_data, env_state, org_masses, &
       is_mixed, srhs, is_separated, kappas, scs, scs_est, has_film, &
       scs_film)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Mass of the organic matter of each particle (kg).
    real(kind=dp), allocatable, intent(inout) :: org_masses(:)
    !> Whether each particle contains organic matter and inorganic salt.
    logical, allocatable, intent(inout) :: is_mixed(:)
    !> Separation relative humidity of each particle (1), 0 for the
    !> particles without organic matter or salt.
    real(kind=dp), allocatable, intent(inout) :: srhs(:)
    !> Whether each particle is separated into two liquid phases.
    logical, allocatable, intent(inout) :: is_separated(:)
    !> Solute kappa of each particle (1).
    real(kind=dp), allocatable, intent(inout) :: kappas(:)
    !> Critical supersaturation of each particle with the surface tension
    !> of water (1).
    real(kind=dp), allocatable, intent(inout) :: scs(:)
    !> Critical supersaturation of each particle with the effective
    !> surface tension (1).
    real(kind=dp), allocatable, intent(inout) :: scs_est(:)
    !> Whether each particle keeps the organic film until it activates,
    !> for each limit of FILM_O_C (particle, limit).
    logical, allocatable, intent(inout) :: has_film(:,:)
    !> Critical supersaturation of each particle with the effective
    !> surface tension if it keeps the film and with the surface tension
    !> of water if not, for each limit of FILM_O_C (particle, limit) (1).
    real(kind=dp), allocatable, intent(inout) :: scs_film(:,:)

    integer :: n_part, i_part, i_salt, i_spec, i_film
    real(kind=dp), allocatable :: org_carbon_masses(:), salt_masses(:), &
         o_cs(:)

    n_part = aero_state_n_part(aero_state)
    org_masses = aero_state_organic_masses(aero_state, aero_data)
    org_carbon_masses = aero_state_organic_carbon_masses(aero_state, &
         aero_data)
    allocate(salt_masses(n_part))
    salt_masses = 0d0
    do i_salt = 1,size(SALT_SPEC)
       i_spec = aero_data_spec_by_name(aero_data, SALT_SPEC(i_salt))
       if (i_spec == 0) cycle
       salt_masses = salt_masses + aero_particle_species_mass( &
            aero_state%apa%particle(1:n_part), i_spec, aero_data)
    end do
    if (allocated(srhs)) deallocate(srhs)
    if (allocated(kappas)) deallocate(kappas)
    if (allocated(is_mixed)) deallocate(is_mixed)
    if (allocated(is_separated)) deallocate(is_separated)
    allocate(srhs(n_part), kappas(n_part), is_mixed(n_part), &
         is_separated(n_part))
    do i_part = 1,n_part
       is_mixed(i_part) = (org_carbon_masses(i_part) > 0d0) &
            .and. (salt_masses(i_part) > 0d0)
       if (is_mixed(i_part)) then
          srhs(i_part) = 1d0 / (1d0 + exp(7.7d0 * org_masses(i_part) &
               / org_carbon_masses(i_part) - 15.8d0))
       else
          srhs(i_part) = 0d0
       end if
       is_separated(i_part) = is_mixed(i_part) &
            .and. (env_state%rel_humid < srhs(i_part))
    end do

    scs = aero_state_crit_rel_humids(aero_state, aero_data, env_state) &
         - 1d0
    scs_est = aero_state_crit_rel_humids_est(aero_state, aero_data, &
         env_state) - 1d0
    do i_part = 1,n_part
       kappas(i_part) = aero_particle_solute_kappa( &
            aero_state%apa%particle(i_part), aero_data)
       if (kappas(i_part) < 1d-30) scs_est(i_part) = scs(i_part)
    end do

    o_cs = aero_state_o_c(aero_state, aero_data)
    if (allocated(has_film)) deallocate(has_film)
    if (allocated(scs_film)) deallocate(scs_film)
    allocate(has_film(n_part, N_FILM), scs_film(n_part, N_FILM))
    do i_film = 1,N_FILM
       has_film(:, i_film) = (org_carbon_masses > 0d0) &
            .and. (o_cs <= FILM_O_C(i_film))
       scs_film(:, i_film) = merge(scs_est, scs, has_film(:, i_film))
    end do

  end subroutine llps_ccn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Internal mixture of the particles of an aerosol state.
  !!
  !! The dry species of every particle are replaced by the average dry
  !! composition of all particles or, if \c sectional, of the particles
  !! in its size section (the volume fractions of the dry species,
  !! weighted by the number concentrations of the particles), scaled to
  !! the dry volume of the particle; its water is kept. So the dry
  !! diameter and the total volume of every particle (and with it its
  !! number concentration) and the total concentration of every species
  !! are those of the original particles, and only the composition is
  !! averaged.
  subroutine make_internal_mixture(aero_state, aero_data, sectional, &
       mixed_aero_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to mix the particles of each size section separately.
    logical, intent(in) :: sectional
    !> Internal mixture of the particles of \c aero_state.
    type(aero_state_t), intent(inout) :: mixed_aero_state

    integer :: n_part, i_part, i_sect
    integer, allocatable :: part_sects(:)
    real(kind=dp) :: dry_diam, sect_vol
    real(kind=dp), allocatable :: num_concs(:), mix_num_concs(:), vols(:), &
         sect_fracs(:,:), spec_concs(:), mix_spec_concs(:)

    n_part = aero_state_n_part(aero_state)
    mixed_aero_state = aero_state
    num_concs = aero_state_num_concs(aero_state, aero_data)

    ! number-weighted dry species volumes of each section, then their
    ! fractions
    allocate(part_sects(n_part), &
         sect_fracs(aero_data_n_spec(aero_data), N_SECT))
    sect_fracs = 0d0
    do i_part = 1,n_part
       i_sect = 1
       if (sectional) then
          dry_diam = aero_particle_dry_diameter( &
               aero_state%apa%particle(i_part), aero_data)
          i_sect = floor(real(N_SECT, kind=dp) &
               * log(dry_diam / SECT_DIAM_MIN) &
               / log(SECT_DIAM_MAX / SECT_DIAM_MIN)) + 1
          i_sect = max(1, min(N_SECT, i_sect))
       end if
       part_sects(i_part) = i_sect
       vols = aero_state%apa%particle(i_part)%vol
       if (aero_data%i_water > 0) vols(aero_data%i_water) = 0d0
       sect_fracs(:, i_sect) = sect_fracs(:, i_sect) &
            + num_concs(i_part) * vols
    end do
    do i_sect = 1,N_SECT
       sect_vol = sum(sect_fracs(:, i_sect))
       if (sect_vol > 0d0) then
          sect_fracs(:, i_sect) = sect_fracs(:, i_sect) / sect_vol
       end if
    end do

    do i_part = 1,n_part
       vols = sect_fracs(:, part_sects(i_part)) &
            * aero_particle_dry_volume(aero_state%apa%particle(i_part), &
            aero_data)
       if (aero_data%i_water > 0) then
          vols(aero_data%i_water) &
               = aero_state%apa%particle(i_part)%vol(aero_data%i_water)
       end if
       call aero_particle_set_vols(mixed_aero_state%apa%particle(i_part), &
            vols)
    end do

    ! the number concentrations and the total concentration of every
    ! species must not change
    mix_num_concs = aero_state_num_concs(mixed_aero_state, aero_data)
    call assert_msg(270554861, &
         all(abs(mix_num_concs - num_concs) <= 1d-10 * num_concs), &
         "internal mixing changed number concentrations")
    allocate(spec_concs(aero_data_n_spec(aero_data)), &
         mix_spec_concs(aero_data_n_spec(aero_data)))
    spec_concs = 0d0
    mix_spec_concs = 0d0
    do i_part = 1,n_part
       spec_concs = spec_concs + num_concs(i_part) &
            * aero_state%apa%particle(i_part)%vol
       mix_spec_concs = mix_spec_concs + num_concs(i_part) &
            * mixed_aero_state%apa%particle(i_part)%vol
    end do
    call assert_msg(841663752, &
         all(abs(mix_spec_concs - spec_concs) <= 1d-9 * spec_concs), &
         "internal mixing changed species concentrations")

  end subroutine make_internal_mixture

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Class of every particle by the origin of its constituents.
  !!
  !! PartMC keeps for every particle the source and the creation time of
  !! each primary particle that it is made of (its components, which
  !! coagulation joins). The particles of the initial air are created at
  !! time 0, those that the exchange with the background air brings in
  !! at the time they enter. A particle is in the class
  !! - 1 (city) if all its components are urban primary particles (the
  !!   sources other than REGIONAL_SOURCE);
  !! - 2 (city_regional) if it has urban and regional components;
  !! - 3 (regional) if all its components are regional particles and one
  !!   of them is from the initial air;
  !! - 4 (background) if all its components are regional particles that
  !!   entered after the start.
  function particle_origins(aero_state, aero_data)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    !> Class of every particle.
    integer :: particle_origins(aero_state_n_part(aero_state))

    integer :: i_part, i_comp, i_regional_source
    logical :: has_city, has_initial, has_entered

    i_regional_source = 0
    if (aero_data_n_source(aero_data) > 0) then
       i_regional_source = string_array_find(aero_data%source_name, &
            REGIONAL_SOURCE)
    end if
    do i_part = 1,aero_state_n_part(aero_state)
       call assert_msg(397166208, aero_particle_n_components( &
            aero_state%apa%particle(i_part)) >= 1, &
            "particle without components")
       has_city = .false.
       has_initial = .false.
       has_entered = .false.
       do i_comp = 1,aero_particle_n_components( &
            aero_state%apa%particle(i_part))
          if (aero_state%apa%particle(i_part)%component(i_comp)%source_id &
               /= i_regional_source) then
             has_city = .true.
          elseif (aero_state%apa%particle(i_part)%component(i_comp) &
               %create_time > 0d0) then
             has_entered = .true.
          else
             has_initial = .true.
          end if
       end do
       if (has_city .and. (has_initial .or. has_entered)) then
          particle_origins(i_part) = 2
       elseif (has_city) then
          particle_origins(i_part) = 1
       elseif (has_initial) then
          particle_origins(i_part) = 3
       else
          particle_origins(i_part) = 4
       end if
    end do

  end function particle_origins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mass concentration of a group of species (kg m^{-3}): the sum over
  !> the species of the group that the aerosol data contain.
  real(kind=dp) function group_mass_conc(aero_data, names, spec_mass_concs)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Names of the species of the group (empty names are skipped).
    character(len=*), intent(in) :: names(:)
    !> Mass concentration of every species (kg m^{-3}).
    real(kind=dp), intent(in) :: spec_mass_concs(:)

    integer :: i_name, i_spec

    group_mass_conc = 0d0
    do i_name = 1,size(names)
       if (len_trim(names(i_name)) == 0) cycle
       i_spec = aero_data_spec_by_name(aero_data, names(i_name))
       if (i_spec == 0) cycle
       group_mass_conc = group_mass_conc + spec_mass_concs(i_spec)
    end do

  end function group_mass_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program process
