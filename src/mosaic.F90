! Copyright (C) 2007-2012, 2016 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_mosaic module.

!> Interface to the MOSAIC aerosol and gas phase chemistry code.
module pmc_mosaic

  use pmc_aero_data
  use pmc_aero_state
  use pmc_constants
  use pmc_env_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_output
  use pmc_spec_file
  use pmc_util

  !> Number of photolysis reactions of MOSAIC.
  integer, parameter :: MOSAIC_N_PHOTOLYSIS = 20
  !> Length of the names of the photolysis reactions.
  integer, parameter :: MOSAIC_PHOTOLYSIS_NAME_LEN = 10
  !> Names of the photolysis reactions of MOSAIC, in the order of the
  !> rates in the \c photolysis_rates argument of mosaic_init().
  !!
  !! Ozone and formaldehyde have two channels each: \c O3_O3P and \c
  !! O3_O1D give O(3P) and O(1D), \c HCHO_R gives the radicals (2 HO2
  !! + CO) and \c HCHO_M the molecular products (CO + H2).
  character(len=MOSAIC_PHOTOLYSIS_NAME_LEN), parameter :: &
       MOSAIC_PHOTOLYSIS_NAME(MOSAIC_N_PHOTOLYSIS) = &
       [character(len=MOSAIC_PHOTOLYSIS_NAME_LEN) :: &
       "NO2",    "NO3",    "HONO",   "HNO3",   "HNO4", &
       "N2O5",   "O3_O3P", "O3_O1D", "H2O2",   "HCHO_R", &
       "HCHO_M", "CH3OOH", "ETHOOH", "ALD2",   "AONE", &
       "MGLY",   "OPEN",   "ROOH",   "ONIT",   "ISOPRD"]

  !> Acidity of the acid-catalyzed IEPOX uptake: MOSAIC's own acidity of
  !> each particle.
  integer, parameter :: MOSAIC_IEPOX_ACIDITY_MOSAIC = 0
  !> Acidity of the acid-catalyzed IEPOX uptake: one ISORROPIA II
  !> calculation on the totals of gases and particles, as in CMAQ.
  integer, parameter :: MOSAIC_IEPOX_ACIDITY_ISORROPIA = 1

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether MOSAIC support is compiled in.
  logical function mosaic_support()

#ifdef PMC_USE_MOSAIC
    mosaic_support = .true.
#else
    mosaic_support = .false.
#endif

  end function mosaic_support

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize all MOSAIC data-structures.
  subroutine mosaic_init(env_state, aero_data, del_t, do_optical, &
       photolysis_rates, iepox_acidity)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: alpha_ASTEM, rtol_eqb_ASTEM, &
         ptol_mol_ASTEM, mGAS_AER_XFER, mDYNAMIC_SOLVER

    use module_data_mosaic_main, only: tbeg_sec, dt_sec, rlon, rlat, &
         zalt_m, RH, te, pr_atm, cair_mlc, cair_molm3, ppb, avogad, &
         mmode, mgas, maer, mcld, maeroptic, mshellcore, &
         msolar, mphoto, lun_aeroptic, naerbin

    use module_data_mosaic_gas, only: nphoto, rk_photo_fixed, &
         jphoto_no2, jphoto_no3, jphoto_hono, jphoto_hno3, jphoto_hno4, &
         jphoto_n2o5, jphoto_o3a, jphoto_o3b, jphoto_h2o2, jphoto_hchoa, &
         jphoto_hchob, jphoto_ch3ooh, jphoto_ethooh, jphoto_ald2, &
         jphoto_aone, jphoto_mgly, jphoto_open, jphoto_rooh, jphoto_onit, &
         jphoto_isoprd

    use module_data_mosaic_soa, only: nsoa, nsoa_sv, soa_name_a, soa_mw, &
         nsoa_prxn, soa_prxn_from, soa_prxn_to, lsv_pcso, lnv_pcso2, &
         soa_iepox_acidity
#endif

    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Timestep for coagulation.
    real(kind=dp), intent(in) :: del_t
    !> Whether to compute optical properties.
    logical, intent(in) :: do_optical
    !> Fixed photolysis rates (s^{-1}) in the order of
    !> MOSAIC_PHOTOLYSIS_NAME, as in a chamber with artificial
    !> light. If not present, the rates follow the solar zenith angle.
    real(kind=dp), intent(in), optional :: &
         photolysis_rates(MOSAIC_N_PHOTOLYSIS)
    !> Acidity of the acid-catalyzed IEPOX uptake
    !> (MOSAIC_IEPOX_ACIDITY_ISORROPIA or MOSAIC_IEPOX_ACIDITY_MOSAIC);
    !> MOSAIC's default (ISORROPIA) if not present.
    integer, intent(in), optional :: iepox_acidity

#ifdef PMC_USE_MOSAIC
    integer :: i_soa, i_spec, i_rxn
    character(len=8) :: rxn_from(nsoa_prxn + 1), rxn_to(nsoa_prxn + 1)

    ! MOSAIC function interfaces
    interface
       subroutine LoadPeroxyParameters()
       end subroutine LoadPeroxyParameters
       subroutine init_data_modules()
       end subroutine init_data_modules
       subroutine AllocateMemory()
       end subroutine AllocateMemory
    end interface

    call init_data_modules  ! initialize indices and vars

    if (present(iepox_acidity)) then
       call assert_msg(218446017, &
            (iepox_acidity == MOSAIC_IEPOX_ACIDITY_MOSAIC) &
            .or. (iepox_acidity == MOSAIC_IEPOX_ACIDITY_ISORROPIA), &
            "unknown acidity of the IEPOX uptake")
       soa_iepox_acidity = iepox_acidity
    end if

    ! allocate one aerosol bin
    naerbin = 1
    call AllocateMemory()

    ! parameters
    mmode = 1               ! 1 = time integration, 2 = parametric analysis
    mgas = 1                ! 1 = gas chem on, 0 = gas chem off
    maer = 1                ! 1 = aer chem on, 0 = aer chem off
    mcld = 0                ! 1 = cld chem on, 0 = cld chem off
    if (do_optical) then
       maeroptic = 1        ! 1 = aer_optical on, 0 = aer_optical off
    else
       maeroptic = 0
    end if
    mshellcore = 1          ! 0 = no shellcore, 1 = core is BC only
                            ! 2 = core is BC and DUST
    if (present(photolysis_rates)) then
       msolar = 2           ! 1 = diurnally varying phot, 2 = fixed phot
       call assert(742360195, nphoto == MOSAIC_N_PHOTOLYSIS)
       rk_photo_fixed(jphoto_no2) = photolysis_rates(1)
       rk_photo_fixed(jphoto_no3) = photolysis_rates(2)
       rk_photo_fixed(jphoto_hono) = photolysis_rates(3)
       rk_photo_fixed(jphoto_hno3) = photolysis_rates(4)
       rk_photo_fixed(jphoto_hno4) = photolysis_rates(5)
       rk_photo_fixed(jphoto_n2o5) = photolysis_rates(6)
       rk_photo_fixed(jphoto_o3a) = photolysis_rates(7)
       rk_photo_fixed(jphoto_o3b) = photolysis_rates(8)
       rk_photo_fixed(jphoto_h2o2) = photolysis_rates(9)
       rk_photo_fixed(jphoto_hchoa) = photolysis_rates(10)
       rk_photo_fixed(jphoto_hchob) = photolysis_rates(11)
       rk_photo_fixed(jphoto_ch3ooh) = photolysis_rates(12)
       rk_photo_fixed(jphoto_ethooh) = photolysis_rates(13)
       rk_photo_fixed(jphoto_ald2) = photolysis_rates(14)
       rk_photo_fixed(jphoto_aone) = photolysis_rates(15)
       rk_photo_fixed(jphoto_mgly) = photolysis_rates(16)
       rk_photo_fixed(jphoto_open) = photolysis_rates(17)
       rk_photo_fixed(jphoto_rooh) = photolysis_rates(18)
       rk_photo_fixed(jphoto_onit) = photolysis_rates(19)
       rk_photo_fixed(jphoto_isoprd) = photolysis_rates(20)
    else
       msolar = 1
    end if
    mphoto = 2              ! 1 = Rick's param, 2 = Yang's param
    mGAS_AER_XFER = 1       ! 1 = gas-aerosol partitioning, 0 = no partition
    mDYNAMIC_SOLVER = 1     ! 1 = astem, 2 = lsodes
    alpha_ASTEM = 0.5d0     ! solver parameter. range: 0.01 - 1.0
    rtol_eqb_ASTEM = 0.01d0 ! relative eqb tolerance. range: 0.01 - 0.03
    ptol_mol_ASTEM = 0.01d0 ! percent mol tolerance.  range: 0.01 - 1.0

    ! time variables
    dt_sec = del_t                                 ! time-step (s)
    tbeg_sec = env_state%start_day*24*3600 + &     ! time since the beg of
         nint(env_state%start_time)                ! year 00:00, UTC (s)

    ! geographic location
    rlon = deg2rad(env_state%longitude)            ! longitude
    rlat = deg2rad(env_state%latitude)             ! latitude
    zalt_m = env_state%altitude                    ! altitude (m)

    ! environmental parameters: map PartMC -> MOSAIC
    RH = env_state%rel_humid * 100.d0              ! relative humidity (%)
    te = env_state%temp                            ! temperature (K)
    pr_atm = env_state%pressure / const%air_std_press ! pressure (atm)
    cair_mlc = avogad*pr_atm/(82.056d0*te)         ! air conc [molec/cc]
    cair_molm3 = 1d6*pr_atm/(82.056d0*te)          ! air conc [mol/m^3]
    ppb = 1d9

    call LoadPeroxyParameters ! Aperox and Bperox only once

    ! get unit for aerosol optical output
    if (lun_aeroptic <= 0 ) lun_aeroptic = get_unit()

    ! ensure H2O is a valid species
    call assert_msg(111041803, aero_data%i_water > 0, &
         "MOSAIC requires H2O as an aerosol species")

    ! the saturation concentrations of the organic species are given per
    ! mass and are converted to moles with the MOSAIC molecular weights,
    ! so the molecular weights in aero_data must be the same
    do i_soa = 1,nsoa
       i_spec = aero_data_spec_by_name(aero_data, trim(soa_name_a(i_soa)))
       if (i_spec > 0) then
          call assert_msg(280416507, almost_equal( &
               aero_data%molec_weight(i_spec), soa_mw(i_soa) * 1d-3), &
               "molecular weight of aerosol species " &
               // trim(aero_data%name(i_spec)) // " (" &
               // trim(real_to_string(aero_data%molec_weight(i_spec))) &
               // " kg/mol) differs from the MOSAIC value (" &
               // trim(real_to_string(soa_mw(i_soa) * 1d-3)) // " kg/mol)")
       end if
    end do

    ! the products of the particle-phase reactions of the organic species
    ! (oligomerization, nitrate hydrolysis, aging of potential combustion
    ! SOA) must be aerosol species too, or the mass that reacts is lost
    rxn_from(1:nsoa_prxn) = soa_name_a(soa_prxn_from)
    rxn_to(1:nsoa_prxn) = soa_name_a(nsoa_sv + soa_prxn_to)
    rxn_from(nsoa_prxn + 1) = soa_name_a(lsv_pcso)
    rxn_to(nsoa_prxn + 1) = soa_name_a(nsoa_sv + lnv_pcso2)
    do i_rxn = 1,nsoa_prxn + 1
       if ((aero_data_spec_by_name(aero_data, trim(rxn_from(i_rxn))) > 0) &
            .and. (aero_data_spec_by_name(aero_data, trim(rxn_to(i_rxn))) &
            <= 0)) then
          call warn_msg(518263309, "aerosol species " &
               // trim(rxn_from(i_rxn)) // " reacts to " &
               // trim(rxn_to(i_rxn)) // " in MOSAIC, which is not an " &
               // "aerosol species: the mass that reacts is lost")
       end if
    end do

#endif

  end subroutine mosaic_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the fixed photolysis rates of MOSAIC from the file.
  subroutine spec_file_read_mosaic_photolysis(file, photolysis_rates)

    !> File to read the rates from.
    type(spec_file_t), intent(inout) :: file
    !> Photolysis rates (s^{-1}) in the order of MOSAIC_PHOTOLYSIS_NAME.
    real(kind=dp), allocatable, intent(inout) :: photolysis_rates(:)

    integer :: n_line, i_line, i_reaction
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: reaction_name(:)
    real(kind=dp), allocatable :: reaction_data(:,:)

    !> \page input_format_mosaic_photolysis Input File Format: MOSAIC Photolysis
    !!
    !! A file with fixed photolysis rates describes artificial light
    !! that does not change with time, as in a chamber experiment. It
    !! must consist of one line per photolysis reaction, with each line
    !! having the reaction name followed by the photolysis rate in
    !! s^{-1}. The valid reaction names are those of the
    !! photolyzed species of MOSAIC: \c NO2, \c NO3, \c HONO, \c HNO3,
    !! \c HNO4, \c N2O5, \c O3_O3P, \c O3_O1D, \c H2O2, \c HCHO_R,
    !! \c HCHO_M, \c CH3OOH, \c ETHOOH, \c ALD2, \c AONE, \c MGLY,
    !! \c OPEN, \c ROOH, \c ONIT and \c ISOPRD. Ozone and formaldehyde
    !! have two channels each: \c O3_O3P and \c O3_O1D give O(3P) and
    !! O(1D), \c HCHO_R gives the radical products (2 HO2 + CO) and \c
    !! HCHO_M the molecular products (CO + H2). Not all reactions have
    !! to be listed. Any missing reactions will have rates of zero. For
    !! example, a photolysis rate file could contain:
    !! <pre>
    !! # reaction  rate (s^{-1})
    !! NO2         5.0e-3
    !! HONO        1.1e-3
    !! H2O2        2.9e-6
    !! </pre>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_real_named_array(file, 0, reaction_name, &
         reaction_data)

    ! check the data size
    n_line = size(reaction_data, 1)
    if (.not. ((size(reaction_data, 2) == 1) .or. (n_line == 0))) then
       call die_msg(473826501, 'each line in ' // trim(file%name) &
            // ' must contain exactly one data value')
    end if

    ! copy over the data
    if (allocated(photolysis_rates)) deallocate(photolysis_rates)
    allocate(photolysis_rates(MOSAIC_N_PHOTOLYSIS))
    photolysis_rates = 0d0
    do i_line = 1,n_line
       i_reaction = string_array_find(MOSAIC_PHOTOLYSIS_NAME, &
            reaction_name(i_line))
       if (i_reaction == 0) then
          call die_msg(605719284, 'unknown photolysis reaction ' // &
               trim(reaction_name(i_line)) // ' in file ' &
               // trim(file%name))
       end if
       if (reaction_data(i_line,1) < 0d0) then
          call die_msg(318840657, 'negative photolysis rate of ' // &
               trim(reaction_name(i_line)) // ' in file ' &
               // trim(file%name))
       end if
       photolysis_rates(i_reaction) = reaction_data(i_line,1)
    end do

  end subroutine spec_file_read_mosaic_photolysis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clean-up after running MOSAIC, deallocating memory.
  subroutine mosaic_cleanup()

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine DeallocateMemory()
       end subroutine DeallocateMemory
    end interface

    call DeallocateMemory()
#endif

  end subroutine mosaic_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data PartMC -> MOSAIC.
  subroutine mosaic_from_partmc(env_state, aero_data, &
       aero_state, gas_data, gas_state)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: nbin_a, aer, num_a, jhyst_leg, &
#ifdef PMC_USE_WRF
         particle_error, &
#endif
         jtotal, water_a

    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, RH, te, pr_atm, cnn, cair_mlc, &
         cair_molm3, ppb, avogad, msolar, naerbin
    use module_data_mosaic_gas, only: photo_factor
#endif

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state

#ifdef PMC_USE_MOSAIC
    ! local variables
    real(kind=dp) :: time_UTC    ! 24-hr UTC clock time (hr).
    real(kind=dp) :: tmar21_sec  ! Time at noon, march 21, UTC (s).
    real(kind=dp) :: conv_fac(aero_data_n_spec(aero_data)), dum_var
    integer :: i_part, i_spec, i_spec_mosaic
    real(kind=dp) :: num_conc

    ! MOSAIC function interfaces
    interface
       subroutine AllocateMemory()
       end subroutine AllocateMemory
       subroutine DeallocateMemory()
       end subroutine DeallocateMemory
    end interface

    ! update time variables
    tmar21_sec = real((79*24 + 12)*3600, kind=dp)    ! noon, mar 21, UTC
    tcur_sec = real(tbeg_sec, kind=dp) + env_state%elapsed_time
    ! current (old) time since the beg of year 00:00, UTC (s)

    time_UTC = env_state%start_time/3600d0  ! 24-hr UTC clock time (hr)
    time_UTC = time_UTC + dt_sec/3600d0

    do while (time_UTC >= 24d0)
       time_UTC = time_UTC - 24d0
    end do

    tmid_sec = tcur_sec + 0.5d0*dt_sec
    if(tmid_sec .ge. tmar21_sec)then
       tmid_sec = tmid_sec - tmar21_sec     ! seconds since noon, march 21
    else
       tmid_sec = tmid_sec &
            + dble(((365-79)*24 - 12)*3600) ! seconds since noon, march 21
    endif

    ! transport timestep (min)
    dt_min = dt_sec/60d0
    ! aerosol optics timestep (min)
    dt_aeroptic_min = 0d0

    ! compute aerosol conversion factors
    do i_spec = 1,aero_data_n_spec(aero_data)
       ! converts m^3(species) to nmol(species)/m^3(air)
       conv_fac(i_spec) = 1.D9 * aero_data%density(i_spec) &
            / aero_data%molec_weight(i_spec)
    enddo

    ! environmental parameters: map PartMC -> MOSAIC
    RH = env_state%rel_humid * 100.d0              ! relative humidity (%)
    te = env_state%temp                            ! temperature (K)
    pr_atm = env_state%pressure / const%air_std_press ! pressure (atm)
    cair_mlc = avogad*pr_atm/(82.056d0*te)   ! air conc [molec/cc]
    cair_molm3 = 1d6*pr_atm/(82.056d0*te)    ! air conc [mol/m^3]
    ppb = 1d9
    photo_factor = env_state%photolysis_factor ! factor on all photolysis rates

    ! aerosol data: map PartMC -> MOSAIC
    nbin_a = aero_state_total_particles(aero_state)
    if (nbin_a > naerbin) then
       call DeallocateMemory()
       naerbin = nbin_a
       call AllocateMemory()
    end if
    aer = 0d0    ! initialize to zero
    ! work backwards for consistency with mosaic_to_partmc(), which
    ! has specific ordering requirements
    do i_part = aero_state_n_part(aero_state),1,-1
       num_conc = aero_weight_array_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
       do i_spec = 1,aero_data_n_spec(aero_data)
          i_spec_mosaic = aero_data%mosaic_index(i_spec)
          if (i_spec_mosaic > 0) then
             ! convert m^3(species) to nmol(species)/m^3(air)
             aer(i_spec_mosaic, 3, i_part) &   ! nmol/m^3(air)
                  = aero_state%apa%particle(i_part)%vol(i_spec) &
                  * conv_fac(i_spec) * num_conc
          end if
       end do
       ! handle water specially
       ! convert m^3(water) to kg(water)/m^3(air)
       water_a(i_part) = &
            aero_state%apa%particle(i_part)%vol(aero_data%i_water) &
            * aero_data%density(aero_data%i_water) * num_conc
       num_a(i_part) = 1d-6 * num_conc ! num conc (#/cc(air))
       jhyst_leg(i_part) = aero_state%apa%particle(i_part)%water_hyst_leg
    end do

    ! gas chemistry: map PartMC -> MOSAIC
    cnn = 0d0
    do i_spec = 1,gas_data_n_spec(gas_data)
       i_spec_mosaic = gas_data%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          ! convert ppbv to molec/cc
          cnn(i_spec_mosaic) = gas_state%mix_rat(i_spec) * cair_mlc / ppb
       end if
    end do
#ifdef PMC_USE_WRF
    particle_error = .false.
#endif
#endif

  end subroutine mosaic_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data MOSAIC -> PartMC.
  subroutine mosaic_to_partmc(env_state, aero_data, aero_state, gas_data, &
       gas_state, uuid)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: nbin_a, aer, num_a, jhyst_leg, &
#ifdef PMC_USE_WRF
         particle_error, &
#endif
         jtotal, water_a

    use module_data_mosaic_main, only: tbeg_sec, tcur_sec, tmid_sec, &
         dt_sec, dt_min, dt_aeroptic_min, RH, te, pr_atm, cnn, cair_mlc, &
         cair_molm3, ppb, avogad, msolar, cos_sza
#endif

    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> UUID for this simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid

#ifdef PMC_USE_MOSAIC
    ! local variables
    real(kind=dp) :: conv_fac(aero_data_n_spec(aero_data)), dum_var, num_conc
    integer :: i_part, i_spec, i_spec_mosaic
    real(kind=dp), allocatable :: reweight_num_conc(:)

#ifdef PMC_USE_WRF
    if (any(particle_error) .eqv. .true.) then
       call output_state('before_error', OUTPUT_TYPE_DIST, aero_data, &
            aero_state, gas_data, gas_state, env_state, &
            int(env_state%elapsed_time / dt_sec), tcur_sec, dt_sec, 1, &
            .false., .false., uuid)
    end if
#endif

    ! compute aerosol conversion factors
    do i_spec = 1,aero_data_n_spec(aero_data)
       ! converts m^3(species) to nmol(species)/m^3(air)
       conv_fac(i_spec) = 1d9 * aero_data%density(i_spec) &
            / aero_data%molec_weight(i_spec)
    enddo

    ! environmental parameters: map MOSAIC -> PartMC
    env_state%rel_humid = RH / 100d0
    env_state%temp = te
    env_state%pressure = pr_atm * const%air_std_press
    if (msolar == 1) then
       env_state%solar_zenith_angle = acos(cos_sza)
    end if
    cair_mlc = avogad*pr_atm/(82.056d0*te)   ! air conc [molec/cc]
    cair_molm3 = 1d6*pr_atm/(82.056d0*te)    ! air conc [mol/m^3]
    ppb = 1d9

    ! gas chemistry: map MOSAIC -> PartMC
    do i_spec = 1,gas_data_n_spec(gas_data)
       i_spec_mosaic = gas_data%mosaic_index(i_spec)
       if (i_spec_mosaic > 0) then
          ! convert molec/cc to ppbv
          gas_state%mix_rat(i_spec) = cnn(i_spec_mosaic) / cair_mlc * ppb
       end if
    end do

    ! We're modifying particle diameters, so the bin sort is now invalid
    aero_state%valid_sort = .false.

    ! aerosol data: map MOSAIC -> PartMC
    allocate(reweight_num_conc(aero_state_n_part(aero_state)))
    call aero_state_num_conc_for_reweight(aero_state, aero_data, &
         reweight_num_conc)
    do i_part = 1,aero_state_n_part(aero_state),1
       num_conc = aero_weight_array_num_conc(aero_state%awa, &
            aero_state%apa%particle(i_part), aero_data)
       do i_spec = 1,aero_data_n_spec(aero_data)
          i_spec_mosaic = aero_data%mosaic_index(i_spec)
          if (i_spec_mosaic > 0) then
             aero_state%apa%particle(i_part)%vol(i_spec) = &
                  ! convert nmol(species)/m^3(air) to m^3(species)
                  aer(i_spec_mosaic, 3, i_part) / (conv_fac(i_spec) * num_conc)
          end if
       end do
       aero_state%apa%particle(i_part)%water_hyst_leg = jhyst_leg(i_part)
       ! handle water specially
       ! convert kg(water)/m^3(air) to m^3(water)
       aero_state%apa%particle(i_part)%vol(aero_data%i_water) = &
            water_a(i_part) / aero_data%density(aero_data%i_water) / num_conc
    end do

#ifdef PMC_USE_WRF
    if (any(particle_error) .eqv. .true.) then
       call output_state('after_error', OUTPUT_TYPE_DIST, aero_data, &
            aero_state, gas_data, gas_state, env_state, &
            int(env_state%elapsed_time / dt_sec), tcur_sec, dt_sec, 1,&
            .false., .false., uuid)
       do i_part = aero_state_n_part(aero_state),1,-1
          if (particle_error(i_part)) then
             ! Move last element to the current spot if it is removed.
             if (i_part < aero_state_n_part(aero_state)) then
                reweight_num_conc(i_part) = reweight_num_conc( &
                     aero_state_n_part(aero_state))
             end if
             call aero_state_remove_particle_no_info(aero_state, i_part)
          end if
       end do
       ! Shrink the reweight_num_conc array to the number of particles
       reweight_num_conc = reweight_num_conc(1:aero_state_n_part(aero_state))
    end if
#endif

    ! adjust particles to account for weight changes
    call aero_state_reweight(aero_state, aero_data, reweight_num_conc)

#endif

  end subroutine mosaic_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do one timestep with MOSAIC.
  !!
  !! We currently also compute aerosol optical properties within this
  !! subroutine. In principle this could be done at data analysis
  !! time, rather than inside the timestepper. It's not clear if this
  !! really matters, however. Because of this mosaic_aero_optical() is
  !! currently disabled.
  subroutine mosaic_timestep(env_state, aero_data, aero_state, gas_data, &
       gas_state, do_optical, uuid)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_main, only: msolar
#endif

    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Whether to compute optical properties.
    logical, intent(in) :: do_optical
    !> UUID for this simulation.
    character(len=PMC_UUID_LEN), intent(in) :: uuid

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine SolarZenithAngle()
       end subroutine SolarZenithAngle
       subroutine IntegrateChemistry()
       end subroutine IntegrateChemistry
#ifdef PMC_USE_MOSAIC_MULTI_OPT
       subroutine aerosol_optical(i_wavelength)
         integer, optional :: i_wavelength
       end subroutine aerosol_optical
#else
       subroutine aerosol_optical()
       end subroutine aerosol_optical
#endif
    end interface

    ! map PartMC -> MOSAIC
    call mosaic_from_partmc(env_state, aero_data, aero_state, gas_data, &
         gas_state)

    if (msolar == 1) then
      call SolarZenithAngle
    end if

    call IntegrateChemistry

    ! map MOSAIC -> PartMC
    if (do_optical) then
       ! must do optical properties first, as mosaic_to_partmc() may
       ! change the number of particles
       call aerosol_optical
       call mosaic_aero_optical(env_state, aero_data, &
            aero_state, gas_data, gas_state)
    end if

    call mosaic_to_partmc(env_state, aero_data, aero_state, gas_data, &
         gas_state, uuid)
#endif

  end subroutine mosaic_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the optical properties of each aerosol particle.
  !> FIXME: currently disabled.
  !!
  !! At the moment we are computing the aerosol optical properties
  !! every timestep from withing mosaic_timestep. This decision should
  !! be re-evaluated at some point in the future.
  subroutine mosaic_aero_optical(env_state, aero_data, &
       aero_state, gas_data, gas_state)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: ri_shell_a, ri_core_a, &
         ext_cross, scat_cross, asym_particle, dp_core_a
#endif

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
#ifdef PMC_USE_MOSAIC_MULTI_OPT
       subroutine aerosol_optical(i_wavelength)
         integer, optional :: i_wavelength
       end subroutine aerosol_optical
#else
       subroutine aerosol_optical()
       end subroutine aerosol_optical
#endif
    end interface

    integer :: i_part

    ! map PartMC -> MOSAIC
!    call mosaic_from_partmc(env_state, aero_data, aero_state, &
!         gas_data, gas_state)

!    call aerosol_optical

    ! map MOSAIC -> PartMC
    ! work backwards for consistency with mosaic_to_partmc(), which
    ! has specific ordering requirements
    do i_part = aero_state_n_part(aero_state),1,-1
#ifdef PMC_USE_MOSAIC_MULTI_OPT
       aero_state%apa%particle(i_part)%absorb_cross_sect = &
            (ext_cross(i_part,:) - scat_cross(i_part,:)) / 1d4 ! (m^2)
       aero_state%apa%particle(i_part)%scatter_cross_sect = &
            scat_cross(i_part,:) / 1d4 ! (m^2)
       aero_state%apa%particle(i_part)%asymmetry = &
            asym_particle(i_part,:) ! (1)
       aero_state%apa%particle(i_part)%refract_shell = &
            cmplx(ri_shell_a(i_part,:), kind=dc) ! (1)
       aero_state%apa%particle(i_part)%refract_core =&
            cmplx(ri_core_a(i_part,:), kind=dc) ! (1)
       aero_state%apa%particle(i_part)%core_vol = &
            aero_data_diam2vol(aero_data, dp_core_a(i_part)) / 1d6 ! (m^3)
#else
       aero_state%apa%particle(i_part)%absorb_cross_sect = (ext_cross(i_part) &
            - scat_cross(i_part)) / 1d4                       ! (m^2)
       aero_state%apa%particle(i_part)%scatter_cross_sect = &
            scat_cross(i_part) / 1d4 ! (m^2)
       aero_state%apa%particle(i_part)%asymmetry = asym_particle(i_part) ! (1)
       aero_state%apa%particle(i_part)%refract_shell = &
            cmplx(ri_shell_a(i_part), kind=dc) ! (1)
       aero_state%apa%particle(i_part)%refract_core =&
            cmplx(ri_core_a(i_part), kind=dc) ! (1)
       aero_state%apa%particle(i_part)%core_vol = &
            aero_data_diam2vol(aero_data, dp_core_a(i_part)) / 1d6 ! (m^3)
#endif
    end do
#endif

  end subroutine mosaic_aero_optical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the optical properties of each aerosol particle for initial
  !> timestep.
  subroutine mosaic_aero_optical_init(env_state, aero_data, &
       aero_state, gas_data, gas_state)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_aero, only: ri_shell_a, ri_core_a, &
         ext_cross, scat_cross, asym_particle, dp_core_a
#endif

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
       subroutine load_mosaic_parameters()
       end subroutine load_mosaic_parameters
#ifdef PMC_USE_MOSAIC_MULTI_OPT
       subroutine aerosol_optical(i_wavelength)
         integer, optional :: i_wavelength
       end subroutine aerosol_optical
#else
       subroutine aerosol_optical()
       end subroutine aerosol_optical
#endif
    end interface

    call load_mosaic_parameters

    ! map PartMC -> MOSAIC
    call mosaic_from_partmc(env_state, aero_data, aero_state, &
         gas_data, gas_state)

    call aerosol_optical

    call mosaic_aero_optical(env_state, aero_data, &
         aero_state, gas_data, gas_state)
#endif

  end subroutine mosaic_aero_optical_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the optical properties of each aerosol particle for a single
  !> wavelength selected by index when MOSAIC has multiple wavelengths enabled.
  subroutine mosaic_aero_optical_single_wavelength(env_state, aero_data, &
       aero_state, gas_data, gas_state, i_wavelength)

#ifdef PMC_USE_MOSAIC
    use module_data_mosaic_main, only: RH, pr_atm, te
    use module_data_mosaic_aero, only: ri_shell_a, ri_core_a, &
         ext_cross, scat_cross, asym_particle, dp_core_a, &
         p_atm, RH_pc, aH2O, T_K
#endif

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Wavelength index to compute properties.
    integer, intent(in) :: i_wavelength

#ifdef PMC_USE_MOSAIC
    ! MOSAIC function interfaces
    interface
#ifdef PMC_USE_MOSAIC_MULTI_OPT
       subroutine aerosol_optical(i_wavelength)
         integer, optional :: i_wavelength
       end subroutine aerosol_optical
#else
       subroutine aerosol_optical()
       end subroutine aerosol_optical
#endif
       subroutine load_mosaic_parameters()
       end subroutine load_mosaic_parameters
    end interface

    integer :: i_part

    call load_mosaic_parameters

    ! map PartMC -> MOSAIC
    call mosaic_from_partmc(env_state, aero_data, aero_state, &
         gas_data, gas_state)

    RH_pc = RH                                ! RH(%)
    aH2O = 0.01*RH_pc                         ! aH2O (aerosol water activity)
    P_atm = pr_atm                            ! P(atm)
    T_K = te                                  ! T(K)
#ifdef PMC_USE_MOSAIC_MULTI_OPT
    call aerosol_optical(i_wavelength)
#endif
    call mosaic_aero_optical(env_state, aero_data, &
         aero_state, gas_data, gas_state)
#endif

  end subroutine mosaic_aero_optical_single_wavelength

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fetch the wavelengths that optical cross-sections are calculated at in
  !> MOSAIC. If not using multiple wavelengths, the wavelength is the default
  !> value of 550 nm.
  subroutine mosaic_optical_wavelengths(aero_data)

#ifdef PMC_USE_MOSAIC
#ifdef PMC_USE_MOSAIC_MULTI_OPT
    use module_data_mosaic_aero, only: wavelength
#endif
#endif
    type(aero_data_t), intent(inout) :: aero_data
    integer :: i

    if(allocated(aero_data%wavelengths)) deallocate(aero_data%wavelengths)
    allocate(aero_data%wavelengths(n_swbands))
    aero_data%wavelengths = 0.0d0
#ifdef PMC_USE_MOSAIC_MULTI_OPT
    do i = 1,n_swbands
       aero_data%wavelengths(i) = wavelength(i) * 1d-9
    end do
#else
    aero_data%wavelengths = 550.0d0 * 1d-9
#endif

  end subroutine mosaic_optical_wavelengths

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mosaic
