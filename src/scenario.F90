! Copyright (C) 2005-2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_scenario module.

!> The scenario_t structure and associated subroutines.
module pmc_scenario

  use pmc_gas_state
  use pmc_aero_dist
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_spec_file
  use pmc_aero_data
  use pmc_gas_data
  use pmc_chamber
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_QUADPACK
  use quadpack_double, only: dqagse
#endif

  !> Type code for an undefined or invalid loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_INVALID  = 0
  !> Type code for a zero loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_NONE     = 1
  !> Type code for a constant loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_CONSTANT = 2
  !> Type code for a loss rate function proportional to particle volume.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_VOLUME   = 3
  !> Type code for a loss rate function based on dry deposition
  integer, parameter :: SCENARIO_LOSS_FUNCTION_DRYDEP  = 4
  !> Type code for a loss rate function for chamber experiments.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_CHAMBER  = 5

  !> Parameter to switch between algorithms for particle loss.
  !! A value of 0 will always use the naive algorithm, and
  !! a value of 1 will always use the accept-reject algorithm.
  real(kind=dp), parameter :: SCENARIO_LOSS_ALG_THRESHOLD = 1.0d0

  !> Parameters for simulating dry deposition
  type drydep_params_t
     !> Reference height
     real(kind=dp) :: z_ref = 20.0d0
     !> Mean wind speed at reference height
     real(kind=dp) :: u_mean = 5.0d0
     !> Roughness length (land-use category dependent)
     ! crops, mixed farming (LUC 7 from Zhang et al., 2001)
     real(kind=dp) :: z_rough = 0.8d0
     !> Characteristic radius of collectors (land-use category dependent)
     real(kind=dp) :: A = 2.0d0 / 1000.0d0
     !> Impaction parameter
     real(kind=dp) :: alpha = 1.0d0
     !> Surface resistance parameter
     real(kind=dp) :: eps_0 = 3.0d0
     !> Diffusivity paramter
     real(kind=dp) :: gamma = 0.56d0 ! LUC-dependent for Zhang et al., 2001
     !> Brownian diffusion coefficient
     real(kind=dp) :: C_B = 1.0d0
     !> Interception coefficient
     real(kind=dp) :: C_IN = 0.5d0
     !> Impaction coefficient
     real(kind=dp) :: C_IM = 1.0d0
     !> Interception exponent
     real(kind=dp) :: nu = 2.0d0
     !> Impaction exponent
     real(kind=dp) :: beta = 2.0d0
  end type drydep_params_t

  !> Scenario data.
  !!
  !! This is everything needed to drive the scenario being simulated.
  !!
  !! The temperature, pressure, emissions and background states are profiles
  !! prescribed as functions of time by giving a number of times and
  !! the corresponding data. Simple data such as temperature and pressure is
  !! linearly interpolated between times, with constant interpolation
  !! outside of the range of times. Gases and aerosols are
  !! interpolated with gas_state_interp_1d() and
  !! aero_dist_interp_1d(), respectively.
  type scenario_t
     !> Temperature set-point times (s).
     real(kind=dp), allocatable :: temp_time(:)
     !> Temperatures at set-points (K).
     real(kind=dp), allocatable :: temp(:)

     !> Pressure set-point times (s).
     real(kind=dp), allocatable :: pressure_time(:)
     !> Pressures at set-points (Pa).
     real(kind=dp), allocatable :: pressure(:)

     !> Height set-point times (s).
     real(kind=dp), allocatable :: height_time(:)
     !> Heights at set-points (m).
     real(kind=dp), allocatable :: height(:)

     !> Gas emission set-point times (s).
     real(kind=dp), allocatable :: gas_emission_time(:)
     !> Gas emisssion rate scales at set-points (1).
     real(kind=dp), allocatable :: gas_emission_rate_scale(:)
     !> Gas emission rates at set-points (mol m^{-2} s^{-1}).
     type(gas_state_t), allocatable :: gas_emission(:)

     !> Gas-background dilution set-point times (s).
     real(kind=dp), allocatable :: gas_dilution_time(:)
     !> Gas-background dilution rates at set-points (s^{-1}).
     real(kind=dp), allocatable :: gas_dilution_rate(:)
     !> Background gas mixing ratios at set-points (ppb).
     type(gas_state_t), allocatable :: gas_background(:)

     !> Aerosol emission set-points times (s).
     real(kind=dp), allocatable :: aero_emission_time(:)
     !> Aerosol emission rate scales at set-points (1).
     real(kind=dp), allocatable :: aero_emission_rate_scale(:)
     !> Aerosol emissions at set-points (# m^{-2} s^{-1}).
     type(aero_dist_t), allocatable :: aero_emission(:)

     !> Aerosol-background dilution set-point times (s).
     real(kind=dp), allocatable :: aero_dilution_time(:)
     !> Aerosol-background dilution rates at set-points (s^{-1}).
     real(kind=dp), allocatable :: aero_dilution_rate(:)
     !> Aerosol background at set-points (# m^{-3}).
     type(aero_dist_t), allocatable :: aero_background(:)

     !> Type of loss rate function.
     integer :: loss_function_type
     !> Dry deposition parameters
     type(drydep_params_t) :: drydep
     !> Chamber parameters for wall loss and sedimentation.
     type(chamber_t) :: chamber

     !> Gas dry deposition set-point times (s); not allocated without gas
     !> dry deposition.
     real(kind=dp), allocatable :: gas_deposition_time(:)
     !> Gas dry deposition velocity scales at set-points (1).
     real(kind=dp), allocatable :: gas_deposition_rate_scale(:)
     !> Gas dry deposition velocities at set-points (m s^{-1}), stored in
     !> the mixing ratio array of a gas state.
     type(gas_state_t), allocatable :: gas_deposition(:)

     !> Whether air leaving the mixed layer is kept in a residual layer
     !> above it and entrained again when the mixed layer grows.
     logical :: do_residual_layer = .false.
  end type scenario_t

  !> Residual layer above the mixed layer.
  !!
  !! When the mixed layer becomes shallower, the air between the new and
  !! the old mixing height is added to the residual layer, which spans
  !! from the mixing height to \c top. When the mixed layer grows, it
  !! entrains residual-layer air up to \c top and background air above
  !! it. With MOSAIC the residual layer has the same chemistry as the
  !! mixed layer, at the temperature, pressure and humidity of the mixed
  !! layer; it has no emissions, coagulation, dilution or deposition (see
  !! scenario_update_residual_layer()).
  type residual_layer_t
     !> Height of the top of the residual layer (m). The residual layer is
     !> empty when \c top is not above the mixing height.
     real(kind=dp) :: top = 0d0
     !> Gas mixing ratios in the residual layer (ppb).
     type(gas_state_t) :: gas_state
     !> Particles of the residual layer (concentrations at the temperature
     !> and pressure of the mixed layer).
     type(aero_state_t) :: aero_state
  end type residual_layer_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the time-dependent contents of the
  !> environment. Thereafter scenario_update_env_state() should be used.
  subroutine scenario_init_env_state(scenario, env_state, time)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time

    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    env_state%pressure = interp_1d(scenario%pressure_time, &
         scenario%pressure, time)
    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time
    ! FIXME: should compute this at some point
    env_state%solar_zenith_angle = 0d0

  end subroutine scenario_init_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update time-dependent contents of the environment.
  !> scenario_init_env_state() should have been called at the start.
  subroutine scenario_update_env_state(scenario, env_state, time)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time

    !> Ambient water vapor pressure (Pa).
    real(kind=dp) :: pmv_old, pmv_new
    !> Ambient pressure (Pa)
    real(kind=dp) :: pressure_old
    !> Ambient temperature (K)
    real(kind=dp) :: temp_old

    ! Update temperature and pressure and adjust relative humidity to maintain
    ! water mixing ratio.

    pmv_old = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    pressure_old = env_state%pressure
    temp_old = env_state%temp

    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    env_state%pressure = interp_1d(scenario%pressure_time, &
         scenario%pressure, time)

    pmv_new = pmv_old * env_state%pressure / pressure_old
    env_state%rel_humid = pmv_new / env_state_sat_vapor_pressure(env_state)

    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time

  end subroutine scenario_update_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do gas emissions and background dilution.
  !!
  !! Emissions are given as an areal rate \f$e(t)\f$, and then divided
  !! by the current box height \f$h(t)\f$ to obtain a volume
  !! rate. There is also a dimensionless rate scaling \f$r(t)\f$. All
  !! input functions are asusumed constant over the timestep, so the
  !! concentration \f$c(t)\f$ change is given by
  !! \f[
  !!     c(t) = c(0) + \frac{r t}{h} e.
  !! \f]
  !!
  !! We model dilution by considering a gas concentration \f$c(t)\f$
  !! in a box of height \f$h(t)\f$, subject to first-order dilution
  !! with a rate \f$\lambda(t)\f$. Then the effective dilution rate is
  !! \f[
  !!     \lambda_{\rm eff}(t) = \lambda(t) + \frac{\dot{h}(t)}{h(t)}
  !! \f]
  !! and the evolution of \f$c(t)\f$ is given by
  !! \f[
  !!     \dot{c}(t) = - \lambda_{\rm eff}(t) c(t).
  !! \f]
  !! Solving this with separation of variables gives
  !! \f[
  !!     \frac{c(t)}{c(0)} = \frac{h(0)}{h(t)}
  !!                         \exp\left( - \int_0^t \lambda(t)\,dt\right).
  !! \f]
  !! If we define \f$p = c(t)/c(0)\f$ to be the remaining proportion
  !! of the initial concentration, and \f$b\f$ to be the constant
  !! background concentration, then we have
  !! \f[
  !!     c(t) = p(t) c(0) + (1 - p(t)) b.
  !! \f]
  !! We assume constant \f$\lambda\f$ and we only do entrainment with
  !! increasing height \f$h(t)\f$, so we have
  !! \f[
  !!     p(t) = \min\left(1, \frac{h(0)}{h(t)}\right) \exp(-\lambda t).
  !! \f]
  !! With a residual layer (\c residual_layer in the spec file) the
  !! factor \f$h(0)/h(t)\f$ is left out here, and the entrainment is done
  !! by scenario_update_residual_layer() instead.
  !!
  !! With gas dry deposition (\c gas_deposition in the spec file), each
  !! gas species \f$i\f$ is then lost at the first-order rate
  !! \f$r(t) v_i(t) / h(t)\f$, with \f$v_i\f$ the dry deposition
  !! velocity and \f$r\f$ the scaling of the deposition profile:
  !! \f[
  !!     c_i(t) = c_i(0) \exp(- r v_i t / h).
  !! \f]
  subroutine scenario_update_gas_state(scenario, delta_t, env_state, &
       old_env_state, gas_data, gas_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state to update.
    type(gas_state_t), intent(inout) :: gas_state

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    real(kind=dp) :: deposition_rate_scale
    type(gas_state_t) :: emissions, background, deposition

    ! emissions
    if (size(scenario%gas_emission) > 0) then
       call gas_state_interp_1d(scenario%gas_emission, &
            scenario%gas_emission_time, scenario%gas_emission_rate_scale, &
            env_state%elapsed_time, emissions, emission_rate_scale)
       call gas_state_molar_conc_to_ppb(emissions, env_state)
       p = emission_rate_scale * delta_t / env_state%height
       call gas_state_add_scaled(gas_state, emissions, p)
    end if
#ifndef PMC_USE_WRF
    ! dilution; with a residual layer, the entrainment when the mixed
    ! layer grows is done by scenario_update_residual_layer()
    call gas_state_interp_1d(scenario%gas_background, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         env_state%elapsed_time, background, dilution_rate)
    p = exp(- dilution_rate * delta_t)
    if ((env_state%height > old_env_state%height) &
         .and. (.not. scenario%do_residual_layer)) then
       p = p * old_env_state%height / env_state%height
    end if
    call gas_state_scale(gas_state, p)
    call gas_state_add_scaled(gas_state, background, 1d0 - p)
    ! dry deposition: first-order loss at the rate v_d / h
    if (allocated(scenario%gas_deposition_time)) then
       call gas_state_interp_1d(scenario%gas_deposition, &
            scenario%gas_deposition_time, &
            scenario%gas_deposition_rate_scale, env_state%elapsed_time, &
            deposition, deposition_rate_scale)
       gas_state%mix_rat = gas_state%mix_rat &
            * exp(- deposition_rate_scale * deposition%mix_rat * delta_t &
            / env_state%height)
    end if
#endif
    call gas_state_ensure_nonnegative(gas_state)

  end subroutine scenario_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution for a particle aerosol
  !> distribution.
  !!
  !! See scenario_update_gas_state() for a description of the
  !! model. We additionally scale the number concentration to account
  !! for temperature changes.
  subroutine scenario_update_aero_state(scenario, delta_t, env_state, &
       old_env_state, aero_data, aero_state, n_emit, n_dil_in, n_dil_out, &
       allow_doubling, allow_halving)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero state to update.
    type(aero_state_t), intent(inout) :: aero_state
    !> Number of emitted particles.
    integer, intent(out) :: n_emit
    !> Number of diluted-in particles.
    integer, intent(out) :: n_dil_in
    !> Number of diluted-out particles.
    integer, intent(out) :: n_dil_out
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    real(kind=dp), parameter :: sample_timescale = 3600.0d0
    real(kind=dp) :: characteristic_factor, sample_timescale_effective
    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(aero_dist_t) :: emissions, background
    type(aero_state_t) :: aero_state_delta

    sample_timescale_effective = max(1.0, min(sample_timescale, &
         env_state%elapsed_time))
    characteristic_factor = sample_timescale_effective / delta_t

    ! emissions
    if (size(scenario%aero_emission) > 0) then
       call aero_dist_interp_1d(scenario%aero_emission, &
            scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
            env_state%elapsed_time, emissions, emission_rate_scale)
       p = emission_rate_scale * delta_t / env_state%height
       call aero_state_add_aero_dist_sample(aero_state, aero_data, &
            emissions, p, characteristic_factor, env_state%elapsed_time, &
            allow_doubling, allow_halving, n_emit)
    end if
#ifndef PMC_USE_WRF
    ! dilution
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         env_state%elapsed_time, background, dilution_rate)
    p = exp(- dilution_rate * delta_t)
    if ((env_state%height > old_env_state%height) &
         .and. (.not. scenario%do_residual_layer)) then
       p = p * old_env_state%height / env_state%height
    end if
    ! loss to background
    call aero_state_zero(aero_state_delta)
    call aero_state_copy_weight(aero_state, aero_state_delta)
    call aero_state_sample_particles(aero_state, aero_state_delta, &
         aero_data, 1d0 - p, AERO_INFO_DILUTION)
    n_dil_out = aero_state_total_particles(aero_state_delta)
    ! addition from background
    call aero_state_add_aero_dist_sample(aero_state, aero_data, &
         background, 1d0 - p, characteristic_factor, env_state%elapsed_time, &
         allow_doubling, allow_halving, n_dil_in)

    ! particle loss function
    call scenario_particle_loss(scenario, delta_t, aero_data, aero_state, &
         env_state)
#endif

#ifdef PMC_USE_WRF
    call aero_weight_array_scale(aero_state%awa, &
         old_env_state%inverse_density * (1.0d0 / env_state%inverse_density))
#else
    ! update computational volume
    call aero_weight_array_scale(aero_state%awa, &
         old_env_state%temp * env_state%pressure &
         / (env_state%temp * old_env_state%pressure))
#endif

  end subroutine scenario_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize an empty residual layer above the mixed layer.
  subroutine residual_layer_init(residual_layer, env_state, gas_data, &
       aero_state)

    !> Residual layer.
    type(residual_layer_t), intent(inout) :: residual_layer
    !> Environment state (gives the initial mixing height).
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol state of the mixed layer (gives the particle weighting).
    type(aero_state_t), intent(in) :: aero_state

    residual_layer%top = env_state%height
    call gas_state_set_size(residual_layer%gas_state, &
         gas_data_n_spec(gas_data))
    residual_layer%gas_state%mix_rat = 0d0
    call aero_state_zero(residual_layer%aero_state)
    call aero_state_copy_weight(aero_state, residual_layer%aero_state)

  end subroutine residual_layer_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add to \c aero_state_to copies of the particles of \c
  !> aero_state_from so that the expected added concentration is \c
  !> sample_prob times the concentration of \c aero_state_from, for any
  !> weighting of the two states.
  !!
  !! Each particle of \c aero_state_from is copied
  !! <tt>prob_round(sample_prob * c_from / c_to)</tt> times, where \c
  !! c_from and \c c_to are the number concentrations it represents with
  !! the weighting of \c aero_state_from and of \c aero_state_to, so the
  !! result is unbiased particle by particle. (aero_state_sample_particles()
  !! counts removals and is unbiased only when this ratio is the same for
  !! all particles.) Weight classes without particles in \c aero_state_to
  !! first take the weights of \c aero_state_from, so that the ratio is
  !! never formed with the unset weights of an unused class. The copies get
  !! new IDs.
  subroutine scenario_add_sample_copies(aero_state_from, aero_state_to, &
       aero_data, sample_prob)

    !> State to copy from.
    type(aero_state_t), intent(in) :: aero_state_from
    !> State to add to.
    type(aero_state_t), intent(inout) :: aero_state_to
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Fraction of the concentration of \c aero_state_from to add.
    real(kind=dp), intent(in) :: sample_prob

    integer :: i_part, i_copy, n_copies, i_group, i_class, n_in_class
    real(kind=dp) :: num_conc_from, num_conc_to
    type(aero_particle_t) :: new_aero_particle

    do i_class = 1,aero_weight_array_n_class(aero_state_to%awa)
       n_in_class = 0
       do i_group = 1,aero_weight_array_n_group(aero_state_to%awa)
          n_in_class = n_in_class &
               + aero_state_total_particles(aero_state_to, i_group, i_class)
       end do
       if (n_in_class == 0) then
          aero_state_to%awa%weight(:, i_class) &
               = aero_state_from%awa%weight(:, i_class)
       end if
    end do

    do i_part = 1,aero_state_n_part(aero_state_from)
       num_conc_from = aero_weight_array_num_conc(aero_state_from%awa, &
            aero_state_from%apa%particle(i_part), aero_data)
       num_conc_to = aero_weight_array_num_conc(aero_state_to%awa, &
            aero_state_from%apa%particle(i_part), aero_data)
       n_copies = prob_round(sample_prob * num_conc_from / num_conc_to)
       do i_copy = 1,n_copies
          new_aero_particle = aero_state_from%apa%particle(i_part)
          call aero_particle_new_id(new_aero_particle)
          call aero_state_add_particle(aero_state_to, new_aero_particle, &
               aero_data)
       end do
    end do

  end subroutine scenario_add_sample_copies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange of the mixed layer with the residual layer above it when
  !> the mixing height changes from \f$h_0\f$ to \f$h_1\f$.
  !!
  !! The residual layer spans from the mixing height to its top
  !! \f$z\f$, so its depth before the step is \f$d = \max(z - h_0,
  !! 0)\f$. When the mixing height falls, the layer between \f$h_1\f$
  !! and \f$h_0\f$ leaves the mixed layer, whose concentrations do not
  !! change, and joins the residual layer:
  !! \f[
  !!     c_{\rm rl} \leftarrow (1 - f) c_{\rm rl} + f c, \quad
  !!     f = \frac{h_0 - h_1}{d + h_0 - h_1}, \quad
  !!     z \leftarrow \max(z, h_0).
  !! \f]
  !! When the mixing height rises, the mixed layer entrains \f$d_{\rm
  !! rl} = \max(\min(h_1, z) - h_0, 0)\f$ of residual-layer air and
  !! \f$d_{\rm bg} = h_1 - h_0 - d_{\rm rl}\f$ of background air:
  !! \f[
  !!     c \leftarrow \frac{h_0 c + d_{\rm rl} c_{\rm rl}
  !!     + d_{\rm bg} b}{h_1},
  !! \f]
  !! while the concentrations of the remaining residual layer do not
  !! change; it is used up once \f$h_1 \ge z\f$. Particles leave a layer
  !! by random removal and enter the other layer as copies with the
  !! expected concentration for any particle weighting (see
  !! scenario_add_sample_copies()), so the particles of the residual
  !! layer return with their composition and mixing state. The residual
  !! layer has no emissions, coagulation, dilution or deposition; with
  !! MOSAIC its chemistry is computed in run_part_timestep() at the
  !! temperature, pressure and humidity of the mixed layer. Its
  !! concentrations follow the temperature and pressure of the mixed
  !! layer, so that its mixing ratios are conserved. Exchange with the
  !! background at the dilution rate of the background profiles applies
  !! to the mixed layer only (scenario_update_gas_state() and
  !! scenario_update_aero_state()).
  subroutine scenario_update_residual_layer(scenario, delta_t, env_state, &
       old_env_state, gas_state, aero_data, aero_state, residual_layer, &
       allow_doubling, allow_halving, n_dil_in, n_dil_out)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Gas state of the mixed layer.
    type(gas_state_t), intent(inout) :: gas_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state of the mixed layer.
    type(aero_state_t), intent(inout) :: aero_state
    !> Residual layer.
    type(residual_layer_t), intent(inout) :: residual_layer
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving
    !> Number of particles entrained into the mixed layer.
    integer, intent(out) :: n_dil_in
    !> Number of particles removed from the mixed layer by dilution.
    integer, intent(out) :: n_dil_out

    real(kind=dp), parameter :: sample_timescale = 3600.0d0
    real(kind=dp) :: h_old, h_new, d_rl, f, d_from_rl, d_from_bg
    real(kind=dp) :: dilution_rate, characteristic_factor
    type(gas_state_t) :: gas_background
    type(aero_dist_t) :: aero_background
    type(aero_state_t) :: aero_state_delta
    integer :: n_before, n_bg

    n_dil_in = 0
    n_dil_out = 0
    h_old = old_env_state%height
    h_new = env_state%height

    ! the residual layer follows the temperature and pressure of the
    ! mixed layer, as the mixed layer does in scenario_update_aero_state()
    call aero_weight_array_scale(residual_layer%aero_state%awa, &
         old_env_state%temp * env_state%pressure &
         / (env_state%temp * old_env_state%pressure))

    if (h_new < h_old) then
       d_rl = max(residual_layer%top - h_old, 0d0)
       f = (h_old - h_new) / (d_rl + h_old - h_new)
       if (d_rl == 0d0) then
          ! empty residual layer: start with the weighting of the mixed layer
          call aero_state_zero(residual_layer%aero_state)
          call aero_state_copy_weight(aero_state, residual_layer%aero_state)
       end if
       residual_layer%gas_state%mix_rat = (1d0 - f) &
            * residual_layer%gas_state%mix_rat + f * gas_state%mix_rat
       ! remove the fraction f of the residual-layer particles ...
       call aero_state_zero(aero_state_delta)
       call aero_state_copy_weight(residual_layer%aero_state, &
            aero_state_delta)
       call aero_state_sample_particles(residual_layer%aero_state, &
            aero_state_delta, aero_data, f, AERO_INFO_NONE)
       ! ... and add copies of the fraction f of the mixed layer
       call scenario_add_sample_copies(aero_state, &
            residual_layer%aero_state, aero_data, f)
       residual_layer%top = max(residual_layer%top, h_old)
    elseif (h_new > h_old) then
       d_from_rl = max(min(h_new, residual_layer%top) - h_old, 0d0)
       d_from_bg = max(h_new - h_old - d_from_rl, 0d0)
       call gas_state_interp_1d(scenario%gas_background, &
            scenario%gas_dilution_time, scenario%gas_dilution_rate, &
            env_state%elapsed_time, gas_background, dilution_rate)
       gas_state%mix_rat = (h_old * gas_state%mix_rat &
            + d_from_rl * residual_layer%gas_state%mix_rat &
            + d_from_bg * gas_background%mix_rat) / h_new
       ! dilution of the mixed-layer particles ...
       call aero_state_zero(aero_state_delta)
       call aero_state_copy_weight(aero_state, aero_state_delta)
       call aero_state_sample_particles(aero_state, aero_state_delta, &
            aero_data, 1d0 - h_old / h_new, AERO_INFO_DILUTION)
       n_dil_out = aero_state_total_particles(aero_state_delta)
       ! ... particles of the residual layer ...
       if (d_from_rl > 0d0) then
          n_before = aero_state_total_particles(aero_state)
          call scenario_add_sample_copies(residual_layer%aero_state, &
               aero_state, aero_data, d_from_rl / h_new)
          n_dil_in = aero_state_total_particles(aero_state) - n_before
       end if
       ! ... and background particles above the residual layer
       if (d_from_bg > 0d0) then
          call aero_dist_interp_1d(scenario%aero_background, &
               scenario%aero_dilution_time, scenario%aero_dilution_rate, &
               env_state%elapsed_time, aero_background, dilution_rate)
          characteristic_factor = max(1d0, min(sample_timescale, &
               env_state%elapsed_time)) / delta_t
          call aero_state_add_aero_dist_sample(aero_state, aero_data, &
               aero_background, d_from_bg / h_new, characteristic_factor, &
               env_state%elapsed_time, allow_doubling, allow_halving, n_bg)
          n_dil_in = n_dil_in + n_bg
       end if
       if (h_new >= residual_layer%top) then
          residual_layer%top = h_new
          residual_layer%gas_state%mix_rat = 0d0
          call aero_state_zero(residual_layer%aero_state)
       end if
    end if

  end subroutine scenario_update_residual_layer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> binned aerosol distribution.
  !!
  !! See scenario_update_gas_state() for a description of the model.
  subroutine scenario_update_aero_binned(scenario, delta_t, env_state, &
       old_env_state, bin_grid, aero_data, do_aero_dilution, aero_binned)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to do dilution.
    logical, intent(in) :: do_aero_dilution
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(aero_dist_t) :: emissions, background
    type(aero_binned_t) :: emissions_binned, background_binned

    ! emissions
    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
         env_state%elapsed_time, emissions, emission_rate_scale)
    call aero_binned_add_aero_dist(emissions_binned, bin_grid, aero_data, &
         emissions)
    p = emission_rate_scale * delta_t / env_state%height
    call aero_binned_add_scaled(aero_binned, emissions_binned, p)

    ! dilution
    if (do_aero_dilution) then
       call aero_dist_interp_1d(scenario%aero_background, &
            scenario%aero_dilution_time, scenario%aero_dilution_rate, &
            env_state%elapsed_time, background, dilution_rate)
       call aero_binned_add_aero_dist(background_binned, bin_grid, aero_data, &
            background)
       p = exp(- dilution_rate * delta_t)
       if (env_state%height > old_env_state%height) then
          p = p * old_env_state%height / env_state%height
       end if
       call aero_binned_scale(aero_binned, p)
       call aero_binned_add_scaled(aero_binned, background_binned, 1d0 - p)
    end if

    ! loss
    call scenario_binned_loss(scenario, bin_grid, delta_t, aero_data, &
         env_state, aero_binned)

  end subroutine scenario_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the modal aerosol distribution to account for particle loss.
  subroutine scenario_update_aero_modes(aero_dist, del_t, env_state, &
       density, scenario)

    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Timestep.
    real(kind=dp), intent(in) :: del_t
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Particle density (kg m^-3), assumed uniform across all modes.
    real(kind=dp), intent(in) :: density
    !> Scenario data.
    type(scenario_t), intent(in) :: scenario

    real(kind=dp) :: N, d_pg, ln_sigma_g
    real(kind=dp) :: m_0_rate, m_3_rate
    real(kind=dp) :: M, new_M
    real(kind=dp) :: new_N, new_d_pg
    integer :: i_mode

    if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_INVALID) then
       return
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_NONE) then
       return
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_DRYDEP) then
       do i_mode = 1,aero_dist_n_mode(aero_dist)
          N = aero_dist%mode(i_mode)%num_conc

          if (N == 0d0) cycle

          d_pg = aero_dist%mode(i_mode)%char_radius * 2.0d0
          ln_sigma_g = aero_dist%mode(i_mode)%log10_std_dev_radius &
               / log10(exp(1.0d0))

          ! Integrated deposition rate for the 0-th moment (num. conc.)
          m_0_rate = -1.0d0 * scenario_integrated_loss_rate_drydep(scenario, &
               aero_dist%mode(i_mode), 0.0d0, density, env_state)
          new_N = N * exp(m_0_rate * del_t)

          aero_dist%mode(i_mode)%num_conc = new_N

          ! Integrated deposition rate for the 3-rd moment
          ! (proportional to volume conc.)
          m_3_rate = -1.0d0 * scenario_integrated_loss_rate_drydep(scenario, &
               aero_dist%mode(i_mode), 3.0d0, density, env_state)
          M = N * d_pg**3.0d0 * exp((3.0d0**2.0d0) / 2.0d0 &
               * (ln_sigma_g**2.0d0))
          new_M = M * exp(m_3_rate * del_t)

          ! New geometric mean diameter
          new_d_pg = (new_M / new_N * exp(-(3.0d0**2.0d0) / 2.0d0 &
               * (ln_sigma_g**2.0d0)))**(1.0d0/3.0d0)
          aero_dist%mode(i_mode)%char_radius = new_d_pg / 2.0d0
       end do
    else
       call die_msg(192500240, "loss_function must be none or drydep " &
            // "for modal runs; unsupported loss function id: " &
            // trim(integer_to_string(scenario%loss_function_type)))
    end if

  end subroutine scenario_update_aero_modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evaluate a loss rate function.
  real(kind=dp) function scenario_loss_rate(scenario, vol, density, &
       aero_data, env_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Volume of particle (m^3).
    real(kind=dp), intent(in) :: vol
    !> Density of particle (kg/m^3).
    real(kind=dp), intent(in) :: density
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    scenario_loss_rate = 0d0
    if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_INVALID) then
       scenario_loss_rate = 0d0
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_NONE) then
       scenario_loss_rate = 0d0
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_CONSTANT) &
         then
       scenario_loss_rate = 1d-3
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_VOLUME) then
       scenario_loss_rate = 1d15*vol
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_DRYDEP) &
         then
       scenario_loss_rate = scenario_loss_rate_drydep(vol, density, &
            aero_data, env_state, scenario)
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_CHAMBER) &
         then
       scenario_loss_rate = chamber_loss_rate_wall(scenario%chamber, vol, &
            aero_data, env_state) &
            + chamber_loss_rate_sedi(scenario%chamber, vol, density, &
            aero_data, env_state)
    else
       call die_msg(201594391, "Unknown loss function id: " &
            // trim(integer_to_string(scenario%loss_function_type)))
    end if

  end function scenario_loss_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute and return the dry deposition rate for a given particle.
  !! All equations used here are written in detail in the file
  !! \c doc/deposition/deposition.tex.
  real(kind=dp) function scenario_loss_rate_drydep(vol, density, aero_data, &
       env_state, scenario)

    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Particle density (kg m^-3).
    real(kind=dp), intent(in) :: density
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Scenario data.
    type(scenario_t), intent(in) :: scenario

    real(kind=dp) :: V_d
    real(kind=dp) :: V_s
    real(kind=dp) :: d_p
    real(kind=dp) :: density_air
    real(kind=dp) :: visc_d, visc_k
    real(kind=dp) :: gas_speed, gas_mean_free_path
    real(kind=dp) :: knud, cunning
    real(kind=dp) :: grav
    real(kind=dp) :: R_s, R_a
    real(kind=dp) :: diff_p
    real(kind=dp) :: von_karman
    real(kind=dp) :: St, Sc, u_star
    real(kind=dp) :: E_B, E_IM, E_IN, R1
    type(drydep_params_t) :: drydep_params

    drydep_params = scenario%drydep

    ! particle diameter
    d_p = aero_data_vol2diam(aero_data, vol)
    ! density of air
    density_air = (const%air_molec_weight * env_state%pressure) &
         / (const%univ_gas_const * env_state%temp)
    ! dynamic viscosity
    visc_d = 1.8325d-5 * (416.16 / (env_state%temp + 120.0d0)) &
         * (env_state%temp / 296.16)**1.5d0
    ! kinematic viscosity
    visc_k = visc_d / density_air
    ! gas speed
    gas_speed = sqrt((8.0d0 * const%boltzmann * env_state%temp &
         * const%avagadro) / (const%pi * const%air_molec_weight))
    ! gas free path
    gas_mean_free_path = (2.0d0 * visc_d) / (density_air * gas_speed)
    ! knudson number
    knud = (2.0d0 * gas_mean_free_path) / d_p
    ! cunningham correction factor
    cunning = 1.0d0 + knud * (1.257d0 + 0.4d0 * exp(-1.1d0 / knud))
    ! gravity
    grav = 9.81d0
    ! Compute V_s
    V_s = (density * d_p**2.0d0 * grav * cunning) / (18.0d0 * visc_d)

    ! Aerodynamic resistance
    ! For neutral stability
    u_star = .4d0 * drydep_params%u_mean &
         / log(drydep_params%z_ref / drydep_params%z_rough)
    R_a = (1.0d0 / (.4d0 * u_star)) &
         * log(drydep_params%z_ref / drydep_params%z_rough)
    ! Brownian diffusion efficiency
    diff_p = (const%boltzmann * env_state%temp * cunning) &
         / (3.d0 * const%pi * visc_d * d_p)
    Sc = visc_k / diff_p
    E_B = drydep_params%C_B * Sc**(-drydep_params%gamma)

    ! Interception efficiency
    ! Characteristic radius of large collectors
    E_IN = drydep_params%C_IN * (d_p / drydep_params%A)**drydep_params%nu

    ! Impaction efficiency
    St = (V_s * u_star) / (grav * drydep_params%A)
    E_IM = drydep_params%C_IM * (St &
         / (drydep_params%alpha + St))**drydep_params%beta

    ! Rebound correction
    R1 = exp(-St**.5d0)

    ! Surface resistance
    R_s = 1.0d0 / (drydep_params%eps_0 * u_star * (E_B + E_IN + E_IM) * R1)

    ! Dry deposition
    V_d = V_s + (1.0d0 / (R_a + R_s))

    ! The loss rate
    scenario_loss_rate_drydep = V_d / env_state%height

  end function scenario_loss_rate_drydep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute and return the integrated dry deposition rate for a given
  !> lognormal aerosol mode (modal approximation).
  real(kind=dp) function scenario_integrated_loss_rate_drydep(scenario, &
       aero_mode, moment, density, env_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol mode.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Moment to calculate loss rate for.
    real(kind=dp), intent(in) :: moment
    !> Particle density assumed for entire mode (kg m^-3).
    real(kind=dp), intent(in) :: density
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: V_d_hat, V_g_hat
    real(kind=dp) :: V_g_bar
    real(kind=dp) :: d_pg, ln_sigma_g
    real(kind=dp) :: density_air
    real(kind=dp) :: visc_d, visc_k
    real(kind=dp) :: gas_speed, gas_mean_free_path
    real(kind=dp) :: knud
    real(kind=dp) :: D_bar, D_hat
    real(kind=dp) :: St, Sc
    real(kind=dp) :: u_star
    real(kind=dp) :: R_a, R_s
    real(kind=dp) :: E_B, E_IN, E_IM, R1
    type(drydep_params_t) :: drydep_params

#ifdef PMC_USE_QUADPACK
     ! Do numerical integration when QUADPACK is available
     scenario_integrated_loss_rate_drydep = &
         scenario_integrated_loss_rate_drydep_quadpack( &
            scenario, aero_mode, moment, density, env_state)
     return
#endif

    drydep_params = scenario%drydep

    ! particle diameter equal to geometric mean diameter
    d_pg = aero_mode%char_radius * 2.0d0
    ! natural log of geometric standard deviation
    ln_sigma_g = aero_mode%log10_std_dev_radius / log10(exp(1.0d0))
    ! density of air
    density_air = (const%air_molec_weight * env_state%pressure) &
         / (const%univ_gas_const * env_state%temp)
    ! dynamic viscosity
    visc_d = 1.8325d-5 * (416.16 / (env_state%temp + 120.0d0)) &
         * (env_state%temp / 296.16)**1.5d0
    ! kinematic viscosity
    visc_k = visc_d / density_air
    ! gas speed
    gas_speed = sqrt((8.0d0 * const%boltzmann * env_state%temp &
         * const%avagadro) / (const%pi * const%air_molec_weight))
    ! gas mean free path
    gas_mean_free_path = (2.0d0 * visc_d) / (density_air * gas_speed)
    ! Knudsen number
    knud = (2.0d0 * gas_mean_free_path) / d_pg
    ! Settling velocity
    V_g_bar = (density * d_pg**2.0d0 * const%std_grav) / (18.0d0 * visc_d)
    ! Compute integrated settling velocity
    V_g_hat = V_g_bar &
         * (exp((4.0d0 * moment + 4.0d0) / 2.0d0 * ln_sigma_g**2.0d0) &
         + 1.246d0 * knud &
         * exp((2.0d0 * moment + 1.0d0) / 2.0d0 * ln_sigma_g**2.0d0))

    ! Aerodynamic resistance (assuming neutral stability)
    u_star = 0.4d0 * drydep_params%u_mean / log(drydep_params%z_ref &
         / drydep_params%z_rough)
    R_a = (1.0d0 / (0.4d0 * u_star)) &
         * log(drydep_params%z_ref / drydep_params%z_rough)

    ! Brownian diffusivity
    D_bar = (const%boltzmann * env_state%temp) &
         / (3.0d0 * const%pi * visc_d * d_pg)
    ! Compute integrated Brownian diffusivity
    D_hat = D_bar &
         * ((exp((-2.0d0 * moment + 1.0d0) / 2.0d0 * ln_sigma_g**2.0d0) &
         + 1.246d0 * knud &
         * exp((-4.0d0 * moment + 4.0d0) / 2.0d0 * ln_sigma_g**2.0d0)))
    ! Schmidt number based on integrated diffusivity
    Sc = visc_k / D_hat
    ! Collection efficiency due to Brownian diffusion
    E_B = drydep_params%C_B * Sc**(-drydep_params%gamma)

    ! Collection efficiency due to interception
    E_IN = drydep_params%C_IN * (d_pg / drydep_params%A)**drydep_params%nu

    ! Stokes number based on integrated settling velocity
    St = (V_g_hat * u_star) / (const%std_grav * drydep_params%A)
    ! Collection efficiency due to impaction
    E_IM = drydep_params%C_IM * (St / (drydep_params%alpha + St)) &
         **drydep_params%beta

    ! Rebound correction
    R1 = exp(-St**0.5d0)

    ! Surface resistance
    R_s = 1.0d0 / (drydep_params%eps_0 * u_star * (E_B + E_IN + E_IM) * R1)

    ! Integrated deposition velocity
    V_d_hat = V_g_hat + (1.0d0 / (R_a + R_s))

    ! Loss rate
    scenario_integrated_loss_rate_drydep = V_d_hat / env_state%height

  end function scenario_integrated_loss_rate_drydep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_QUADPACK
  real(kind=dp) function scenario_integrated_loss_rate_drydep_quadpack( &
       scenario, aero_mode, moment, density, env_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol mode.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Moment to calculate loss rate for.
    real(kind=dp), intent(in) :: moment
    !> Particle density assumed for entire mode (kg m^-3).
    real(kind=dp), intent(in) :: density
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: d_pg, ln_sigma_g
    real(kind=dp) :: density_air
    real(kind=dp) :: visc_d, visc_k
    real(kind=dp) :: gas_speed, gas_mean_free_path
    real(kind=dp) :: M_k
    real(kind=dp) :: lower, upper, epsabs, epsrel, result, abserr, ref
    integer :: limit, neval, ier, last
    real(kind=dp), allocatable :: alist(:), blist(:), rlist(:), elist(:)
    integer, allocatable :: iord(:)

    type(drydep_params_t) :: drydep_params

    drydep_params = scenario%drydep

    ! particle diameter equal to geometric mean diameter
    d_pg = aero_mode%char_radius * 2.0d0
    ! natural log of geometric standard deviation
    ln_sigma_g = aero_mode%log10_std_dev_radius / log10(exp(1.0d0))
    ! density of air
    density_air = (const%air_molec_weight * env_state%pressure) &
         / (const%univ_gas_const * env_state%temp)
    ! dynamic viscosity
    visc_d = 1.8325d-5 * (416.16 / (env_state%temp + 120.0d0)) &
         * (env_state%temp / 296.16)**1.5d0
    ! kinematic viscosity
    visc_k = visc_d / density_air
    ! gas speed
    gas_speed = sqrt((8.0d0 * const%boltzmann * env_state%temp &
         * const%avagadro) / (const%pi * const%air_molec_weight))
    ! gas mean free path
    gas_mean_free_path = (2.0d0 * visc_d) / (density_air * gas_speed)

    ref = exp(log(d_pg) + moment*ln_sigma_g**2)
    lower = ref / (exp(ln_sigma_g)*10.0)
    upper = ref * exp(ln_sigma_g)*10.0

    ! Set integration parameters
    limit = 1000
    epsabs = 1.0d-10
    epsrel = 1.0d-6

    ! Allocate arrays for QUADPACK integration
    allocate(alist(limit), blist(limit), rlist(limit), &
         elist(limit), iord(limit))

    ! Call QUADPACK integration routine
    call dqagse(dep_vel_integrand, lower, upper, epsabs, epsrel, limit, &
         result, abserr, neval, ier, alist, blist, rlist, elist, iord, last)

    call assert_msg(909106718, ier == 0, &
       "QUADPACK integration failed, error code: " &
       // trim(integer_to_string(ier)))

    M_k = d_pg**moment * exp(moment**2 * ln_sigma_g**2 / 2.0d0)

    ! Integration result
    scenario_integrated_loss_rate_drydep_quadpack = 1.0d0 / M_k * result &
         / env_state%height

    ! Clean up
    deallocate(alist, blist, rlist, elist, iord)

  contains

  real(kind=dp) function dep_vel_integrand(d_p)
    real(kind=dp), intent(in) :: d_p

    ! Local variables
    real(kind=dp) :: V_d, V_s
    real(kind=dp) :: knud_local, cunning
    real(kind=dp) :: diff_p, Sc, St
    real(kind=dp) :: u_star, R_a, R_s
    real(kind=dp) :: E_B, E_IN, E_IM, R1
    real(kind=dp) :: ln_dp, ln_dp_g, n_ddp

    ! Knudsen number for this diameter
    knud_local = (2.0d0 * gas_mean_free_path) / d_p

    ! Cunningham correction factor
    cunning = 1.0d0 + knud_local * (1.257d0 + 0.4d0 * exp(-1.1d0 / knud_local))

    ! Settling velocity
    V_s = (density * d_p**2.0d0 * const%std_grav * cunning) / (18.0d0 * visc_d)

    ! Friction velocity and aerodynamic resistance
    u_star = 0.4d0 * drydep_params%u_mean &
         / log(drydep_params%z_ref / drydep_params%z_rough)
    R_a = (1.0d0 / (0.4d0 * u_star)) &
         * log(drydep_params%z_ref / drydep_params%z_rough)

    ! Particle diffusivity and Schmidt number
    diff_p = (const%boltzmann * env_state%temp * cunning) &
         / (3.0d0 * const%pi * visc_d * d_p)
    Sc = visc_k / diff_p

    ! Collection efficiencies
    E_B = drydep_params%C_B * Sc**(-drydep_params%gamma)
    E_IN = drydep_params%C_IN * (d_p / drydep_params%A)**drydep_params%nu

    ! Stokes number and impaction efficiency
    St = (V_s * u_star) / (const%std_grav * drydep_params%A)
    E_IM = drydep_params%C_IM &
         * (St / (drydep_params%alpha + St))**drydep_params%beta

    ! Rebound correction
    R1 = exp(-St**0.5d0)

    ! Surface resistance
    R_s = 1.0d0 / (drydep_params%eps_0 * u_star * (E_B + E_IN + E_IM) * R1)

    ! Deposition velocity
    V_d = V_s + (1.0d0 / (R_a + R_s))

    ! Log-normal size distribution
    ln_dp = log(d_p)
    ln_dp_g = log(d_pg)
    n_ddp = (1.0d0/(sqrt(2.0d0 * const%pi) * d_p * ln_sigma_g)) &
         * exp(-((ln_dp - ln_dp_g)**2) / (2.0d0 * ln_sigma_g**2))

    ! Final integrand
    dep_vel_integrand = d_p**moment * V_d * n_ddp
  end function dep_vel_integrand
  end function scenario_integrated_loss_rate_drydep_quadpack
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Updates an array to contain the integrated deposition velocities for the
  !> given moment of each aerosol mode in the distribution.
  subroutine scenario_modal_drydep_velocities(scenario, aero_dist, moment, &
       density, env_state, velocities)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Moment to compute deposition velocities for.
    real(kind=dp), intent(in) :: moment
    !> Density for the mode.
    real(kind=dp), intent(in) :: density
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Modal dry deposition velocities.
    real(kind=dp), intent(inout) :: velocities(:)

    integer :: i_mode

    do i_mode = 1,size(velocities)
       velocities(i_mode) = scenario_integrated_loss_rate_drydep( &
            scenario, aero_dist%mode(i_mode), moment, density, env_state) &
            * env_state%height
    end do

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute and return the max loss rate function for a given volume.
  real(kind=dp) function scenario_loss_rate_max(scenario, vol, aero_data, &
       env_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Particle volume (m^3).
    real(kind=dp), intent(in) :: vol
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    !> Number of density sample points.
    integer, parameter :: n_sample = 3

    real(kind=dp) :: d, d_min, d_max, loss
    integer :: i

    d_min = minval(aero_data%density)
    d_max = maxval(aero_data%density)

    scenario_loss_rate_max = 0d0
    do i = 1,n_sample
       d = interp_linear_disc(d_min, d_max, n_sample, i)
       loss = scenario_loss_rate(scenario, vol, d, aero_data, env_state)
       scenario_loss_rate_max = max(scenario_loss_rate_max, loss)
    end do

  end function scenario_loss_rate_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute an upper bound on the maximum kernel value for each
  !> bin.
  !! Value over_scale is multiplied to the maximum sampled value
  !! to get the upper bound.  A tighter bound may be reached if over_scale
  !! is smaller, but that also risks falling below a kernel value.
  subroutine scenario_loss_rate_bin_max(scenario, bin_grid, aero_data, &
       env_state, loss_max)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Maximum loss vals.
    real(kind=dp), intent(out) :: loss_max(bin_grid_size(bin_grid))

    !> Number of sample points per bin.
    integer, parameter :: n_sample = 3
    !> Over-estimation scale factor parameter.
    real(kind=dp), parameter :: over_scale = 2d0

    real(kind=dp) :: v_low, v_high, vol, r, r_max
    integer :: b, i

    do b = 1,bin_grid_size(bin_grid)
       v_low = aero_data_rad2vol(aero_data, bin_grid%edges(b))
       v_high = aero_data_rad2vol(aero_data, bin_grid%edges(b + 1))
       r_max = 0d0
       do i = 1,n_sample
          vol = interp_linear_disc(v_low, v_high, n_sample, i)
          r = scenario_loss_rate_max(scenario, vol, aero_data, env_state)
          r_max = max(r_max, r)
       end do
       loss_max(b) = r_max * over_scale
    end do

  end subroutine scenario_loss_rate_bin_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Performs stochastic particle loss for one time-step.
  !! If a particle \c i_part has a scenario_loss_rate() value of rate, then the
  !! probability p will be removed by this function is
  !! <tt>1 - exp(-delta_t*rate)</tt>.
  !! Uses an accept-reject algorithm for efficiency, in which a particle
  !! is first sampled with rate <tt>1 - exp(-delta_t*over_rate) </tt>
  !! and then accepted with rate
  !! <tt>(1 - exp(-delta_t*rate))/(1 - exp(-delta_t*over_rate))</tt>.
  subroutine scenario_particle_loss(scenario, delta_t, aero_data, aero_state, &
       env_state)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    integer :: c, b, s, i_part
    real(kind=dp) :: over_rate, over_prob, rand_real, rand_geom

    if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_NONE .or. &
        scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_INVALID) return

    if (SCENARIO_LOSS_ALG_THRESHOLD <= 0d0) then
       ! use naive algorithm for everything
       do i_part = aero_state%apa%n_part, 1, -1
          call scenario_try_single_particle_loss(scenario, delta_t, &
               aero_data, aero_state, env_state, i_part, 1d0)
       end do
       return
    end if

    call aero_state_sort(aero_state, aero_data)

    if (.not. aero_state%aero_sorted%removal_rate_bounds_valid) then
       call scenario_loss_rate_bin_max(scenario, &
            aero_state%aero_sorted%bin_grid, aero_data, env_state, &
            aero_state%aero_sorted%removal_rate_max)
       aero_state%aero_sorted%removal_rate_bounds_valid = .true.
    end if

    do c = 1,aero_sorted_n_class(aero_state%aero_sorted)
       do b = 1,aero_sorted_n_bin(aero_state%aero_sorted)
          over_rate = aero_state%aero_sorted%removal_rate_max(b)
          if (delta_t * over_rate <= 0d0) cycle
          over_prob = 1d0 - exp(-delta_t * over_rate)
          if (over_prob >= SCENARIO_LOSS_ALG_THRESHOLD) then
             ! use naive algorithm over bin
             do s = aero_state%aero_sorted%size_class%inverse(b, c)%n_entry, &
                  1,-1
                i_part = &
                     aero_state%aero_sorted%size_class%inverse(b, c)%entry(s)
                call scenario_try_single_particle_loss(scenario, delta_t, &
                     aero_data, aero_state, env_state, i_part, 1d0)
             end do
          else
             ! use accept-reject algorithm over bin
             s = aero_state%aero_sorted%size_class%inverse(b, c)%n_entry + 1
             do while (.true.)
                rand_real = pmc_random()
                if (rand_real <= 0d0) exit
                rand_geom = -log(rand_real) / (delta_t * over_rate) + 1d0
                if (rand_geom >= real(s, kind=dp)) exit
                s = s - floor(rand_geom)

                ! note: floor(rand_geom) is a random geometric variable
                ! with accept probability 1 - exp(-delta_t*over_rate)

                i_part = &
                     aero_state%aero_sorted%size_class%inverse(b, c)%entry(s)
                call scenario_try_single_particle_loss(scenario, delta_t, &
                     aero_data, aero_state, env_state, i_part, over_prob)
             end do
          end if
       end do
    end do

    !call aero_state_check_sort(aero_state)

  end subroutine scenario_particle_loss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test a candidate particle to see if it should be removed,
  !> and remove if necessary.
  !! Particle is removed with probability
  !! (1d0 - exp(-delta_t*rate))/over_prob, where rate is the loss function
  !! evaluated for the given particle.
  subroutine scenario_try_single_particle_loss(scenario, delta_t, &
       aero_data, aero_state, env_state, i_part, over_prob)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Index of particle to attempt removal
    integer, intent(in) :: i_part
    !> Overestimated removal probability used previously
    real(kind=dp), intent(in) :: over_prob

    real(kind=dp) :: prob, rate, vol, density
    type(aero_info_t) :: aero_info

    vol = aero_particle_volume(aero_state%apa%particle(i_part))
    density = aero_particle_density(aero_state%apa%particle(i_part), aero_data)
    rate = scenario_loss_rate(scenario, vol, density, aero_data, env_state)
    prob = 1d0 - exp(-delta_t * rate)
    call warn_assert_msg(295846288, prob <= over_prob, &
         "particle loss upper bound estimation is too tight: " &
         // trim(real_to_string(prob)) // " > " &
         // trim(real_to_string(over_prob)) )
    if (pmc_random() * over_prob > prob) return

    aero_info%id = aero_state%apa%particle(i_part)%id
    aero_info%action = AERO_INFO_DILUTION
    aero_info%other_id = 0
    call aero_state_remove_particle_with_info(aero_state, i_part, aero_info)

  end subroutine scenario_try_single_particle_loss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Performs loss for a binned aerosol distribution for one time step.
  subroutine scenario_binned_loss(scenario, bin_grid, delta_t, aero_data, &
       env_state, aero_binned)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environmental state.
    type(env_state_t), intent(in) :: env_state
    !> Binned aerosol data.
    type(aero_binned_t), intent(inout) :: aero_binned

    integer :: i_bin
    real(kind=dp) :: density, vol, loss_rate, p

    if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_NONE .or. &
         scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_INVALID) then
       return
    else if (scenario%loss_function_type == SCENARIO_LOSS_FUNCTION_DRYDEP) then

       ! Assumes a single aerosol species, as run_sect() requires. Both the
       ! density and the number diagnosis below use species 1 only, and
       ! must be generalized if more species are allowed.
       call assert_msg(258732274, aero_data_n_spec(aero_data) == 1, &
            "sectional loss assumes a single aerosol species")
       density = aero_data%density(1)

       do i_bin = 1,bin_grid_size(bin_grid)

          if (aero_binned%num_conc(i_bin) <= 0d0) cycle

          vol = aero_data_rad2vol(aero_data, bin_grid%centers(i_bin))
          loss_rate = scenario_loss_rate(scenario, vol, density, &
               aero_data, env_state)
          p = exp(- loss_rate * delta_t)
          aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) * p
          aero_binned%num_conc(i_bin) = aero_binned%vol_conc(i_bin,1) / vol

       end do
    else
       call die_msg(769879198, "loss_function must be none or drydep " &
            // "for sectional runs; unsupported loss function id: " &
            // trim(integer_to_string(scenario%loss_function_type)))
    end if

  end subroutine scenario_binned_loss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether any of the contained aerosol modes are of the given type.
  elemental logical function scenario_contains_aero_mode_type(scenario, &
       aero_mode_type)

    !> Scenario data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol mode type to test for.
    integer, intent(in) :: aero_mode_type

    scenario_contains_aero_mode_type &
         = any(aero_dist_contains_aero_mode_type(scenario%aero_emission, &
         aero_mode_type)) &
         .or. any(aero_dist_contains_aero_mode_type(scenario%aero_background, &
         aero_mode_type))

  end function scenario_contains_aero_mode_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment data from an spec file.
  subroutine spec_file_read_scenario(file, gas_data, aero_data, &
       read_aero_weight_classes, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Whether the weight classes for each source are specified in inputs.
    logical, intent(in) :: read_aero_weight_classes
    !> Scenario data.
    type(scenario_t), intent(inout) :: scenario

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file
    character(len=SPEC_LINE_MAX_VAR_LEN) :: function_name
    type(spec_line_t) :: line
    logical :: eof
    integer :: i

    ! note that we have to hard-code the list for doxygen below

    !> \page input_format_scenario Input File Format: Scenario
    !!
    !! The scenario parameters are:
    !! <ul>
    !! <li> \b temp_profile (string): the name of the file from which to
    !!      read the temperature profile --- the file format should be
    !!      \subpage input_format_temp_profile
    !! <li> \b pressure_profile (string): the name of the file from which to
    !!      read the pressure profile --- the file format should be
    !!      \subpage input_format_pressure_profile
    !! <li> \b height_profile (string): the name of the file from which
    !!      to read the mixing layer height profile --- the file format
    !!      should be \subpage input_format_height_profile
    !! <li> \b gas_emissions (string): the name of the file from which to
    !!      read the gas emissions profile --- the file format should be
    !!      \subpage input_format_gas_profile
    !! <li> \b gas_background (string): the name of the file from which
    !!      to read the gas background profile --- the file format should
    !!      be \subpage input_format_gas_profile
    !! <li> \b aero_emissions (string): the name of the file from which
    !!      to read the aerosol emissions profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! <li> \b aero_background (string): the name of the file from which
    !!      to read the aerosol background profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! <li> \b loss_function (string): the type of loss function ---
    !!      must be one of: \c none for no particle loss, \c constant
    !!      for constant loss rate, \c volume for particle loss proportional
    !!      to particle volume, \c drydep for particle loss proportional
    !!      to dry deposition velocity, or \c chamber for a chamber model.
    !!      Sectional and modal runs support only \c none and \c drydep.
    !!      If \c loss_function is \c drydep, it may be followed by
    !!      \b drydep_params (string), the name of the file from which to
    !!      read the dry deposition parameters --- the file format should
    !!      be \subpage input_format_drydep_params. If \c drydep_params is
    !!      omitted, the default parameters are used and a warning is
    !!      printed.
    !!      If \c loss_function is \c chamber, then the following
    !!      parameters must also be provided:
    !!      - \subpage input_format_chamber
    !! <li> \b gas_deposition (string, optional): the name of the file
    !!      from which to read the gas dry deposition profile --- the file
    !!      format should be \subpage input_format_gas_profile, with dry
    !!      deposition velocities (m/s) as the species values and
    !!      dimensionless scalings as the rates. Each gas species is lost
    !!      at the rate (rate times velocity) / (mixing layer height). If
    !!      this line is omitted there is no gas dry deposition.
    !! <li> \b residual_layer (logical, optional): whether the air that
    !!      leaves the mixed layer when the mixing height falls is kept in
    !!      a residual layer above it and entrained again, with its
    !!      particles, when the mixing height rises; above the top of the
    !!      residual layer the mixed layer entrains background air
    !!      (particle-resolved runs only). If this line is omitted or \c
    !!      no, the mixed layer always entrains background air.
    !! </ul>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    ! temperature profile
    call spec_file_read_string(file, "temp_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "temp", &
         scenario%temp_time, scenario%temp)
    call spec_file_close(sub_file)

    ! pressure profile
    call spec_file_read_string(file, "pressure_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "pressure", &
         scenario%pressure_time, scenario%pressure)
    call spec_file_close(sub_file)

    ! height profile
    call spec_file_read_string(file, "height_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "height", &
         scenario%height_time, scenario%height)
    call spec_file_close(sub_file)

    ! gas emissions profile
    call spec_file_read_string(file, "gas_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_emission_time, scenario%gas_emission_rate_scale, &
         scenario%gas_emission)
    call spec_file_close(sub_file)

    ! gas background profile
    call spec_file_read_string(file, "gas_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         scenario%gas_background)
    call spec_file_close(sub_file)

    ! aerosol emissions profile
    call spec_file_read_string(file, "aero_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         read_aero_weight_classes, scenario%aero_emission_time, &
         scenario%aero_emission_rate_scale, scenario%aero_emission)
    call spec_file_close(sub_file)

    ! aerosol background profile
    call spec_file_read_string(file, "aero_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         read_aero_weight_classes, scenario%aero_dilution_time, &
         scenario%aero_dilution_rate, scenario%aero_background)
    call spec_file_close(sub_file)

    ! loss function
    call spec_file_read_string(file, 'loss_function', function_name)
    if (trim(function_name) == 'none') then
       scenario%loss_function_type = SCENARIO_LOSS_FUNCTION_NONE
    else if (trim(function_name) == 'constant') then
       scenario%loss_function_type = SCENARIO_LOSS_FUNCTION_CONSTANT
    else if (trim(function_name) == 'volume') then
       scenario%loss_function_type = SCENARIO_LOSS_FUNCTION_VOLUME
    else if (trim(function_name) == 'drydep') then
       scenario%loss_function_type = SCENARIO_LOSS_FUNCTION_DRYDEP
       call spec_file_read_line_no_eof(file, line)
       call spec_file_unread_line(file)
       if (line%name /= 'drydep_params') then
          call warn_msg(735291468, "using default dry deposition parameters")
       else
          call spec_file_read_string(file, 'drydep_params', sub_filename)
          call spec_file_open(sub_filename, sub_file)
          call spec_file_read_drydep_params(sub_file, scenario%drydep)
          call spec_file_close(sub_file)
       end if
    else if (trim(function_name) == 'chamber') then
       scenario%loss_function_type = SCENARIO_LOSS_FUNCTION_CHAMBER
       call spec_file_read_chamber(file, scenario%chamber)
    else
       call spec_file_die_msg(518248400, file, &
            "Unknown loss function type: " // trim(function_name))
    end if

    ! optional gas dry deposition profile
    call spec_file_read_line(file, line, eof)
    if (.not. eof) call spec_file_unread_line(file)
    if ((.not. eof) .and. (line%name == 'gas_deposition')) then
       call spec_file_read_string(file, 'gas_deposition', sub_filename)
       call spec_file_open(sub_filename, sub_file)
       call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
            scenario%gas_deposition_time, &
            scenario%gas_deposition_rate_scale, scenario%gas_deposition)
       call spec_file_close(sub_file)
       do i = 1,size(scenario%gas_deposition)
          call assert_msg(410857123, &
               all(scenario%gas_deposition(i)%mix_rat >= 0d0), &
               "negative gas dry deposition velocity in " &
               // trim(sub_filename))
       end do
       call spec_file_read_line(file, line, eof)
       if (.not. eof) call spec_file_unread_line(file)
    end if

    ! optional residual layer
    if ((.not. eof) .and. (line%name == 'residual_layer')) then
       call spec_file_read_logical(file, 'residual_layer', &
            scenario%do_residual_layer)
    end if

  end subroutine spec_file_read_scenario

  ! the following blocks belong in the subroutine above, but they are
  ! outside because Doxygen 1.8.7 doesn't resolve references when
  ! multiple \page blocks are in one subroutine

  !> \page input_format_temp_profile Input File Format: Temperature Profile
  !!
  !! A temperature profile input file must consist of two lines:
  !! - the first line must begin with \c time and should be followed
  !!   by \f$N\f$ space-separated real scalars, giving the times (in
  !!   s after the start of the simulation) of the temperature set
  !!   points --- the times must be in increasing order
  !! - the second line must begin with \c temp and should be followed
  !!   by \f$N\f$ space-separated real scalars, giving the
  !!   temperatures (in K) at the corresponding times
  !!
  !! The temperature profile is linearly interpolated between the
  !! specified times, while before the first time it takes the first
  !! temperature value and after the last time it takes the last
  !! temperature value.
  !!
  !! Example:
  !! <pre>
  !! time  0    600  1800  # time (in s) after simulation start
  !! temp  270  290  280   # temperature (in K)
  !! </pre>
  !! Here the temperature starts at 270&nbsp;K at the start of the
  !! simulation, rises to 290&nbsp;K after 10&nbsp;min, and then
  !! falls again to 280&nbsp;K at 30&nbsp;min. Between these times
  !! the temperature is linearly interpolated, while after
  !! 30&nbsp;min it is held constant at 280&nbsp;K.
  !!
  !! See also:
  !!   - \ref spec_file_format --- the input file text format
  !!   - \ref input_format_scenario --- the environment data
  !!     containing the temperature profile

  !> \page input_format_pressure_profile Input File Format: Pressure Profile
  !!
  !! A pressure profile input file must consist of two lines:
  !! - the first line must begin with \c time and should be followed
  !!   by \f$N\f$ space-separated real scalars, giving the times (in
  !!   s after the start of the simulation) of the pressure set
  !!   points --- the times must be in increasing order
  !! - the second line must begin with \c pressure and should be followed
  !!   by \f$N\f$ space-separated real scalars, giving the
  !!   pressures (in Pa) at the corresponding times
  !!
  !! The pressure profile is linearly interpolated between the
  !! specified times, while before the first time it takes the first
  !! pressure value and after the last time it takes the last
  !! pressure value.
  !!
  !! Example:
  !! <pre>
  !! time      0    600  1800  # time (in s) after simulation start
  !! pressure  1e5  9e4  7.5e4 # pressure (in Pa)
  !! </pre>
  !! Here the pressure starts at 1e5&nbsp;Pa at the start of the
  !! simulation, decreases to 9e4&nbsp;Pa after 10&nbsp;min, and then
  !! decreases further to 7.5e4&nbsp;Pa at 30&nbsp;min. Between these times
  !! the pressure is linearly interpolated, while after
  !! 30&nbsp;min it is held constant at 7.5e4&nbsp;Pa.
  !!
  !! See also:
  !!   - \ref spec_file_format --- the input file text format
  !!   - \ref input_format_scenario --- the environment data
  !!     containing the pressure profile

  !> \page input_format_height_profile Input File Format: Mixing Layer Height Profile
  !!
  !! A mixing layer height profile input file must consist of two
  !! lines:
  !! - the first line must begin with \c time and should be followed
  !!   by \f$N\f$ space-separated real scalars, giving the times (in
  !!   s after the start of the simulation) of the height set
  !!   points --- the times must be in increasing order
  !! - the second line must begin with \c height and should be
  !!   followed by \f$N\f$ space-separated real scalars, giving the
  !!   mixing layer heights (in m) at the corresponding times
  !!
  !! The mixing layer height profile is linearly interpolated
  !! between the specified times, while before the first time it
  !! takes the first height value and after the last time it takes
  !! the last height value.
  !!
  !! Example:
  !! <pre>
  !! time    0    600   1800  # time (in s) after simulation start
  !! height  500  1000  800   # mixing layer height (in m)
  !! </pre>
  !! Here the mixing layer height starts at 500&nbsp;m at the start
  !! of the simulation, rises to 1000&nbsp;m after 10&nbsp;min, and
  !! then falls again to 800&nbsp;m at 30&nbsp;min. Between these
  !! times the mixing layer height is linearly interpolated, while
  !! after 30&nbsp;min it is held constant at 800&nbsp;m.
  !!
  !! See also:
  !!   - \ref spec_file_format --- the input file text format
  !!   - \ref input_format_scenario --- the environment data
  !!     containing the mixing layer height profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read dry deposition parameters from a spec file.
  subroutine spec_file_read_drydep_params(file, drydep_params)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Dry deposition parameters.
    type(drydep_params_t), intent(inout) :: drydep_params

    !> \page input_format_drydep_params Input File Format: Dry Deposition Parameters
    !!
    !! Dry deposition is simulated using the specified parameters:
    !! - \b z_ref (real, unit m): the reference height \f$z_{\rm ref}\f$ used
    !!   in the calculation of aerodynamic resistance \f$R_a\f$
    !! - \b u_mean (real, unit m s^{-1}): the wind speed at the reference
    !!   height
    !! - \b z_rough (real, unit m): the roughness length associated with
    !!   the surface
    !! - \b A (real, unit m): the characteristic radius of collectors
    !!   associated with the surface
    !! - \b alpha (real, dimensionless): the parameter \f$\alpha\f$ used
    !!   in the calculation of impaction efficiency \f$E_{\rm IM}\f$
    !! - \b eps_0 (real, dimensionless): the empirical constant used in the
    !!   calculation of surface resistance \f$R_s\f$
    !! - \b gamma (real, dimensionless): the exponent \f$\gamma\f$ used in the
    !!   calculation of Brownian diffusion collection
    !!   efficiency \f$E_{\rm B}\f$
    !! - \b C_B (real, dimensionless): the coefficient for Brownian diffusion
    !!   collection efficiency \f$E_{\rm B}\f$
    !! - \b C_IN (real, dimensionless): the coefficient for interception
    !!   collection efficiency \f$E_{\rm IN}\f$
    !! - \b C_IM (real, dimensionless): the coefficient for impaction
    !!   collection efficiency \f$E_{\rm IM}\f$
    !! - \b nu (real, dimensionless): the exponent \f$\nu\f$ used in the
    !!   calculation of interception collection efficiency \f$E_{\rm IN}\f$
    !! - \b beta (real, dimensionless): the exponent \f$\beta\f$ used in the
    !!   calculation of impaction collection efficiency \f$E_{\rm IM}\f$
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the prescribed profiles of
    !!     other environment data

    call spec_file_read_real(file, "z_ref", drydep_params%z_ref)
    call spec_file_read_real(file, "u_mean", drydep_params%u_mean)
    call spec_file_read_real(file, "z_rough", drydep_params%z_rough)
    call spec_file_read_real(file, "A", drydep_params%A)
    call spec_file_read_real(file, "alpha", drydep_params%alpha)
    call spec_file_read_real(file, "eps_0", drydep_params%eps_0)
    call spec_file_read_real(file, "gamma", drydep_params%gamma)
    call spec_file_read_real(file, "C_B", drydep_params%C_B)
    call spec_file_read_real(file, "C_IN", drydep_params%C_IN)
    call spec_file_read_real(file, "C_IM", drydep_params%C_IM)
    call spec_file_read_real(file, "nu", drydep_params%nu)
    call spec_file_read_real(file, "beta", drydep_params%beta)

  end subroutine spec_file_read_drydep_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write dry deposition parameters to a NetCDF file.
  subroutine drydep_params_output_netcdf(drydep_params, ncid)

    !> Dry deposition parameters.
    type(drydep_params_t), intent(in) :: drydep_params
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    associate (d => drydep_params)
      call pmc_nc_write_real(ncid, d%z_ref,   "drydep_z_ref",   unit="m")
      call pmc_nc_write_real(ncid, d%u_mean,  "drydep_u_mean",  &
                             unit="m s^{-1}")
      call pmc_nc_write_real(ncid, d%z_rough, "drydep_z_rough", unit="m")
      call pmc_nc_write_real(ncid, d%A,       "drydep_A",       unit="m")
      call pmc_nc_write_real(ncid, d%alpha,   "drydep_alpha",   unit="1")
      call pmc_nc_write_real(ncid, d%eps_0,   "drydep_eps_0",   unit="1")
      call pmc_nc_write_real(ncid, d%gamma,   "drydep_gamma",   unit="1")
      call pmc_nc_write_real(ncid, d%C_B,     "drydep_C_B",     unit="1")
      call pmc_nc_write_real(ncid, d%C_IN,    "drydep_C_IN",    unit="1")
      call pmc_nc_write_real(ncid, d%C_IM,    "drydep_C_IM",    unit="1")
      call pmc_nc_write_real(ncid, d%nu,      "drydep_nu",      unit="1")
      call pmc_nc_write_real(ncid, d%beta,    "drydep_beta",    unit="1")
    end associate

  end subroutine drydep_params_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read dry deposition parameters from a NetCDF file.
  subroutine drydep_params_input_netcdf(drydep_params, ncid, must_be_present)

    !> Dry deposition parameters.
    type(drydep_params_t), intent(inout) :: drydep_params
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Whether the parameters must be present in the file (default
    !> \c .true.). They are only written when the loss function is dry
    !> deposition, so readers that may be given either kind of file
    !> should pass \c .false., which leaves \c drydep_params unchanged.
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    use_must_be_present = .true.
    if (present(must_be_present)) use_must_be_present = must_be_present

    if (.not. use_must_be_present) then
       status = nf90_inq_varid(ncid, "drydep_z_ref", varid)
       if (status == NF90_ENOTVAR) return
    end if

    associate (d => drydep_params)
      call pmc_nc_read_real(ncid, d%z_ref,   "drydep_z_ref")
      call pmc_nc_read_real(ncid, d%u_mean,  "drydep_u_mean")
      call pmc_nc_read_real(ncid, d%z_rough, "drydep_z_rough")
      call pmc_nc_read_real(ncid, d%A,       "drydep_A")
      call pmc_nc_read_real(ncid, d%alpha,   "drydep_alpha")
      call pmc_nc_read_real(ncid, d%eps_0,   "drydep_eps_0")
      call pmc_nc_read_real(ncid, d%gamma,   "drydep_gamma")
      call pmc_nc_read_real(ncid, d%C_B,     "drydep_C_B")
      call pmc_nc_read_real(ncid, d%C_IN,    "drydep_C_IN")
      call pmc_nc_read_real(ncid, d%C_IM,    "drydep_C_IM")
      call pmc_nc_read_real(ncid, d%nu,      "drydep_nu")
      call pmc_nc_read_real(ncid, d%beta,    "drydep_beta")
    end associate

  end subroutine drydep_params_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_scenario(val)

    !> Value to pack.
    type(scenario_t), intent(in) :: val

    integer :: total_size, i

    total_size = &
         pmc_mpi_pack_size_real_array(val%temp_time) &
         + pmc_mpi_pack_size_real_array(val%temp) &
         + pmc_mpi_pack_size_real_array(val%pressure_time) &
         + pmc_mpi_pack_size_real_array(val%pressure) &
         + pmc_mpi_pack_size_real_array(val%height_time) &
         + pmc_mpi_pack_size_real_array(val%height) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_time) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_rate_scale) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_time) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_rate_scale) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_rate) &
         + pmc_mpi_pack_size_integer(val%loss_function_type) &
         + pmc_mpi_pack_size_drydep(val%drydep) &
         + pmc_mpi_pack_size_chamber(val%chamber) &
         + pmc_mpi_pack_size_real_array(val%gas_deposition_time) &
         + pmc_mpi_pack_size_real_array(val%gas_deposition_rate_scale) &
         + pmc_mpi_pack_size_logical(val%do_residual_layer)
    if (allocated(val%gas_deposition_time)) then
       do i = 1,size(val%gas_deposition)
          total_size = total_size &
               + pmc_mpi_pack_size_gas_state(val%gas_deposition(i))
       end do
    end if
    if (allocated(val%gas_emission_time)) then
       do i = 1,size(val%gas_emission)
          total_size = total_size &
               + pmc_mpi_pack_size_gas_state(val%gas_emission(i))
       end do
    end if
    if (allocated(val%gas_dilution_time)) then
       do i = 1,size(val%gas_background)
          total_size = total_size &
               + pmc_mpi_pack_size_gas_state(val%gas_background(i))
       end do
    end if
    if (allocated(val%aero_emission_time)) then
       do i = 1,size(val%aero_emission)
          total_size = total_size &
               + pmc_mpi_pack_size_aero_dist(val%aero_emission(i))
       end do
    end if
    if (allocated(val%aero_dilution_time)) then
       do i = 1,size(val%aero_background)
          total_size = total_size &
               + pmc_mpi_pack_size_aero_dist(val%aero_background(i))
       end do
    end if

    pmc_mpi_pack_size_scenario = total_size

  end function pmc_mpi_pack_size_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_pack_real_array(buffer, position, val%temp)
    call pmc_mpi_pack_real_array(buffer, position, val%pressure_time)
    call pmc_mpi_pack_real_array(buffer, position, val%pressure)
    call pmc_mpi_pack_real_array(buffer, position, val%height_time)
    call pmc_mpi_pack_real_array(buffer, position, val%height)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_rate_scale)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, &
         val%aero_emission_rate_scale)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_rate)
    call pmc_mpi_pack_integer(buffer, position, val%loss_function_type)
    call pmc_mpi_pack_drydep(buffer, position, val%drydep)
    call pmc_mpi_pack_chamber(buffer, position, val%chamber)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_deposition_time)
    call pmc_mpi_pack_real_array(buffer, position, &
         val%gas_deposition_rate_scale)
    call pmc_mpi_pack_logical(buffer, position, val%do_residual_layer)
    if (allocated(val%gas_deposition_time)) then
       do i = 1,size(val%gas_deposition)
          call pmc_mpi_pack_gas_state(buffer, position, &
               val%gas_deposition(i))
       end do
    end if
    if (allocated(val%gas_emission_time)) then
       do i = 1,size(val%gas_emission)
          call pmc_mpi_pack_gas_state(buffer, position, val%gas_emission(i))
       end do
    end if
    if (allocated(val%gas_dilution_time)) then
       do i = 1,size(val%gas_background)
          call pmc_mpi_pack_gas_state(buffer, position, val%gas_background(i))
       end do
    end if
    if (allocated(val%aero_emission_time)) then
       do i = 1,size(val%aero_emission)
          call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emission(i))
       end do
    end if
    if (allocated(val%aero_dilution_time)) then
       do i = 1,size(val%aero_background)
          call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background(i))
       end do
    end if
    call assert(639466930, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_pack_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp)
    call pmc_mpi_unpack_real_array(buffer, position, val%pressure_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%pressure)
    call pmc_mpi_unpack_real_array(buffer, position, val%height_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%height)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, &
         val%gas_emission_rate_scale)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, &
         val%aero_emission_rate_scale)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_rate)
    call pmc_mpi_unpack_integer(buffer, position, val%loss_function_type)
    call pmc_mpi_unpack_drydep(buffer, position, val%drydep)
    call pmc_mpi_unpack_chamber(buffer, position, val%chamber)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_deposition_time)
    call pmc_mpi_unpack_real_array(buffer, position, &
         val%gas_deposition_rate_scale)
    call pmc_mpi_unpack_logical(buffer, position, val%do_residual_layer)
    if (allocated(val%gas_deposition)) deallocate(val%gas_deposition)
    if (allocated(val%gas_deposition_time)) then
       allocate(val%gas_deposition(size(val%gas_deposition_time)))
       do i = 1,size(val%gas_deposition)
          call pmc_mpi_unpack_gas_state(buffer, position, &
               val%gas_deposition(i))
       end do
    end if
    if (allocated(val%gas_emission)) deallocate(val%gas_emission)
    if (allocated(val%gas_background)) deallocate(val%gas_background)
    if (allocated(val%aero_emission)) deallocate(val%aero_emission)
    if (allocated(val%aero_background)) deallocate(val%aero_background)
    if (allocated(val%gas_emission_time)) then
       allocate(val%gas_emission(size(val%gas_emission_time)))
       do i = 1,size(val%gas_emission)
          call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emission(i))
       end do
    end if
    if (allocated(val%gas_dilution_time)) then
       allocate(val%gas_background(size(val%gas_dilution_time)))
       do i = 1,size(val%gas_background)
          call pmc_mpi_unpack_gas_state(buffer, position, &
               val%gas_background(i))
       end do
    end if
    if (allocated(val%aero_emission_time)) then
       allocate(val%aero_emission(size(val%aero_emission_time)))
       do i = 1,size(val%aero_emission)
          call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emission(i))
       end do
    end if
    if (allocated(val%aero_dilution_time)) then
       allocate(val%aero_background(size(val%aero_dilution_time)))
       do i = 1,size(val%aero_background)
          call pmc_mpi_unpack_aero_dist(buffer, position, &
               val%aero_background(i))
       end do
    end if
    call assert(611542570, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_unpack_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_drydep(val)

    !> Value to pack.
    type(drydep_params_t), intent(in) :: val

    pmc_mpi_pack_size_drydep = &
         pmc_mpi_pack_size_real(val%z_ref) &
         + pmc_mpi_pack_size_real(val%u_mean) &
         + pmc_mpi_pack_size_real(val%z_rough) &
         + pmc_mpi_pack_size_real(val%A) &
         + pmc_mpi_pack_size_real(val%alpha) &
         + pmc_mpi_pack_size_real(val%eps_0) &
         + pmc_mpi_pack_size_real(val%gamma) &
         + pmc_mpi_pack_size_real(val%C_B) &
         + pmc_mpi_pack_size_real(val%C_IN) &
         + pmc_mpi_pack_size_real(val%C_IM) &
         + pmc_mpi_pack_size_real(val%nu) &
         + pmc_mpi_pack_size_real(val%beta)

   end function pmc_mpi_pack_size_drydep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_drydep(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(drydep_params_t), intent(in) :: val

#ifdef PMC_USE_MPI
      integer :: prev_position

      prev_position = position
      call pmc_mpi_pack_real(buffer, position, val%z_ref)
      call pmc_mpi_pack_real(buffer, position, val%u_mean)
      call pmc_mpi_pack_real(buffer, position, val%z_rough)
      call pmc_mpi_pack_real(buffer, position, val%A)
      call pmc_mpi_pack_real(buffer, position, val%alpha)
      call pmc_mpi_pack_real(buffer, position, val%eps_0)
      call pmc_mpi_pack_real(buffer, position, val%gamma)
      call pmc_mpi_pack_real(buffer, position, val%C_B)
      call pmc_mpi_pack_real(buffer, position, val%C_IN)
      call pmc_mpi_pack_real(buffer, position, val%C_IM)
      call pmc_mpi_pack_real(buffer, position, val%nu)
      call pmc_mpi_pack_real(buffer, position, val%beta)
      call assert(371592468, &
            position - prev_position <= pmc_mpi_pack_size_drydep(val))
#endif

   end subroutine pmc_mpi_pack_drydep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_drydep(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(drydep_params_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%z_ref)
    call pmc_mpi_unpack_real(buffer, position, val%u_mean)
    call pmc_mpi_unpack_real(buffer, position, val%z_rough)
    call pmc_mpi_unpack_real(buffer, position, val%A)
    call pmc_mpi_unpack_real(buffer, position, val%alpha)
    call pmc_mpi_unpack_real(buffer, position, val%eps_0)
    call pmc_mpi_unpack_real(buffer, position, val%gamma)
    call pmc_mpi_unpack_real(buffer, position, val%C_B)
    call pmc_mpi_unpack_real(buffer, position, val%C_IN)
    call pmc_mpi_unpack_real(buffer, position, val%C_IM)
    call pmc_mpi_unpack_real(buffer, position, val%nu)
    call pmc_mpi_unpack_real(buffer, position, val%beta)
    call assert(582741936, &
          position - prev_position <= pmc_mpi_pack_size_drydep(val))
#endif

  end subroutine pmc_mpi_unpack_drydep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_scenario
