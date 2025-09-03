from rocketpy import Environment, SolidMotor, Rocket, Flight
from RocketPy.azula.azula_config import config
from simple_pid import PID

def create_environment():
    env = Environment(latitude=config.latitude, longitude=config.longitude, elevation=config.elevation)

    env.set_date(
        (config.year, config.month, config.day, 12)
    )  # Hour given in UTC time

    env.set_atmospheric_model(type=config.atmosphere_model_type, file=config.atmosphere_model_file)

    return env

def create_rocket(env):
    # Setup Motor
    motor = SolidMotor(
        thrust_source=config.thrust_source,
        dry_mass=config.dry_mass,
        dry_inertia=config.dry_inertia,
        center_of_dry_mass_position=config.center_of_dry_mass_position,
        grains_center_of_mass_position=config.grains_center_of_mass_position,
        burn_time=config.burn_time,
        grain_number=config.grain_number,
        grain_separation=config.grain_separation,
        grain_density=config.grain_density,
        grain_outer_radius=config.grain_outer_radius,
        grain_initial_inner_radius=config.grain_initial_inner_radius,
        grain_initial_height=config.grain_initial_height,
        nozzle_radius=config.nozzle_radius,
        throat_radius=config.throat_radius,
        interpolation_method=config.interpolation_method,
        nozzle_position=config.nozzle_position,
        coordinate_system_orientation=config.coordinate_system_orientation,
    )

    # motor.info()


    # Setup Rocket
    rocket = Rocket(
        radius=config.radius,
        mass=config.mass,
        inertia=config.inertia,
        power_off_drag=config.power_off_drag,
        power_on_drag=config.power_on_drag,
        center_of_mass_without_motor=config.center_of_mass_without_motor,
        coordinate_system_orientation=config.coordinate_system_orientation,
    )

    rocket.add_motor(motor, position=config.motor_position)

    rail_buttons = rocket.set_rail_buttons(
        upper_button_position=config.upper_button_position,
        lower_button_position=config.lower_button_position,
        angular_position=config.angular_position,
    )

    nose_cone = rocket.add_nose(
        length=config.nose_cone_length, kind=config.nose_cone_kind, position=config.nose_cone_position
    )

    fin_set = rocket.add_trapezoidal_fins(
        n=config.num_fins,
        root_chord=config.root_chord,
        tip_chord=config.tip_chord,
        span=config.fin_span,
        position=config.fin_position,
        cant_angle=config.cant_ange,
        airfoil=config.airfoil,
    )

    tail = rocket.add_tail(
        top_radius=config.tail_top_radius,
        bottom_radius=config.tail_bottom_radius,
        length=config.tail_length,
        position=config.tail_position
    )

    main = rocket.add_parachute(
        name=config.main_name,
        cd_s=config.main_cd_s,
        trigger=config.main_trigger,      # ejection altitude in meters
        sampling_rate=config.main_sampling_rate,
        lag=config.main_lag,
        noise=config.main_noise,
    )

    drogue = rocket.add_parachute(
        name=config.drogue_name,
        cd_s=config.drogue_cd_s,
        trigger=config.drogue_trigger,  # ejection at apogee
        sampling_rate=config.drogue_sampling_rate,
        lag=config.drogue_lag,
        noise=config.drogue_noise,
    )

    # air_brakes = rocket.add_air_brakes(
    #     drag_coefficient_curve=config.drag_coefficient_curve,
    #     controller_function=controller_function,
    #     sampling_rate=config.air_brakes_sampling_rate,
    #     reference_area=config.air_brakes_reference_area,
    #     clamp=config.air_brakes_clamp,
    #     initial_observed_variables=config.air_brakes_initial_observed_variables,
    #     override_rocket_drag=config.air_brakes_override_rocket_drag,
    #     name=config.airbrakes_name,
    # )

    TARGET_APOGEE_FT = 10000000
    TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048
    pid = PID(Pk=0.01, Ki=0.0, Kd=0.001, setpoint=TARGET_APOGEE_M)
    pid.output_limits = (0.0, 1.0)

    def controller_function(
        time, sampling_rate, state, state_history, observed_variables, air_brakes
    ):
        # state = [x, y, z, vx, vy, vz, e0, e1, e2, e3, wx, wy, wz]
        altitude_ASL = state[2]
        altitude_AGL = altitude_ASL - env.elevation
        vx, vy, vz = state[3], state[4], state[5]

        # Get winds in x and y directions
        wind_x, wind_y = env.wind_velocity_x(altitude_ASL), env.wind_velocity_y(altitude_ASL)

        # Calculate Mach number
        free_stream_speed = (
            (wind_x - vx) ** 2 + (wind_y - vy) ** 2 + (vz) ** 2
        ) ** 0.5
        mach_number = free_stream_speed / env.speed_of_sound(altitude_ASL)

        # Get previous state from state_history
        previous_state = state_history[-1]
        previous_vz = previous_state[5]

        # If we wanted to we could get the returned values from observed_variables:
        # returned_time, deployment_level, drag_coefficient = observed_variables[-1]

        # Check if the rocket has reached burnout
        if time < rocket.burn_out_time:
            return None

        # If below 1500 meters above ground level, air_brakes are not deployed
        if altitude_AGL < 1500:
            air_brakes.deployment_level = 0

        # Else calculate the deployment level
        else:
            # Controller logic
            new_deployment_level = (
                air_brakes.deployment_level + 0.1 * vz + 0.01 * previous_vz**2
            )

            # Limiting the speed of the air_brakes to 0.2 per second
            # Since this function is called every 1/sampling_rate seconds
            # the max change in deployment level per call is 0.2/sampling_rate
            max_change = 0.2 / sampling_rate
            lower_bound = air_brakes.deployment_level - max_change
            upper_bound = air_brakes.deployment_level + max_change
            new_deployment_level = min(max(new_deployment_level, lower_bound), upper_bound)

            air_brakes.deployment_level = new_deployment_level

        # Return variables of interest to be saved in the observed_variables list
        return (
            time,
            air_brakes.deployment_level,
            air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number),
        )
    
    return rocket