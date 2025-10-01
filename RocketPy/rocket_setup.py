from rocketpy import Environment, SolidMotor, Rocket, Flight
import math

class PID:
    def __init__(self, Kp, Ki, Kd, setpoint=0.0):
        self.Kp = Kp  # proportional gain
        self.Ki = Ki  # integral gain
        self.Kd = Kd  # derivative gain
        self.setpoint = setpoint

        self.integral = 0.0
        self.prev_error = 0.0
        self.prev_time = None

    def update(self, measurement, current_time):
        """Compute PID output given a new measurement."""
        
        # Calculate delta time
        if self.prev_time is None:
            dt = 0.0
        else:
            dt = current_time - self.prev_time

        # Error term
        error = measurement - self.setpoint

        # Integral term
        self.integral += error * dt

        # Derivative term
        derivative = (error - self.prev_error) / dt if dt > 0 else 0.0

        # PID output
        output = (self.Kp * error) + (self.Ki * self.integral) + (self.Kd * derivative)

        # Save for next step
        self.prev_error = error
        self.prev_time = current_time

        return output

def create_environment(config):
    env = Environment(latitude=config.latitude, longitude=config.longitude, elevation=config.elevation)

    env.set_date(
        (config.year, config.month, config.day, 12)
    )  # Hour given in UTC time

    env.set_atmospheric_model(type=config.atmosphere_model_type, file=config.atmosphere_model_file)

    return env

def create_rocket(config, env):
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
        coordinate_system_orientation=config.motor_coordinate_system_orientation,
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
        coordinate_system_orientation=config.rocket_coordinate_system_orientation,
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
        cant_angle=config.cant_angle,
        # airfoil=config.airfoil,
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

    TARGET_APOGEE_FT = 10000
    TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048

    pid = PID(Kp=0.1, Ki=0.01, Kd=0.1, setpoint=TARGET_APOGEE_M)

    # def predict_apogee(cur_alt, vz):
        # Kinematic Equation:
        # return cur_alt + (vz ** 2) / (2 * 9.8)
    
        # Ballistic Coasting Method
        # m   = config.mass     # current mass [kg]

        # if vz <= 0:
        #     return alt  # already past apogee

        # g_eff = 9.80665 + D / m
        # delta_h = vz**2 / (2 * g_eff)

    # Using Drag
    # def predict_apogee(flight, alt, vz, airbrake_Cd):
        # rho = 1.2
        # A = (config.radius ** 2) * math.pi
        # Cd = airbrake_Cd # + rocket_Cd
        # k = 0.5*rho*Cd*A
        # m = config.mass
        # g = 9.8
        # delta_h = (m / (2*k)) * math.log((m* g + k * (vz ** 2)) / (m * g))

        # return alt + delta_h
    
    
    # Simulation
    deltaT = 0.01
    gamma = 1.4
    R = 287
    g = 9.81
    L = 0.0065   # K/m lapse rate
    def predict_apogee(alt, vz, air_brakes, T0, pressure0):
        v_sim = vz
        alt_sim = alt
        A = (config.radius ** 2) * math.pi
        for i in range(0, 3000):

            # Calculate Mach Number
            # T_local = T0
            T_local = max(T0 - L * alt_sim, 1.0)
            speed_of_sound = (gamma * R * T_local) ** 0.5
            mach_number = v_sim / speed_of_sound

            airbrake_Cd = air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number)

            p_local = pressure0 * (T_local / T0) ** (g / (R * L))
            rho_sim = p_local / (R * T_local)

            F = -0.5*airbrake_Cd*rho_sim*A*(v_sim**2) - g*config.mass
            a_sim = F / config.mass
            v_sim = v_sim + a_sim * deltaT
            alt_sim += v_sim * deltaT

            if v_sim < 0:
                break

        return alt_sim

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
        if time < motor.burn_time[1] or vz <= 0:
            air_brakes.deployment_level = 0
            return None
        
        # predicted_apogee = predict_apogee(rocket, altitude_AGL, vz, air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number))
        predicted_apogee = predict_apogee(altitude_AGL, vz, air_brakes, env.temperature(altitude_ASL), env.pressure(altitude_ASL))
        # predicted_apogee = predict_apogee(altitude_AGL, vz)
        new_deployment_level = min(1, max(0, air_brakes.deployment_level + pid.update(predicted_apogee, time)))
        # new_deployment_level = 0

        # Controller logic
        # new_deployment_level = (
        #     air_brakes.deployment_level + 0.1 * vz + 0.01 * previous_vz**2
        # )

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
            predicted_apogee,
        )
    
    # air_brakes = rocket.add_air_brakes(
    #     drag_coefficient_curve=config.drag_coefficient_curve,
    #     controller_function=controller_function,
    #     sampling_rate=config.air_brakes_sampling_rate,
    #     reference_area=config.air_brakes_reference_area,
    #     clamp=config.air_brakes_clamp,
    #     initial_observed_variables=config.air_brakes_initial_observed_variables,
    #     override_rocket_drag=config.air_brakes_override_rocket_drag,
    #     name=config.air_brakes_name,
    # )

    # air_brakes.all_info()
    
    return rocket