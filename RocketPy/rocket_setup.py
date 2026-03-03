from rocketpy import Environment, SolidMotor, Rocket, Flight
import math
import numpy as np
import matplotlib.pyplot as plt

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

    
    def angle_from_vertical(e0, e1, e2, e3):
        w, x, y, z = e0, e1, e2, e3
        # Compute tilt angle from vertical (radians)
        theta = np.arccos(1 - 2 * (x**2 + y**2))
        # return np.degrees(theta)  # convert to degrees if you want
        return theta

    TARGET_APOGEE_FT = 10000
    TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048

    def binary_search_deployment(alt, vz, air_brakes, T0, pressure0, angle_of_attack, speed0):
        num_sims = 10
        
        low = 0
        high = 1

        for i in range(num_sims):
            mid = (high + low) / 2
            if predict_apogee(alt, vz, air_brakes, T0, pressure0, mid, angle_of_attack, speed0) > TARGET_APOGEE_M:
                low = mid
            else:
                high = mid
        return (high + low) / 2

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
    # deltaT = 0.1
    # gamma = 1.4
    # R = 287.05287
    # g = 9.80665
    # L = 0.0065   # K/m lapse rate
    # A = (config.radius ** 2) * math.pi
    # wanted_time = 10 # ms
    # time_per_call = 0.0125 # ms
    # deltaT_coefficient = (time_per_call / wanted_time) / g
    # def predict_apogee(alt, vz, air_brakes, T0, pressure0, deployment_level, angle_of_attack, speed0):

    #     deltaT = max(0.01, min(speed0 * deltaT_coefficient * math.cos(angle_of_attack), 0.1))
    #     print(deltaT)
        
    #     # v_sim = vz
    #     alt_sim = alt

    #     vz_sim = speed0 * math.cos(angle_of_attack)
    #     vx_sim = speed0 * math.sin(angle_of_attack)

    #     for i in range(0, 3000):

    #         # Calculate Mach Number
    #         # T_local = max(T0 - (L * (alt_sim - alt)), 1.0)
    #         # speed_of_sound = math.sqrt(gamma * R * T_local)
    #         # mach_number = v_sim / speed_of_sound

    #         # # airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number) + rocket.power_off_drag(mach_number)
    #         # # airbrake_Cd = rocket.power_on_drag(mach_number)
    #         # airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number)

    #         # p_local = pressure0 * (T_local / T0) ** (g / (R * L))
    #         # rho_sim = p_local / (R * T_local)

    #         # F = -0.5*airbrake_Cd*rho_sim*A*(v_sim**2) - g*(config.mass + config.dry_mass)
    #         # a_sim = F / (config.mass + config.dry_mass)
    #         # v_sim += a_sim * deltaT
    #         # alt_sim += v_sim * deltaT



    #         # Calculate Mach Number
    #         vz_sim_before = vz_sim
    #         T_local = max(T0 - (L * (alt_sim - alt)), 1.0)
    #         speed_of_sound = math.sqrt(gamma * R * T_local)
    #         mach_number = math.sqrt(vz_sim ** 2 + vx_sim ** 2) / speed_of_sound

    #         # airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number) + rocket.power_off_drag(mach_number)
    #         # airbrake_Cd = rocket.power_off_drag(mach_number)
    #         airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number)

    #         p_local = pressure0 * (T_local / T0) ** (g / (R * L))
    #         rho_sim = p_local / (R * T_local)

    #         Fd = 0.5*airbrake_Cd*rho_sim*A*(vx_sim**2 + vz_sim ** 2)
    #         angle_sim = math.atan2(vx_sim, vz_sim)
    #         Fx = -Fd * math.sin(angle_sim)
    #         Fz = -Fd * math.cos(angle_sim) - g*(config.mass + config.dry_mass)
    #         vx_sim += (Fx / (config.mass + config.dry_mass)) * deltaT
    #         vz_sim += (Fz / (config.mass + config.dry_mass)) * deltaT
    #         alt_sim += ((vz_sim + vz_sim_before) / 2) * deltaT



    #         # Angle, no change
    #         # vz_sim_before = vz_sim
    #         # T_local = max(T0 - (L * (alt_sim - alt)), 1.0)
    #         # speed_of_sound = math.sqrt(gamma * R * T_local)
    #         # mach_number = (vz_sim / math.cos(angle_of_attack)) / speed_of_sound

    #         # # airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number) + rocket.power_off_drag(mach_number)
    #         # airbrake_Cd = air_brakes.drag_coefficient(deployment_level, mach_number)

    #         # p_local = pressure0 * (T_local / T0) ** (g / (R * L))
    #         # rho_sim = p_local / (R * T_local)

    #         # Fd = 0.5*airbrake_Cd*rho_sim*A*((vz_sim / math.cos(angle_of_attack)) ** 2)
    #         # Fz = -Fd * math.cos(angle_of_attack) - g*(config.mass + config.dry_mass)
    #         # vz_sim += (Fz / (config.mass + config.dry_mass)) * deltaT
    #         # alt_sim += ((vz_sim + vz_sim_before) / 2) * deltaT

    #         if vz_sim < 0:
    #             break

    #     return alt_sim



    # TARGET_APOGEE_FT = 10000
    # TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048
    # GAMMA = 1.4
    # R = 287.05287
    # g = 9.80665          # Gravity
    # L = 0.0065           # Temperature Lapse Rate
    # MASS = 20
    # WANTED_AIRBRAKE_ALG_TIME = 30      # ms
    # TIME_PER_AIRBRAKE_CALL = 0.0125    # ms

    # deltaT = 0.01
    # A = (0.07886715773 ** 2) * math.pi
    # deltaT_coefficient = (TIME_PER_AIRBRAKE_CALL / WANTED_AIRBRAKE_ALG_TIME) / g
    # deployment = 0

    # NUM_RECORDED_DEPLOYMENT_LEVELS = 11
    # NUM_RECORDED_MACH_NUMS = 8
    # NUM_DEPLOYMENT_LEVELS = 1024

    # air_brakes_drag_coefficient = [
    #     [0.400, 0.400, 0.400, 0.400, 0.400, 0.400, 0.400, 0.400],  # 0.0
    #     [0.430, 0.460, 0.490, 0.520, 0.550, 0.580, 0.610, 0.640],  # 0.1
    #     [0.460, 0.490, 0.520, 0.550, 0.580, 0.610, 0.640, 0.670],  # 0.2
    #     [0.490, 0.520, 0.550, 0.580, 0.610, 0.640, 0.670, 0.700],  # 0.3
    #     [0.520, 0.550, 0.580, 0.610, 0.640, 0.670, 0.700, 0.730],  # 0.4
    #     [0.550, 0.580, 0.610, 0.640, 0.670, 0.700, 0.730, 0.760],  # 0.5
    #     [0.580, 0.610, 0.640, 0.670, 0.700, 0.730, 0.760, 0.790],  # 0.6
    #     [0.610, 0.640, 0.670, 0.700, 0.730, 0.760, 0.790, 0.820],  # 0.7
    #     [0.640, 0.670, 0.700, 0.730, 0.760, 0.790, 0.820, 0.850],  # 0.8
    #     [0.670, 0.700, 0.730, 0.760, 0.790, 0.820, 0.850, 0.880],  # 0.9
    #     [0.700, 0.730, 0.760, 0.790, 0.820, 0.850, 0.880, 0.910],  # 1.0
    # ]

    # def get_drag_coefficient(deployment_level: int, mach_number: float) -> float:
    #     if (mach_number >= 0.7):
    #         return air_brakes_drag_coefficient[NUM_RECORDED_DEPLOYMENT_LEVELS - 1][NUM_RECORDED_MACH_NUMS - 1]

    #     mach_idx = mach_number * (NUM_RECORDED_MACH_NUMS - 1) / 0.7
    #     deployment_idx = deployment_level * (NUM_RECORDED_DEPLOYMENT_LEVELS - 1) / (NUM_DEPLOYMENT_LEVELS - 1)

    #     mach_i = int(mach_idx)
    #     mach_frac = mach_idx - mach_i

    #     deployment_i = int(deployment_idx)
    #     deployment_frac = deployment_idx - deployment_i

    #     low_low = air_brakes_drag_coefficient[deployment_i][mach_i]
    #     low_high = air_brakes_drag_coefficient[deployment_i][mach_i + 1]

    #     if deployment_idx == math.floor(deployment_idx):
    #         return low_low + (low_high - low_low) * mach_frac

    #     high_low = air_brakes_drag_coefficient[deployment_i + 1][mach_i]
    #     high_high = air_brakes_drag_coefficient[deployment_i + 1][mach_i + 1]

    #     x = mach_frac
    #     y = deployment_frac

    #     # Bilinear interpolation
    #     result = (
    #         (1 - x) * (1 - y) * low_low +
    #         x * (1 - y) * low_high +
    #         (1 - x) * y * high_low +
    #         x * y * high_high
    #     )

    #     return result


    # # ---------------------------------------------------------
    # #   MACH NUMBER
    # # ---------------------------------------------------------
    # def get_mach_number(velocity: float, temp: float) -> float:
    #     speed_of_sound = math.sqrt(R * GAMMA * temp)
    #     return velocity / speed_of_sound


    # # ---------------------------------------------------------
    # #   PREDICT APOGEE
    # # ---------------------------------------------------------
    # def predict_apogee(alt: float, temp0: float, pressure0: float,
    #                     angle_of_attack: float, speed0: float,
    #                     deployment_level: int) -> float:
    #     global deltaT

    #     deltaT = max(0.01, min(speed0 * deltaT_coefficient * math.cos(angle_of_attack), 0.1))

    #     alt_sim = alt
    #     vz_sim = speed0 * math.cos(angle_of_attack)
    #     vx_sim = speed0 * math.sin(angle_of_attack)

    #     for _ in range(100000):
    #         vz_before = vz_sim

    #         T_local = max(temp0 - (L * (alt_sim - alt)), 1)

    #         mach_number = get_mach_number(math.sqrt(vz_sim**2 + vx_sim**2), T_local)
    #         airbrake_Cd = get_drag_coefficient(deployment_level, mach_number)

    #         p_local = pressure0 * (T_local / temp0) ** (g / (R * L))
    #         rho_sim = p_local / (R * T_local)

    #         Fd = 0.5 * airbrake_Cd * rho_sim * A * (vx_sim**2 + vz_sim**2)

    #         angle_sim = math.atan2(vx_sim, vz_sim)
    #         Fx = -Fd * math.sin(angle_sim)
    #         Fz = -Fd * math.cos(angle_sim) - g * MASS

    #         vx_sim += (Fx / MASS) * deltaT
    #         vz_sim += (Fz / MASS) * deltaT

    #         alt_sim += ((vz_sim + vz_before) / 2) * deltaT

    #         if vz_sim < 0:
    #             break

    #     return alt_sim


    # # ---------------------------------------------------------
    # #   BINARY SEARCH FOR OPTIMAL DEPLOYMENT
    # # ---------------------------------------------------------
    # def optimal_deployment(alt: float, temp0: float, pressure0: float,
    #                     angle_of_attack: float, speed0: float) -> int:

    #     if angle_of_attack > math.radians(30):
    #         return 0

    #     low = 0
    #     high = NUM_DEPLOYMENT_LEVELS - 1

    #     num_sims = int(math.log2(NUM_DEPLOYMENT_LEVELS))

    #     for _ in range(num_sims):
    #         mid = (high + low) // 2

    #         if predict_apogee(alt, temp0, pressure0, angle_of_attack, speed0, mid) > TARGET_APOGEE_M:
    #             low = mid
    #         else:
    #             high = mid

    #     return low
    




    # def controller_function(
    #     time, sampling_rate, state, state_history, observed_variables, air_brakes
    # ):
    #     # state = [x, y, z, vx, vy, vz, e0, e1, e2, e3, wx, wy, wz]
    #     altitude_ASL = state[2]
    #     altitude_AGL = altitude_ASL - env.elevation
    #     vx, vy, vz = state[3], state[4], state[5]
    #     speed = math.sqrt(vx ** 2 + vz ** 2 + vy ** 2)

    #     angle_of_attack = angle_from_vertical(state[6], state[7], state[8], state[9])

    #     # Get winds in x and y directions
    #     wind_x, wind_y = env.wind_velocity_x(altitude_ASL), env.wind_velocity_y(altitude_ASL)

    #     # Calculate Mach number
    #     free_stream_speed = (
    #         (wind_x - vx) ** 2 + (wind_y - vy) ** 2 + (vz) ** 2
    #     ) ** 0.5
    #     mach_number = free_stream_speed / env.speed_of_sound(altitude_ASL)

    #     speed_of_sound = (GAMMA * R * env.temperature(altitude_ASL)) ** 0.5
    #     my_mach_number = speed / speed_of_sound

    #     # Check if the rocket has reached burnout
    #     if time < motor.burn_time[1] or vz <= 0 or angle_of_attack > 30 / 180 * math.pi:
    #         air_brakes.deployment_level = 0
    #         return None
        
        
    #     # Binary Search Method
    #     # new_deployment_level = binary_search_deployment(altitude_AGL, vz, air_brakes, env.temperature(altitude_ASL), env.pressure(altitude_ASL), angle_of_attack, speed)
    #     new_deployment_level = optimal_deployment(altitude_AGL, env.temperature(altitude_ASL), env.pressure(altitude_ASL), angle_of_attack, speed) / NUM_DEPLOYMENT_LEVELS

    #     # PID Method
    #     predicted_apogee = predict_apogee(altitude_AGL, env.temperature(altitude_ASL), env.pressure(altitude_ASL), angle_of_attack, speed, air_brakes.deployment_level * NUM_DEPLOYMENT_LEVELS)
    #     # new_deployment_level = min(1, max(0, air_brakes.deployment_level + pid.update(predicted_apogee, time)))
    #     # new_deployment_level = 0.3

    #     # Limiting the speed of the air_brakes to 0.2 per second
    #     # Since this function is called every 1/sampling_rate seconds
    #     # the max change in deployment level per call is 0.2/sampling_rate
    #     max_change = 0.2 / sampling_rate
    #     lower_bound = air_brakes.deployment_level - max_change
    #     upper_bound = air_brakes.deployment_level + max_change
    #     new_deployment_level = min(max(new_deployment_level, lower_bound), upper_bound)

    #     air_brakes.deployment_level = new_deployment_level

    #     # Return variables of interest to be saved in the observed_variables list
    #     return (
    #         time,
    #         air_brakes.deployment_level,
    #         air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number),
    #         predicted_apogee,
    #     )
    
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












    # Night Fury test algorithm

    # ==============================
    # Constants
    # ==============================

    GAMMA = 1.4
    R = 287.05287
    g = 9.80665
    L = 0.0065

    DESIRED_SEARCH_TIME = 200
    TIME_PER_SIM_STEP = 0.015

    NUM_RECORDED_DEPLOYMENT_LEVELS = 11
    NUM_RECORDED_MACH_NUMS = 14

    MASS = 23.28
    TARGET_APOGEE_FT = 8000
    TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048

    NUM_DEPLOYMENT_LEVELS = 64
    NUM_RECORDED_DEPLOYMENT_LEVELS = 11
    NUM_RECORDED_MACH_NUMS = 14

    deltaT_coefficient = (
        TIME_PER_SIM_STEP
        * math.log2(NUM_DEPLOYMENT_LEVELS)
        / DESIRED_SEARCH_TIME
        / g
    )

    ground_temp = 300.0

    # ==============================
    # Airbrake CdA Table
    # ==============================

    air_brakes_CdA = [
        [0.0070731,0.00721944,0.007317,0.00741456,0.00743895,0.00743895,0.00746334,0.00748773,0.00748773,0.00748773,0.00751212,0.00753651,0.0075609,0.0075609],
        [0.007467,0.00761634,0.0077159,0.00776568,0.00786524,0.00789013,0.00793991,0.0079648,0.00798969,0.00801458,0.00803947,0.00806436,0.00806436,0.00808925],
        [0.00765355,0.00786328,0.00792168,0.00794961,0.00797754,0.00800547,0.0080334,0.00806133,0.00808925,0.00811972,0.00814765,0.00817558,0.00820351,0.00823144],
        [0.00800113,0.00823566,0.0083261,0.00839045,0.00841109,0.00843173,0.00845238,0.00847304,0.00849369,0.00851433,0.00853497,0.00855561,0.00857626,0.0085969],
        [0.0082368,0.008316,0.0083952,0.0084744,0.0085536,0.0086328,0.008712,0.0087912,0.0088704,0.0089496,0.0090288,0.009108,0.0091872,0.0092664],
        [0.008608,0.0087963,0.0089039,0.0089577,0.0090115,0.0090384,0.0090922,0.0090922,0.0090922,0.0090922,0.0090922,0.0090922,0.0090922,0.0091191],
        [0.0088228,0.00887486,0.00892618,0.0089775,0.00902882,0.00908014,0.00913146,0.00918278,0.0092341,0.00928542,0.00933674,0.00938806,0.00943938,0.0094907],
        [0.0092907,0.0095418,0.0096534,0.0097092,0.009765,0.0098208,0.0098766,0.0099324,0.0099882,0.010044,0.0100719,0.0100998,0.0101277,0.0101556],
        [0.00940371,0.00968781,0.00980145,0.00985827,0.00991509,0.00997191,0.01002873,0.01008555,0.01014237,0.01019919,0.0102276,0.01025601,0.01028442,0.01031283],
        [0.0098294,0.01003177,0.01023314,0.01034878,0.0104066,0.01046442,0.01052224,0.01055115,0.01058006,0.01060897,0.0106957,0.01078243,0.01078243,0.01078243],
        [0.0102935,0.01055719,0.01067583,0.01076406,0.01082288,0.01070524,0.01091011,0.01093952,0.01093952,0.01096893,0.01105716,0.01108657,0.01111598,0.01114539],
    ]

    # ==============================
    # Utility Functions
    # ==============================

    def clampf(x, minimum, maximum):
        return max(minimum, min(x, maximum))


    def get_CdA(deployment_level, mach):
        deployment = deployment_level / (NUM_DEPLOYMENT_LEVELS - 1)
        deployment = clampf(deployment, 0.0, 1.0)
        mach = clampf(mach, 0.05, 0.7)

        dep_index_f = deployment * 10.0
        mach_index_f = (mach - 0.05) * 20.0

        dep_index = int(dep_index_f)
        mach_index = int(mach_index_f)

        if dep_index >= NUM_RECORDED_DEPLOYMENT_LEVELS - 1:
            dep_index = NUM_RECORDED_DEPLOYMENT_LEVELS - 2
        if mach_index >= NUM_RECORDED_MACH_NUMS - 1:
            mach_index = NUM_RECORDED_MACH_NUMS - 2

        dep_frac = dep_index_f - dep_index
        mach_frac = mach_index_f - mach_index

        CdA00 = air_brakes_CdA[dep_index][mach_index]
        CdA10 = air_brakes_CdA[dep_index+1][mach_index]
        CdA01 = air_brakes_CdA[dep_index][mach_index+1]
        CdA11 = air_brakes_CdA[dep_index+1][mach_index+1]

        CdA0 = CdA00 + dep_frac * (CdA10 - CdA00)
        CdA1 = CdA01 + dep_frac * (CdA11 - CdA01)

        return CdA0 + mach_frac * (CdA1 - CdA0)


    def get_mach_number(velocity, temp):
        speed_of_sound = math.sqrt(R * GAMMA * temp)
        return velocity / speed_of_sound


    def get_mag2(x, y):
        return math.sqrt(x*x + y*y)


    def get_mag3(x, y, z):
        return math.sqrt(x*x + y*y + z*z)


    # ==============================
    # Telemetry Class
    # ==============================

    class Telemetry:
        def __init__(self):
            self.altitude = 0.0
            self.velocity_world_x = 0.0
            self.velocity_world_y = 0.0
            self.velocity_world_z = 0.0
            self.pressure = 0.0
            self.temperature = 0.0
            self.predicted_apogee = 0.0
            self.airbrake_deployment = 0


    # ==============================
    # Apogee Prediction
    # ==============================
    def predict_apogee(telemetry, deployment_level):
        global has_done
        global ground_temp

        deltaT = clampf(
            telemetry.velocity_world_z * deltaT_coefficient,
            0.01,
            0.1
        )

        alt_sim = telemetry.altitude
        vz_sim = telemetry.velocity_world_z
        vx_sim = get_mag2(
            telemetry.velocity_world_x,
            telemetry.velocity_world_y
        )

        pressure_Pa = telemetry.pressure
        temperature_K = max(ground_temp - (L * telemetry.altitude), 1)

        for _ in range(100000):

            vz_before = vz_sim
            T_local = max(temperature_K - (L * (alt_sim - telemetry.altitude)), 1)

            mach = get_mach_number(get_mag2(vz_sim, vx_sim), T_local)
            CdA = get_CdA(deployment_level, mach)

            p_local = pressure_Pa * (T_local / temperature_K) ** (g / (R * L))
            rho = p_local / (R * T_local)

            Fd = 0.5 * CdA * rho * (vx_sim**2 + vz_sim**2)

            angle = math.atan2(vx_sim, vz_sim)

            Fx = -Fd * math.sin(angle)
            Fz = -Fd * math.cos(angle) - g * MASS

            vx_sim += (Fx / MASS) * deltaT
            vz_sim += (Fz / MASS) * deltaT

            alt_sim += ((vz_sim + vz_before) / 2) * deltaT

            if vz_sim < 0:
                break

        return alt_sim


    # ==============================
    # Optimal Deployment (Binary Search)
    # ==============================

    def set_optimal_deployment(flight_state, telemetry):

        horizontal_speed = get_mag2(
            telemetry.velocity_world_x,
            telemetry.velocity_world_y
        )

        angle_of_attack = math.atan2(
            horizontal_speed,
            telemetry.velocity_world_z
        )

        local_temp = max(ground_temp - (L * telemetry.altitude), 1)

        if (flight_state != "GLIDING_ASCENT"
            or get_mach_number(
                get_mag3(
                    telemetry.velocity_world_x,
                    telemetry.velocity_world_y,
                    telemetry.velocity_world_z
                ),
                local_temp
            ) > 0.7
            or angle_of_attack > math.radians(30)
        ):
            telemetry.predicted_apogee = 0
            telemetry.airbrake_deployment = 0
            return

        low = 0
        high = NUM_DEPLOYMENT_LEVELS - 1

        num_sims = int(math.log2(NUM_DEPLOYMENT_LEVELS))

        for _ in range(num_sims):
            mid = ((high + low) // 2) + 1
            pred_apogee = predict_apogee(telemetry, mid)

            if pred_apogee >= TARGET_APOGEE_M:
                low = mid
            else:
                high = mid - 1

        telemetry.predicted_apogee = pred_apogee
        telemetry.airbrake_deployment = low

    def set_airbrakes_initial_temp(telemetry):
        global ground_temp
        ground_temp = telemetry.temperature




    def controller_function(
        time, sampling_rate, state, state_history, observed_variables, air_brakes
    ):
        global ground_temp

        # state = [x, y, z, vx, vy, vz, e0, e1, e2, e3, wx, wy, wz]
        telemetry = Telemetry()
        telemetry.altitude = state[2] - env.elevation
        telemetry.velocity_world_x = state[3]
        telemetry.velocity_world_y = state[4]
        telemetry.velocity_world_z = state[5]
        telemetry.temperature = env.temperature(env.elevation)
        telemetry.pressure = env.pressure(state[2])

        set_airbrakes_initial_temp(telemetry)
        set_optimal_deployment("GLIDING_ASCENT", telemetry)

        air_brakes.deployment_level = telemetry.airbrake_deployment / (NUM_DEPLOYMENT_LEVELS - 1)
        # air_brakes.deployment_level = 63 /  (NUM_DEPLOYMENT_LEVELS - 1)
        # telemetry.predicted_apogee = predict_apogee(telemetry, 63)

        # Return variables of interest to be saved in the observed_variables list
        local_temp = max(ground_temp - (L * telemetry.altitude), 1)
        mach_number = get_mach_number(
                get_mag3(
                    telemetry.velocity_world_x,
                    telemetry.velocity_world_y,
                    telemetry.velocity_world_z
                ),
                local_temp
            )

        return (
            time,
            air_brakes.deployment_level,
            air_brakes.drag_coefficient(air_brakes.deployment_level, mach_number),
            telemetry.predicted_apogee,
            mach_number,
        )
    
    air_brakes = rocket.add_air_brakes(
        drag_coefficient_curve=config.drag_coefficient_curve,
        controller_function=controller_function,
        sampling_rate=config.air_brakes_sampling_rate,
        reference_area=config.air_brakes_reference_area,
        clamp=config.air_brakes_clamp,
        initial_observed_variables=config.air_brakes_initial_observed_variables,
        override_rocket_drag=config.air_brakes_override_rocket_drag,
        name=config.air_brakes_name,
    )












    # air_brakes.all_info()
    
    return rocket