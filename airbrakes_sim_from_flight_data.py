import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ==============================
# Constants
# ==============================

GAMMA = 1.4
R = 287.05287
g = 9.80665
L = 0.0065

DESIRED_SEARCH_TIME = 20
TIME_PER_SIM_STEP = 0.0132

NUM_RECORDED_DEPLOYMENT_LEVELS = 11
NUM_RECORDED_MACH_NUMS = 14

MASS = 24.2
TARGET_APOGEE_FT = 9000
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
    [0.0089884,0.0085976,0.0084022,0.0084022,0.0084022,0.0084022,0.0084022,0.0084022,0.0084022,0.0085976,0.0085976,0.008793,0.0089884,0.0091838,0.0093792,0.0095746],

    [0.00943484,0.00904404,0.008853104,0.008853104,0.008857569,0.008866498,0.008870962,0.008875426,0.008875426,0.009075291,0.009079755,0.009284084,0.009488413,0.009688277,0.009892606,0.010101399],

    [0.00988128,0.00949048,0.009304009,0.009304009,0.009312938,0.009330795,0.009339724,0.009348653,0.009348653,0.009552982,0.00956191,0.009775168,0.009988426,0.010192754,0.010406012,0.010628198],

    [0.0103277,0.0099369,0.009754893,0.009754893,0.009768286,0.009795072,0.009808465,0.009821858,0.009821858,0.010030651,0.010044044,0.01026623,0.010488416,0.010697209,0.010919395,0.011154974],

    [0.0107741,0.0103833,0.010205757,0.010205757,0.010223614,0.010259328,0.010277185,0.010295042,0.010295042,0.010508299,0.010526156,0.01075727,0.010988384,0.011201641,0.011432755,0.011681726],

    [0.0112206,0.0108298,0.010656722,0.010656722,0.010679044,0.010723688,0.01074601,0.010768332,0.010768332,0.010986054,0.011008376,0.01124842,0.011488464,0.011706186,0.01194623,0.012208596],

    [0.011667,0.0112762,0.011107586,0.011107586,0.011134372,0.011187944,0.01121473,0.011241516,0.011241516,0.011463702,0.011490488,0.01173946,0.011988432,0.012210618,0.01245959,0.012735348],

    [0.0121135,0.0117227,0.011558551,0.011558551,0.011589802,0.011652304,0.011683555,0.011714806,0.011714806,0.011941457,0.011972708,0.01223061,0.012488512,0.012715163,0.012973065,0.013262218],

    [0.0125599,0.0121691,0.012009415,0.012009415,0.01204513,0.01211656,0.012152275,0.01218799,0.01218799,0.012419105,0.01245482,0.01272165,0.01298848,0.013219595,0.013486425,0.01378897],

    [0.0130063,0.0126155,0.012460279,0.012460279,0.012500458,0.012580816,0.012620995,0.012661174,0.012661174,0.012896753,0.012936932,0.01321269,0.013488448,0.013724027,0.013999785,0.014315722],

    [0.0134528,0.013062,0.012911244,0.012911244,0.012955888,0.013045176,0.01308982,0.013134464,0.013134464,0.013374508,0.013419152,0.01370384,0.013988528,0.014228572,0.01451326,0.014842592]
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

def angle_from_vertical(telemetry):
    w = telemetry.q0
    x = telemetry.q1
    y = telemetry.q2
    z = telemetry.q3
    # Compute tilt angle from vertical (radians)
    theta = np.arccos(1 - 2 * (x**2 + y**2))
    # return np.degrees(theta)  # convert to degrees if you want
    return theta

class Telemetry:
    def __init__(self):
        self.time = 0.0
        self.altitude = 0.0
        self.velocity_world_z = 0.0
        self.pressure = 0.0
        self.temperature = 0.0
        self.predicted_apogee = 0.0
        self.airbrake_deployment = 0
        self.q0 = 1.0
        self.q1 = 0.0
        self.q2 = 0.0
        self.q3 = 0.0

def predict_apogee(telemetry, deployment_level):
    deltaT = clampf(
        telemetry.velocity_world_z * deltaT_coefficient,
        0.01,
        0.1
    )

    # Angle from quaternion
    theta = angle_from_vertical(telemetry)
    theta = min(theta, math.radians(80.0))

    # Initial conditions
    alt_sim = telemetry.altitude
    vz_sim = telemetry.velocity_world_z

    # derive horizontal velocity from tilt
    vx_sim = vz_sim * math.tan(theta)

    pressure_Pa = telemetry.pressure

    global ground_temp
    temperature_K = max(ground_temp - (L * telemetry.altitude), 1)

    global sim_angles_from_vert
    sim_angle_from_vert = []
    global sim_times
    sim_time = []
    global sim_pressures
    sim_pressure = []
    global sim_temperatures
    sim_temperature = []
    global sim_velocities
    sim_velocity = []
    global sim_alts
    sim_alt = []
    global sim_velocities_x
    sim_velocity_x = []

    for i in range(100000):
        sim_time.append(telemetry.time + i * deltaT * 1000)

        vz_before = vz_sim

        angle = math.atan2(vx_sim, vz_sim)
        sim_angle_from_vert.append(angle * 180 / 3.14159265358)

        T_local = max(
            temperature_K - (L * (alt_sim - telemetry.altitude)),
            1
        )
        sim_temperature.append(T_local)

        mach = get_mach_number(get_mag2(vz_sim, vx_sim), T_local)
        if angle < 30 * 3.141592653 / 180:
            CdA = get_CdA(deployment_level, mach)
        else:
            CdA = get_CdA(0, mach)

        p_local = pressure_Pa * (T_local / temperature_K) ** (g / (R * L))
        sim_pressure.append(p_local)
        rho = p_local / (R * T_local)

        Fd = 0.5 * CdA * rho * (vx_sim**2 + vz_sim**2)

        Fx = -Fd * math.sin(angle)
        Fz = -Fd * math.cos(angle) - g * MASS

        vx_sim += (Fx / MASS) * deltaT
        vz_sim += (Fz / MASS) * deltaT
        sim_velocity.append(vz_sim)
        sim_velocity_x.append(vx_sim)

        alt_sim += ((vz_sim + vz_before) / 2) * deltaT
        sim_alt.append(alt_sim)

        if vz_sim < 0:
            break

    if telemetry.time - launch_time > motor_burn_time:
        sim_angles_from_vert.append(sim_angle_from_vert)
        sim_times.append(sim_time)
        sim_pressures.append(sim_pressure)
        sim_temperatures.append(sim_temperature)
        sim_velocities.append(sim_velocity)
        sim_alts.append(sim_alt)
        sim_velocities_x.append(sim_velocity_x)
    return alt_sim

def set_optimal_deployment(flight_state, telemetry):
    angle_of_attack = angle_from_vertical(telemetry)

    global ground_temp
    local_temp = max(ground_temp - (L * telemetry.altitude), 1)

    # NOTE: This matches your C code EXACTLY (but is probably wrong physically)
    if (
        flight_state != "GLIDING_ASCENT"
        or get_mach_number(telemetry.velocity_world_z, local_temp) > 0.7
        or angle_of_attack > math.radians(30)
    ):
        telemetry.predicted_apogee = 0
        telemetry.airbrake_deployment = 0
        telemetry.predicted_apogee = predict_apogee(telemetry, 0)
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

    telemetry.predicted_apogee = predict_apogee(telemetry, low)
    telemetry.airbrake_deployment = low

def set_airbrakes_initial_temp(telemetry):
    global ground_temp
    ground_temp = telemetry.temperature

if __name__ == "__main__":
    # read data file
    df = pd.read_excel("night_fury_3-22-26_airbrakes_input.xlsx")
    arr = df.to_numpy()

    time = df["time"].to_numpy()
    baro_alt = df["altitude"].to_numpy()
    startAlt = df["startAlt"].to_numpy()
    alt_fused = df["alt_fused"].to_numpy()
    initial_temp = [x + 273.71 for x in df["initial_temp"].to_numpy()]
    velocity_world_z = df["velocity_world_z"].to_numpy()
    # velocity_world_z = []
    # for i in range(len(baro_alt)):
    #     if i < len(baro_alt) - 1 and i != 0:
    #         velocity_world_z.append((baro_alt[i+1] - baro_alt[i-1]) / 2)
    #     else:
    #         velocity_world_z.append(0)
    q0 = df["q0"].to_numpy()
    q1 = df["q1"].to_numpy()
    q2 = df["q2"].to_numpy()
    q3 = df["q3"].to_numpy()
    pressure = [x*100 for x in df["pressure"].to_numpy()]
    temperature = [x + 273.71 for x in df["temperature"].to_numpy()]

    altitude = []
    for i in range(len(baro_alt)):
        altitude.append(baro_alt[i] - startAlt[i])

    launch_time = 652889
    motor_burn_time = 4700
    
    # idea: try using measured temp and see if changes much

    # create arrays to plot
    predicted_apogee = []
    deployment_level = []
    angle_from_vert = []
    global sim_angles_from_vert
    sim_angles_from_vert = []
    global sim_times
    sim_times = []
    global sim_pressures
    sim_pressures = []
    global sim_temperatures
    sim_temperatures = []
    global sim_velocities
    sim_velocities = []
    global sim_alts
    sim_alts = []
    global sim_velocities_x
    sim_velocities_x = []

    tau_velocity = 0.3
    tau_altitude = 0.5
    alpha_alt = 0.3
    # alpha_vel = 1 - tau_altitude
    alpha_vel = 1
    # alpha_quat = 1 - tau_altitude - tau_velocity
    alpha_quat = 1
    # alpha_pres = 1 - tau_altitude - tau_velocity
    alpha_pres = 1
    altitude_old = []
    for i in range(len(altitude)):
        if i != 0:
            altitude_old.append(altitude[i] * alpha_alt + altitude[i-1] * (1-alpha_alt))
        else:
            altitude_old.append(0)
    alt_fused_old = []
    for i in range(len(altitude)):
        # alt_fused_old.append(alt_fused[i] - 0.7 * velocity_world_z[i])
        if i != 0:
            alt_fused_old.append(alt_fused[i] * alpha_alt + alt_fused[i-1] * (1-alpha_alt))
        else:
            alt_fused_old.append(0)
    velocity_world_z_old = []
    for i in range(len(altitude)):
        if i != 0:
            velocity_world_z_old.append(velocity_world_z[i] * alpha_vel + velocity_world_z[i-1] * (1-alpha_vel))
        else:
            velocity_world_z_old.append(0)
    q0_old = []
    for i in range(len(q0)):
        if i != 0:
            q0_old.append(q0[i] * alpha_quat + q0[i-1] * (1-alpha_quat))
        else:
            q0_old.append(1)
    q1_old = []
    for i in range(len(q0)):
        if i != 0:
            q1_old.append(q1[i] * alpha_quat + q1[i-1] * (1-alpha_quat))
        else:
            q1_old.append(0)
    q2_old = []
    for i in range(len(q2)):
        if i != 0:
            q2_old.append(q2[i] * alpha_quat + q2[i-1] * (1-alpha_quat))
        else:
            q2_old.append(0)
    q3_old = []
    for i in range(len(q3)):
        if i != 0:
            q3_old.append(q3[i] * alpha_quat + q3[i-1] * (1-alpha_quat))
        else:
            q3_old.append(0)
    pressure_old = []
    for i in range(len(altitude)):
        if i != 0:
            pressure_old.append(pressure[i] * alpha_pres + pressure[i-1] * (1-alpha_pres))
        else:
            pressure_old.append(0)
    

    for i in range(len(time)):
        telemetry = Telemetry()
        telemetry.time = time[i]
        # telemetry.altitude = altitude[i]
        # telemetry.altitude = altitude_old[i]
        # telemetry.altitude = alt_fused[i]
        telemetry.altitude = alt_fused_old[i]
        # telemetry.velocity_world_z = velocity_world_z[i]
        telemetry.velocity_world_z = velocity_world_z_old[i]
        telemetry.temperature = initial_temp[i]
        # telemetry.pressure = pressure[i]
        telemetry.pressure = pressure_old[i]
        # telemetry.q0 = q0[i]
        # telemetry.q1 = q1[i]
        # telemetry.q2 = q2[i]
        # telemetry.q3 = q3[i]
        telemetry.q0 = q0_old[i]
        telemetry.q1 = q1_old[i]
        telemetry.q2 = q2_old[i]
        telemetry.q3 = q3_old[i]

        set_airbrakes_initial_temp(telemetry)
        if time[i] - launch_time > motor_burn_time:
            set_optimal_deployment("GLIDING_ASCENT", telemetry)
        else:
            set_optimal_deployment("MOTOR_BURN", telemetry)

        if time[i] - launch_time > motor_burn_time:
            deployment_level.append(min(telemetry.airbrake_deployment / (NUM_DEPLOYMENT_LEVELS - 1), 0.851))
            # air_brakes.deployment_level = 0 /  (NUM_DEPLOYMENT_LEVELS - 1)
            # telemetry.predicted_apogee = predict_apogee(telemetry, 0)
        else:
            deployment_level.append(0)
        
        telemetry.predicted_apogee = predict_apogee(telemetry, 0)
        predicted_apogee.append(telemetry.predicted_apogee)
        angle_from_vert.append(angle_from_vertical(telemetry) * 180 / 3.14159265358)

    idx_airbrakes_off_bc_angle = next((i for i, x in enumerate(angle_from_vert) if x > 30), None)
    idx_motor_burn_end = next((i for i, x in enumerate(time) if x >= launch_time + motor_burn_time), None)
    idx_apogee = next((i for i, x in enumerate(altitude) if x == max(altitude)), None)
        
    # Plot predicted apogee
    plt.plot(time, [x*3.2808399 for x in predicted_apogee], label="Predicted Apogee")
    plt.plot(time, [x*3.2808399 for x in altitude], label="Altitude")
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.axvline(x=time[idx_airbrakes_off_bc_angle], color='g', linestyle='--', linewidth=2, label="Approx. Airbrakes off bc angle > 30")
    plt.xlabel("Time (ms)")
    plt.ylabel("Altitude (ft)")
    plt.title("Predicted Apogee by Time")
    plt.legend()
    plt.grid()
    plt.show()

    # # Plot deployment level
    # plt.plot(time, deployment_level, label="Deployment Level")
    # plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    # plt.xlabel("Time (ms)")
    # plt.ylabel("Deployment Level")
    # plt.title("Deployment Level by Time")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # Plot overestimate
    overestimate = [(predicted_apogee[i] - max(altitude)) * 3.2808399 for i in range(len(time))]
    print("Max Overestimate: " + str(max(overestimate)))
    plt.plot(time[idx_motor_burn_end:], overestimate[idx_motor_burn_end:], label="Predicted Apogee - Real Apogee")
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.axvline(x=time[idx_apogee], color='g', linestyle='--', linewidth=2, label="Apogee")
    plt.xlabel("Time (ms)")
    plt.ylabel("Altitude (ft)")
    plt.title("Predicted Apogee Overestimate by Time")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot angle from vert
    plt.plot(time, angle_from_vert, label="Measured Angle from Vertical")
    for i in range(len(sim_times)):
        plt.plot(sim_times[i], sim_angles_from_vert[i])
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.axhline(y=30, color='g', linestyle='--', linewidth=2, label="30 degrees (airbrakes off when above)")
    plt.xlabel("Time (ms)")
    plt.ylabel("Angle (deg)")
    plt.title("Angle from vertical by Time")
    plt.legend()
    plt.grid()
    plt.show()

    # # Plot pressure
    # plt.plot(time, pressure, label="Measured Pressure")
    # for i in range(len(sim_times)):
    #     plt.plot(sim_times[i], sim_pressures[i])
    # plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    # plt.xlabel("Time (ms)")
    # plt.ylabel("Pressure (hPa)")
    # plt.title("Pressure by Time")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # # Plot temperature
    # plt.plot(time, temperature, label="Measured Temperature")
    # for i in range(len(sim_times)):
    #     plt.plot(sim_times[i], sim_temperatures[i])
    # plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    # plt.xlabel("Time (ms)")
    # plt.ylabel("Temperature (K)")
    # plt.title("Temperature by Time")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # Plot velocity world z
    for i in range(len(sim_times)):
        plt.plot([x + 700 for x in sim_times[i]], sim_velocities[i])
    plt.plot(time, velocity_world_z, label="Measured Vecicty World Z")
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.xlabel("Time (ms)")
    plt.ylabel("Veclocity (m/s)")
    plt.title("Velocity by Time")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot velocity world x
    for i in range(len(sim_times)):
        plt.plot(sim_times[i], sim_velocities_x[i])
    plt.plot(time, angle_from_vert, label="Measured Angle from Vertical")
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.xlabel("Time (ms)")
    plt.ylabel("Veclocity (m/s)")
    plt.title("Velocity by Time")
    plt.legend()
    plt.grid()
    plt.show()

    # Plot altitude
    for i in range(len(sim_times)):
        plt.plot(sim_times[i], [x*3.2808399 for x in sim_alts[i]])
    plt.plot(time[idx_motor_burn_end:], [x*3.2808399 for x in altitude[idx_motor_burn_end:]], label="Measured Altitude")
    plt.axvline(x=launch_time + motor_burn_time, color='r', linestyle='--', linewidth=2, label="Approx. Motor Burn End")
    plt.xlabel("Time (ms)")
    plt.ylabel("Altitude (ft)")
    plt.title("Altitude by Time")
    plt.legend()
    plt.grid()
    plt.show()