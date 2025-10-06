from rocket_setup import *
from azula.azula_config import config
import matplotlib.pyplot as plt

if __name__ == "__main__":
    conf = config()
    env = create_environment(conf)
    rocket = create_rocket(conf, env)

    # rocket.plots.static_margin()
    # rocket.draw()

    # Run simulation
    test_flight = Flight(
        rocket=rocket, environment=env, rail_length=5.2, inclination=90, heading=0, terminate_on_apogee=True
        )
    
    test_flight.aerodynamic_drag()

    # # Print data
    # test_flight.info()
    # test_flight.prints.initial_conditions()
    # test_flight.prints.surface_wind_conditions()
    # test_flight.prints.launch_rail_conditions()
    # test_flight.prints.out_of_rail_conditions()
    # test_flight.prints.burn_out_conditions()
    test_flight.prints.apogee_conditions()
    # test_flight.prints.events_registered()
    # test_flight.prints.impact_conditions()
    # test_flight.prints.maximum_values()

    # test_flight.all_info()
    test_flight.plots.trajectory_3d()
    # test_flight.plots.linear_kinematics_data()
    # test_flight.plots.flight_path_angle_data()
    # test_flight.plots.attitude_data()
    # test_flight.plots.angular_kinematics_data()
    # test_flight.plots.aerodynamic_forces()
    # test_flight.plots.rail_buttons_forces()
    # test_flight.plots.energy_data()
    # test_flight.plots.fluid_mechanics_data()
    # test_flight.plots.stability_and_control_data()

    # Airbrakes data
    time_list, deployment_level_list, drag_coefficient_list, predicted_apogee_list = [], [], [], []
    obs_vars = test_flight.get_controller_observed_variables()
    for time, deployment_level, drag_coefficient, predicted_apogee in obs_vars:
        time_list.append(time)
        deployment_level_list.append(deployment_level)
        drag_coefficient_list.append(drag_coefficient)
        predicted_apogee_list.append(predicted_apogee)

    # Plot deployment level by time
    plt.plot(time_list, deployment_level_list)
    plt.xlabel("Time (s)")
    plt.ylabel("Deployment Level")
    plt.title("Deployment Level by Time")
    plt.grid()
    plt.show()

    # Plot drag coefficient by time
    plt.plot(time_list, drag_coefficient_list)
    plt.xlabel("Time (s)")
    plt.ylabel("Drag Coefficient")
    plt.title("Drag Coefficient by Time")
    plt.grid()
    plt.show()

    # Plot predicted apogee by time
    plt.plot(time_list, [x / 0.3048 for x in predicted_apogee_list])
    plt.xlabel("Time (s)")
    plt.ylabel("Predicted Apogee (ft)")
    plt.title("Predicted Apogee by Time")
    plt.grid()
    plt.show()