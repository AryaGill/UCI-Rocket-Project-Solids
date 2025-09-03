from rocket_setup import *

if __name__ == "__main__":
    env = create_environment()
    rocket = create_rocket(env)

    # rocket.plots.static_margin()
    rocket.draw()


    # Run simulation
    # test_flight = Flight(
    #     rocket=rocket, environment=env, rail_length=5.2, inclination=85, heading=0
    #     )

    # # Print data
    # test_flight.info()
    # # test_flight.prints.initial_conditions()
    # # test_flight.prints.surface_wind_conditions()
    # # test_flight.prints.launch_rail_conditions()
    # # test_flight.prints.out_of_rail_conditions()
    # # test_flight.prints.burn_out_conditions()
    # # test_flight.prints.apogee_conditions()
    # # test_flight.prints.events_registered()
    # # test_flight.prints.impact_conditions()
    # # test_flight.prints.maximum_values()

    # test_flight.all_info()
    # test_flight.plots.trajectory_3d()
    # test_flight.plots.linear_kinematics_data()
    # test_flight.plots.flight_path_angle_data()
    # test_flight.plots.attitude_data()
    # test_flight.plots.angular_kinematics_data()
    # test_flight.plots.aerodynamic_forces()
    # test_flight.plots.rail_buttons_forces()
    # test_flight.plots.energy_data()
    # test_flight.plots.fluid_mechanics_data()
    # test_flight.plots.stability_and_control_data()