import datetime

class config:
    def __init__(self):
        # Environment
        self.latitude = 31.0498
        self.longitude = -103.5473
        self.elevation = 0

        tomorrow = datetime.date.today() + datetime.timedelta(days=1)
        self.year = tomorrow.year
        self.month = tomorrow.month
        self.day = tomorrow.day

        # self.year = 2025
        # self.month = 6
        # self.day = 13

        self.atmosphere_model_type = "Forecast"
        # self.atmosphere_model_type = "standard_atmosphere"
        self.atmosphere_model_file = "GFS"

        # Motor
        self.thrust_source = "RocketPy/azula/azula_thrust_curve.csv"
        self.dry_mass = 7.07377301
        self.dry_inertia = (6.300, 6.300, 0.042) # Guess
        self.center_of_dry_mass_position = 1.833197866 # Guess
        self.grains_center_of_mass_position = 0.527051054 # Guess
        self.burn_time = 4.625
        self.grain_number = 7
        self.grain_separation = 0.0127000254
        self.grain_density = 1750
        self.grain_outer_radius = 0.04113538227
        self.grain_initial_inner_radius = 0.01817 # Guess
        self.grain_initial_height = 0.1397002794
        self.nozzle_radius = 0.01231902464
        self.throat_radius = 0.01231902464
        self.interpolation_method = "linear"
        self.nozzle_position = 0
        self.motor_coordinate_system_orientation = "nozzle_to_combustion_chamber"

        # Rocket
        self.radius = 0.07886715773
        self.mass = 15.591
        self.inertia = (10.020, 10.020, 0.067) # Guess
        self.power_off_drag = "RocketPy/azula/azula_power_on_drag.csv"
        self.power_on_drag = "RocketPy/azula/azula_power_on_drag.csv"
        self.center_of_mass_without_motor = -1.833197866
        self.rocket_coordinate_system_orientation = "tail_to_nose"

        self.motor_position = -3.606807214

        # Rail Buttons
        self.upper_button_position = -2.387604775
        self.lower_button_position = -2.84480569
        self.angular_position = 45

        # Nose Cone
        self.nose_cone_length = 0.9144018288
        self.nose_cone_kind = "von karman"
        self.nose_cone_position = 0

        # Fins
        self.num_fins = 4
        self.root_chord = 0.2794005588
        self.tip_chord = 0.1524003048
        self.fin_span = 0.1778003556
        self.fin_position = -3.149606299
        self.cant_angle = 0
        # self.airfoil = ("../data/airfoils/NACA0012-radians.txt","radians") # Guess

        # Tail
        self.tail_top_radius = 0.07886715773
        self.tail_bottom_radius = 0.05270510541
        self.tail_length = 0.1651003302
        self.tail_position = -3.606807214

        # Main Parachute
        self.main_name = "main"
        self.main_cd_s = 6.6902
        self.main_trigger = 304.7851265      # ejection altitude in meters
        self.main_sampling_rate = 105
        self.main_lag = 1.5 # Guess
        self.main_noise = (0, 8.3, 0.5) # Guess

        # Drogue Parachute
        self.drogue_name = "drogue"
        self.drogue_cd_s = 1.79955
        self.drogue_trigger = "apogee"  # ejection at apogee
        self.drogue_sampling_rate = 105
        self.drogue_lag = 1.5 # Guess
        self.drogue_noise = (0, 8.3, 0.5) # Guess

        # Airbrakes
        self.drag_coefficient_curve = "RocketPy/azula/azula_air_brakes_drag_coefficient_curve.csv" # Guess
        self.air_brakes_sampling_rate = 10
        self.air_brakes_reference_area = None
        self.air_brakes_clamp = True
        self.air_brakes_initial_observed_variables = [0, 0, 0, 0]
        self.air_brakes_override_rocket_drag = False
        self.air_brakes_name = "Air Brakes"