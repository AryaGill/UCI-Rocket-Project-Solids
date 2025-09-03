import datetime

class config:
    def __init__(self):
        # Environment
        self.latitude = 32.990254
        self.longitude = -106.974998
        self.elevation = 1400

        tomorrow = datetime.date.today() + datetime.timedelta(days=1)
        self.year = tomorrow.year
        self.month = tomorrow.month
        self.day = tomorrow.day

        self.atmosphere_model_type = "Forecast"
        self.atmosphere_model_file = "GFS"

        # Motor
        # self.thrust_source = "../data/motors/projeto-jupiter/keron_thrust_curve.csv"
        # self.dry_mass = 1.815
        # self.dry_inertia = (0.125, 0.125, 0.002)
        # self.center_of_dry_mass_position = 0.317
        # self.grains_center_of_mass_position = 0.397
        self.burn_time = 4.625
        self.grain_number = 7
        self.grain_separation = 0.0127000254
        self.grain_density = 1750
        self.grain_outer_radius = 0.04113538227
        # self.grain_initial_inner_radius = 15 / 1000
        self.grain_initial_height = 0.1397002794
        self.nozzle_radius = 0.01231902464
        self.throat_radius = 0.01231902464
        self.interpolation_method = "linear"
        self.nozzle_position = 0
        self.coordinate_system_orientation = "nozzle_to_combustion_chamber"

        # Rocket
        self.radius = 0.07886715773
        self.mass = 21.25985704
        # self.inertia = (6.321, 6.321, 0.034)
        self.power_off_drag = "azula_power_on_drag.csv"
        self.power_on_drag = "azula_power_on_drag.csv"
        self.center_of_mass_without_motor = 1.833197866
        self.coordinate_system_orientation = "tail_to_nose"

        self.motor_position = -1.255

        # Rail Buttons
        self.upper_button_position = 2.387604775
        self.lower_button_position = 2.84480569
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
        self.fin_position = 3.149606299
        self.cant_angle = 0
        # self.airfoil = ("../data/airfoils/NACA0012-radians.txt","radians")

        # Tail
        self.tail_top_radius = 0.07886715773
        self.tail_bottom_radius = 0.05270510541
        self.tail_length = 0.1651003302
        self.tail_position = 3.606807214

        # Main Parachute
        self.main_name = "main"
        self.main_cd_s = 6.6902
        self.main_trigger = 304.7851265      # ejection altitude in meters
        self.main_sampling_rate = 105
        # self.main_lag = 1.5
        # self.main_noise = (0, 8.3, 0.5)

        # Drogue Parachute
        self.drogue_name = "drogue"
        self.drogue_cd_s = 1.79955
        self.drogue_trigger = "apogee"  # ejection at apogee
        self.drogue_sampling_rate = 105
        # self.drogue_lag = 1.5
        # self.drogue_noise = (0, 8.3, 0.5)

        # Airbrakes
        self.drag_coefficient_curve = "../data/rockets/calisto/air_brakes_cd.csv"
        self.sampling_rate = 10
        self.reference_area = None
        self.clamp = True
        self.initial_observed_variables = [0, 0, 0]
        self.override_rocket_drag = False
        self.name = "Air Brakes"