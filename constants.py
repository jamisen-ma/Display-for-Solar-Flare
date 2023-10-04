from PyQt5 import QtCore


sun_center_x_list_placeholder = [2056.060059]
sun_center_y_list_placeholder = [2043.719971]
pixel_arcsec_x_list_placeholder = [0.599489]
pixel_arcsec_y_list_placeholder = [0.599489]
GOES_graph_placeholder_list = ["logarithmic"]
wavelength_list = [1600, 94, 131, 171, 193, 211, 304, 335]

not_main_finalized_frame_list = []
finalized_frame_list = []
auto_verification = []
AEC_list = []
key = []
manual_selected_time_range = []
changing_solar_movies_dict = {}

flare_selected = False
one_frame_state = False
graph_type_regular_bool = True
graph_type_changed_bool = False
slider_changed = False
animation_interval_changed = False
user_clicked_map_bool = False
coord_type_changed_bool = True
playing = True

min_value = 0
max_value_vmax = 10000
max_threshold = 4000
current_frame_number = 0
label_distance = 20
initial_index = 0
GOES_scale_value = 0.6982
cropped_dimension = 15
frame_timer = 20

translate = QtCore.QCoreApplication.translate
