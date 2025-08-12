from __future__ import annotations

from typing import Dict, List

from PyQt5 import QtCore


# Placeholder values for solar geometry
sun_center_x_list_placeholder: List[float] = [2056.060059]
sun_center_y_list_placeholder: List[float] = [2043.719971]
pixel_arcsec_x_list_placeholder: List[float] = [0.599489]
pixel_arcsec_y_list_placeholder: List[float] = [0.599489]
GOES_graph_placeholder_list: List[str] = ["logarithmic"]
wavelength_list: List[int] = [1600, 94, 131, 171, 193, 211, 304, 335]

# Frame selection lists
not_main_finalized_frame_list: List[int] = []
finalized_frame_list: List[int] = []
auto_verification: List[bool] = []
AEC_list: List[float] = []
key: List[str] = []
manual_selected_time_range: List[str] = []
changing_solar_movies_dict: Dict[str, object] = {}

# UI state flags
flare_selected: bool = False
one_frame_state: bool = False
graph_type_regular_bool: bool = True
graph_type_changed_bool: bool = False
slider_changed: bool = False
animation_interval_changed: bool = False
user_clicked_map_bool: bool = False
coord_type_changed_bool: bool = True
playing: bool = True

# Numeric configuration
min_value: int = 0
max_value_vmax: int = 10000
max_threshold: int = 4000
current_frame_number: int = 0
label_distance: int = 20
initial_index: int = 0
GOES_scale_value: float = 0.6982
cropped_dimension: int = 15
frame_timer: int = 20

translate = QtCore.QCoreApplication.translate
