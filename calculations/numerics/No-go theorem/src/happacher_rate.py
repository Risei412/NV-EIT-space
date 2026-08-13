"""Alias shim: this project's phonon-rate module was renamed to phonon_rates.py;
bperp_kernel_map_v2.py still imports it as happacher_rate (its original name)."""
from phonon_rates import *
from phonon_rates import k_orb, ALPHA_PH, ALPHA_PH_ERR
