"""
..module:: QTAG_config
  :platform: Unix, Windows
  :synopsis: This module contains various dictionaries to be used as input
             selections for QTAG calculations.
..moduleauthors :: Matthew Dutra
"""

univ = {
"ndof" : 1,
"ntraj" : 15,
"dt" : 1,
"niter" : 5000,
"mass" : [1743.0],
"n_snapshots" : 4,
"n_data_out" : 10
}

wf0 = {
"q" : [1.76],
"p" : [0.0],
"a" : [16.16],
"s" : [0.0]
}

traj0 = {
"placement" : "grid",
"grid_dims" : [15],
"rho" : 1e-12,
"a0" : [12.0]
}

basis = {
"qtype" : "adpt",
"ptype" : "adpt",
"atype" : "adpt",
"stype" : "frzn"
}

mss = {
"prop_method" : "sync",
"decpl" : 0.9
}

mom_params = {
"adjust" : "lin_fit",
"beta" : 1e-2
}

model = {
"pot_type" : "nafh",
"rep" : "adiab",
"calc_type" : "BAT",
"coupling" : "BAT"
}

model_params = {
"rNaF" : 3.779,
"a" : 0.01,
"b" : 1.147,
"d1" : [0.005],
"d2" : [1.0],
"d3" : [0.0]
}

