"""
..module:: QTAG_config
  :platform: Unix, Windows
  :synopsis: This module contains various dictionaries to be used as input
             selections for QTAG calculations.
..moduleauthors :: Matthew Dutra
"""

univ = {
"ntraj" : 25,
"dt" : 1,
"niter" : 2000,
"mass" : 2000.0,
"nout" : 1
}

wf0 = {
"q" : -5.0,
"p" : 10.0,
"a" : 2.0,
"s" : 0.0
}

traj0 = {
"type" : "grid",
"rho" : 1e-8,
"a0" : 6.0
}

basis = {
"qtype" : "adpt",
"ptype" : "adpt",
"atype" : "adpt",
"stype" : "frzn"
}

mss = {
"type" : "sync",
"decpl" : 0.9
}

mom_params = {
"type" : "lin_fit",
"beta" : 1e-1
}

model = {
"pottype" : "T1",
"rep" : "diab",
"calctype" : "LHA",
"cplg" : "LHA"
}

model_params = {
"a" : 0.01,
"b" : 1.147,
"d1" : 0.005,
"d2" : 1.0,
"d3" : 0.0
}
