@everywhere m					=  60
@everywhere n					=  60

@everywhere n_samples			=   5
@everywhere n_zwischensamples	=  40   # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

const dx			= 0.5
const dt			= 0.27 #0.99833

@everywhere alpha	= 0.01

maxsteps 			= 1000

#@everywhere rootdir = "$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)/"
@everywhere rootdir = "/tmp/out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)/"

include("beispiele.jl")

include("verfahren.jl")
