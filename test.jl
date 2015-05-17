@everywhere m					=  60
@everywhere n					=  60

@everywhere n_samples			=   5
@everywhere n_zwischensamples	=  40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

const dx			= 0.5
const dt			= 0.95 #0.99833

@everywhere alpha	= 0.01
@everywhere beta	= 0.01

maxsteps 			= 1000

#@everywhere rootdir = "$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)/"
@everywhere rootdir = "/tmp/out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)/"

include("beispiele.jl")

include("verfahren.jl")
include("view.jl")

#include("marcel_matrizen.jl")
#save_value(full(L2), "L2")


#s		= inits(quadrat)
s		= inits(rot_circle)


u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(maxsteps, alpha, s, u, v, L2norm, H1_norm, sample_err_L2)

#I, u, v, p, L2_err, H1_err, J, steps = verfahren_direkt(maxsteps, alpha, s, u, v, L2norm, H1_norm, sample_err_L2)
#save_all()
#save_images_(s, "s")

_="fertig"
