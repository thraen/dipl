#@everywhere const m					= 256+4
#@everywhere const n					= 256+4
@everywhere const m					= 140
@everywhere const n					= 140
#@everywhere const m					= 60
#@everywhere const n					= 60

#@everywhere const m					= 30
#@everywhere const n					= 30

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples			= 5
@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
@everywhere const T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

#multigrid solver tolerance
@everywhere const mg_tol = 1e-1

#@everywhere const dt			= 1/T
#@everywhere const dx			= 1/m

@everywhere const dt			= 0.01
@everywhere const dx			= 0.01

@everywhere const alpha	= 0.001
@everywhere const beta	= 0.001

maxsteps 			= 5
save_every			= 0


# vielleicht lieber so
#method	= "time_regularization__center_velocies"
#method	= "no_time_regularization__center_velocies"
#method	= "no_time_regularization__interf_velocies"
#method	= "no_time_regularization__interf_velocies_project"

velocities_at		= "centers"
#velocities_at		= "interfaces"

transport_paralell	= false
grad_parallel		= true

time_regularization	= false

project_divfree		= false

#poisson_solver	= "multig" 
#poisson_solver	= "gmres"
poisson_solver	= "lufact" #fur gegebene Probleme am besten

#stokes_solver	= "multig"
#stokes_solver	= "gmres"
stokes_solver	= "lufact" #fur gegebene Probleme am besten

timereg_solver	= "multig" #fur gegebene Probleme am besten
#timereg_solver	= "gmres"
#timereg_solver	= "lufact"


#include("view.jl")

include("beispiele.jl")

include("transport.jl")

include("verfahren.jl")

#s		= inits(quadrat)
#s		= inits(rot_circle)
s		= inits(rot_circle_ex)[:,:,1:5]
#s		= load_taxi(m,n,T)[:,:, 1:5:end]

u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

# bei Zellzwischenwerten
#u		= 0* ones( m, n-1, T-1 )
#v		= 0* ones( m-1, n, T-1 )
#u       = SharedArray(Float64, (m, n-1, T-1))#, init= S -> S[localindexes(S)] = 0.0)
#v       = SharedArray(Float64, (m-1, n, T-1))

# load old
@everywhere rootdir = "../out/new/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)_mgtol$(mg_tol)/"
run(`mkdir -p $rootdir`)
run(`sh -c "cp *jl $rootdir"`)

steps=1
#u, v	= load("$(rootdir)zwischenergebnis_$steps.jld", "u", "v")
#change alpha, beta and run
#@everywhere const alpha= 0.001
#@everywhere const beta	= 0.001
#@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"

@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)

_="fertig"
