#@everywhere const m					= 256+4
#@everywhere const n					= 256+4

@everywhere const m					= 60
@everywhere const n					= 60

@everywhere const n_samples			= 5
@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
#@everywhere const n_samples			= 5
#@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
@everywhere const T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.4

#goldstein sigma <0.5
sig					= 0.0000002

#@everywhere const dt			= 0.2
#@everywhere const dx			= 0.4

# super mit r=dt/dx^2 und 60x60x5_40
@everywhere const dt			= 0.28
@everywhere const dx			= 0.25

# die hier haben mit falscher Gitterkonstante ganz gut funktioniert
#@everywhere const dt			= 0.28 
#@everywhere const dx			= 0.5

@everywhere const alpha	= 0.01
@everywhere const beta	= 0.01

#@everywhere const alpha	= 0.00005
#@everywhere const beta	= 0.00005


maxsteps 			= 100000
save_every			= 0

include("beispiele.jl")
include("verfahren.jl")
include("view.jl")

beta > 0 && @everywhere using IterativeSolvers
# das gmres aus KrylowMethods ist leider schlechter als das von IterativeSolvers
# aber das Gauss-Seidel-Verfahren ist besser.
#beta > 0 && using KrylovMethods  

#H1_norm	= H1_norm_beta
H1_norm		= beta == 0 ? H1_norm_nobeta	: H1_norm_beta
grad_J		= beta == 0 ? grad_J_nobeta		: grad_J_beta_parallel

L2norm		= function(s) return Xnorm(s, B) end
sample_err	= sample_err_L2

s		= inits(quadrat)
#s		= inits(rot_circle)
#s		= readtaxi()[:,:, 1:5:end]

u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

# load old
@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"
run(`mkdir -p $rootdir`)

include("pyamgtest.jl")
ml			= construct_mgsolv(WaveOp)
#grad_J		= grad_J_beta_multig_parallel
grad_J		= grad_J_beta_multig


steps=1
#u, v	= load("$(rootdir)zwischenergebnis_$steps.jld", "u", "v")

#change alpha, beta and run

#@everywhere const alpha= 0.001
#@everywhere const beta	= 0.001
#@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"

@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)
#@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad_goldstein(s, u, v, steps)

_="fertig"
