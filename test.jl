@everywhere const m					= 256+4
@everywhere const n					= 256+4
#@everywhere const m					= 60
#@everywhere const n					= 60

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples			= 5
@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
@everywhere const T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

#goldstein sigma <0.5
sig					= 0.4

#@everywhere const dt			= 1/T
#@everywhere const dx			= 1/m

#@everywhere const dt			= 1/m # thr. warum ist das besser als richtig?
#@everywhere const dx			= 1/T

#@everywhere const dt			= 0.005
#@everywhere const dx			= 0.0001

@everywhere const dt			= 0.01
@everywhere const dx			= 0.01

# super mit r=dt/dx^2 und 60x60x5_40
#@everywhere const dt			= 0.0028
#@everywhere const dx			= 0.0025

@everywhere const alpha	= 0.001
@everywhere const beta	= 0.001

maxsteps 			= 10000
save_every			= 250

include("beispiele.jl")
include("verfahren.jl")
include("view.jl")

#include("transport.jl")
#include("transport_neu.jl")
include("transport_interfaces.jl")

# ohne beta
grad_J		= grad_J_nobeta_interf
H1_norm_grd	= H1_norm_nobeta_interf
H1_norm_w	= H1_norm_nobeta_interf

#grad_J		= grad_J_beta
#H1_norm_grd	= H1_norm_beta_grd
#H1_norm_w	= H1_norm_beta_w

L2norm		= function(s) return Xnorm(s, B) end
sample_err	= sample_err_L2

#s		= inits(quadrat)
#s		= inits(rot_circle)
s		= readtaxi()[:,:, 1:5:end]

#u		= 0* ones( m, n, T-1 )
#v		= 0* ones( m, n, T-1 )

# bei Zellzwischenwerten
u		= 0* ones( m, n-1, T-1 )
v		= 0* ones( m-1, n, T-1 )

# load old
@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"
run(`mkdir -p $rootdir`)

steps=1
#u, v	= load("$(rootdir)zwischenergebnis_$steps.jld", "u", "v")
#change alpha, beta and run
#@everywhere const alpha= 0.001
#@everywhere const beta	= 0.001
#@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"

@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)
#@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad_goldstein(s, u, v, steps)

_="fertig"
