#rmprocs( procs()[2:end] )
#addprocs(1)
#addprocs(5)

# @everywhere const m					= 256+4
# @everywhere const n					= 256+4
# @everywhere const m					= 140
# @everywhere const n					= 140
@everywhere const m					= 60
@everywhere const n					= 60

#@everywhere const m					= 30
#@everywhere const n					= 30

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples			= 5
#@everywhere const n_zwischensamples	= 9    # duerfen nicht zu wenige sein? abhaengig von dt?
@everywhere const n_zwischensamples	= 39    # duerfen nicht zu wenige sein? abhaengig von dt?

# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
@everywhere const T					= (n_samples-1)*(n_zwischensamples+1) +1
@show T
# Zuordnung Samplenummer zu Zeitpunkt 
@show sample_times		= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples-1 ]

# what, das ist doch falsch!
# @everywhere const T					= (n_samples-1)*n_zwischensamples+1
#sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

#@everywhere const dt			= 1/T
#@everywhere const dx			= 1/m

@everywhere const dt			= 0.01
# @everywhere const dt			= 0.04
@everywhere const dx			= 0.01

@everywhere const alpha	= 0.001  #*0.001 #warum ist das nicht dasselbe, wie die norm noch mal mit alpha zu multiplizieren? siehe CostNormOp --> thr ahh, wegen der ruecksubstitution nach ellipt gleichung. aber warum funktioniert es so gut, wenn man die norm noch mal mit alpha multipliziert. teste das auch mal ohne zeitreg
@everywhere const beta	= 0.0001 #*0.001

# maxsteps 			= 2
maxsteps 			= 100000

save_every			= 0

time_regularization	= true  # geht nicht mit velocities_at interfaces

# velocities_at		= "interfaces"
velocities_at		= "centers"

transport_parallel	= false # geht nicht gut, verbesserungswuerdig

# das Verfahren mit Zeitregularisierung parallelisiert 
# automatisch die Dimensionen, wenn mehr als ein Worker existiert

grad_parallel		= false # betrifft nur die Verfahren ohne Zeitregularisierung

project_divfree		= false

#thr diese Optionen funktionieren nicht alle und die meisten sind sowieso unsinnvoll.
#poisson_solver	= "multig"  #geht nicht parallel
#poisson_solver	= "gmres"	#ungenau
poisson_solver	= "lufact" #fur gegebene Probleme am besten. Eigentlich Cholesky-Faktorisierung fuer die interfaces und LU-Faktorisierung fuer center

#stokes_solver	= "multig"	#schlecht geeignet, langsam
#stokes_solver	= "gmres"	#ungenau
stokes_solver	= "lufact"#fur gegebene Probleme am besten

timereg_solver	= "multig"#fur gegebene Probleme am besten
#timereg_solver	= "gmres"	#ungenau
#timereg_solver	= "lufact"	#nur fuer sehr kleine Probleme benutzbar

#multigrid solver tolerance
@everywhere const mg_tol = 1e-1 

@everywhere with_cfl_check = false

# Zeitregularisierung funktioniert nur mit Flussdiskretisierung an Zellmittelpunkten
# diese Zeile ist zu Sicherheit, damit man nichts falsch einstellt
velocities_at		= ~time_regularization ? velocities_at : "centers"

#include("view.jl")
include("beispiele.jl")

s		= inits(quadrat)
# s		= inits(rot_circle)
# s		= inits(rot_circle_ex)[:,:,1:n_samples]
# s		= load_taxi(m,n,41)[:,:, 1:5:end]
# s		= readtaxi_alt()[:,:, 1:5:5*n_samples+1]

velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= 0* ones( m-1, n, T-1 )
end

include("verfahren.jl") 

# bei Zellzwischenwerten
#u       = SharedArray(Float64, (m, n-1, T-1))#, init= S -> S[localindexes(S)] = 0.0)
#v       = SharedArray(Float64, (m-1, n, T-1))

@everywhere rootdir = "../out/new/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)_mgtol$(mg_tol)/"
run(`mkdir -p $rootdir/src`)
run(`sh -c "cp *jl $rootdir/src"`)
run(`sh -c "git log -1 > $rootdir/this_git_commit"`) #thr

steps=1
#u, v	= load("$(rootdir)zwischenergebnis_$steps.jld", "u", "v")
#change alpha, beta and run
#@everywhere const alpha= 0.001
#@everywhere const beta	= 0.001
#@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"

@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)

_="fertig"
