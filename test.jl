@everywhere const m					= 60
@everywhere const n					= 60

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples			= 5
@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
@everywhere const T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_bas			= 0.5
armijo_sig			= 0.0

@everywhere const dx			= 0.25
@everywhere const dt			= 0.28

# die hier haben mit falscher Gitterkonstante ganz gut funktioniert
#@everywhere const dx			= 0.5
#@everywhere const dt			= 0.28 

@everywhere const alpha	= 0.01
@everywhere const beta	= 0.01

maxsteps 			= 100000

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

sample_err	= sample_err_L2
L2norm		= function(s) return Xnorm(s, B) end

#s		= inits(quadrat)
s		= inits(rot_circle)

u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

@everywhere rootdir = "../out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"
run(`mkdir -p $rootdir`)

@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v)

#=
testabs = [ (0.01,  0.01),
			(0.01,  0.011),
		   	(0.011, 0.01),
		   	(0.011, 0.011) ]

for testab in testabs
	alpha	= testab[1]
	beta	= testab[2]

	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )

	@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v)
	#@time I, u, v, p, L2_err, H1_err, J, steps = verfahren_direkt(s, u, v)

	save_all()
	#save_images_(s, "s")
end
=#

_="fertig"
