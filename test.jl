armijo_bas			= 0.5
armijo_sig			= 0.0

@everywhere const alpha	= 0.0009 
@everywhere const beta	= 0.0009

# maxsteps 			= 2
maxsteps 			= 100000

save_every			= 0

time_regularization	= false  # geht nicht mit velocities_at interfaces

# velocities_at		= "interfaces"
velocities_at		= "centers"

transport_parallel	= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll

# das Verfahren mit Zeitregularisierung parallelisiert 
# automatisch die Dimensionen, wenn mehr als ein Worker existiert

grad_parallel		= false # betrifft nur die Verfahren ohne Zeitregularisierung

project_divfree		= false # betrifft nur velocities_at = "interfaces"

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

@everywhere with_cfl_check = true

# Zeitregularisierung funktioniert nur mit Flussdiskretisierung an Zellmittelpunkten
# diese Zeile ist zu Sicherheit, damit man nichts falsch einstellt
# velocities_at		= ~time_regularization ? velocities_at : "centers"

include("view.jl")

@everywhere const m					= 100
@everywhere const n					= 100

include("beispiele.jl")

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples				= 2

@everywhere const auslassen				= 7 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 5 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

# die Anzahl zwischen den Referenzframes zu generierenden Frames. 
@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

# ...................... T, alle Frames/Zeitpunkte, also T-1 Zeitschritte von einem Frame auf den naechsten
@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1

# Zuordnung Samplenummer zu Zeitpunkt 
@everywhere		  samples_to_frames		= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples-1 ]

@everywhere const vorgabe_used_indices	= (1:(auslassen+1):(auslassen+1)*n_samples) 
@everywhere const T_vorgabe				= vorgabe_used_indices[end]
@everywhere		  samples_to_vorgabe	= [(k, vorgabe_used_indices[k]) for k in 1:n_samples]

@everywhere const vorgabe_frames		= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
# @everywhere const vorgabe_to_frames		= [(k,vorgabe_frames[k]) for k in 1:T_vorgabe] #wird nicht wirklich gebraucht

@everywhere const dt	= 1/(T-1)
@everywhere const dx	= 1/(max(m,n) -1)

# Zuordnung Samplenummer zu Zeitpunkt 

# I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)

#s      = inits(rot_circle_ex)[:,:,1:5]
I_vorgabe   = init_vorgabe(rot_circle_ex, 2*m,2*n, T_vorgabe)
# s      = readtaxi()[:,:, 1:5:end]

s			= I_vorgabe[:,:,vorgabe_used_indices] 

velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= 0* ones( m-1, n, T-1 )
end

pygui(true)

# include("verfahren.jl") 
# 
# @everywhere rootdir = "../out/new/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)_mgtol$(mg_tol)/"
# run(`mkdir -p $rootdir/src`)
# run(`sh -c "cp *jl $rootdir/src"`)
# run(`sh -c "git log -1 > $rootdir/this_git_commit"`) #thr
# 
# steps=1
# 
# @time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)
# #@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad_altnormalization(s, u, v, steps)
# 
# # Differenz zur Vorgabe
# diff_vorgabe	= zeros( size(I_vorgabe) )
# for t in 1:T_vorgabe
# 	j					= vorgabe_frames[t]
# 	diff_vorgabe[:,:,t]	= I_vorgabe[:,:,t] - I[:,:,j]
# end
# 
# echo("L2( I-I_vorgabe )", L2norm(diff_vorgabe))
# 
# _="fertig"
