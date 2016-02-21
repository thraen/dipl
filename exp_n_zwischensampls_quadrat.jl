#=
	hier wurde ausprobiert, wieviele Zwischenframes noetig sind 
	fuer das Quadrat: Bei 9 ausgelassenen Vorgabeframes werden 
	sukzessive die zwischen_ausgelassen Frames eingefuegt.
	sprich: zwischen_ausgelassen = 0,1,2,3,..

	s und I_vorgabe wurden nur einmal abgespeichert fuer zwischen_ausgelassen = 0
=#

# ergebnisse: zwischen_ausgelassen zu L2norm(diff_bild)
# 0 0.002051071429918122 geht nicht, die typische CFL-Kante
# 1 0.037483320364157524 geht tatsaechlich auch. CFL ist mindestens zwischenzeitlich verletzt
# 2 0.014741600705809448 geht
# 3 0.013988767632127173 geht
# 4 0.014011044831468963 geht
# 5 0.01822167291787852
# 6 0.020832244878910472

armijo_bas			= 0.5
armijo_sig			= 0.0

@everywhere const alpha	= 0.0009 
@everywhere const beta	= 0.0009

# maxsteps 			= 10
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

@everywhere with_cfl_check = false

# Zeitregularisierung funktioniert nur mit Flussdiskretisierung an Zellmittelpunkten
# diese Zeile ist zu Sicherheit, damit man nichts falsch einstellt
# velocities_at		= ~time_regularization ? velocities_at : "centers"

include("view.jl")

@everywhere const m					= 60
@everywhere const n					= 60

include("beispiele.jl")

@everywhere const n_samples				= 2

@everywhere const auslassen				= 9
@everywhere const zwischen_ausgelassen	= 0

# zwischen_auslassen soll veraendert werden und alles ab hier copy paste..

@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen
@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1
@everywhere		  samples_to_frames		= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples-1 ]

@everywhere const vorgabe_used_indices	= (1:(auslassen+1):(auslassen+1)*n_samples) 
@everywhere const T_vorgabe				= vorgabe_used_indices[end]

@everywhere		  samples_to_vorgabe	= [(k, vorgabe_used_indices[k]) for k in 1:n_samples]

@everywhere const dt	= 1/(T-1)
@everywhere const dx	= 1/(max(m,n) -1)
I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)
s			= I_vorgabe[:,:,vorgabe_used_indices] 
velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= 0* ones( m-1, n, T-1 )
end
include("verfahren.jl") 
@everywhere rootdir = "../out/exp_n_zwischen_quadrat/$(zwischen_ausgelassen)/"
make_output_dir(rootdir)
steps=1
@time I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)

vorgabe_fehler	= diff_vorgabe(I_vorgabe, I, auslassen, zwischen_ausgelassen)

echo("L2( I-I_vorgabe )", L2norm(vorgabe_fehler))

# @everywhere dpi=1200
# @everywhere isuff=".eps"
# save_images_(s, "s")
# save_surfs_(I_vorgabe, "I_given")
# save_surfs_(I,"I")
# save_surfs_(diff_vorgabe, "diff_vorgabe")


