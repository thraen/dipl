#=
	hier wurde ausprobiert, wieviele Zwischenframes noetig sind 
	fuer das Quadrat: Bei 9 ausgelassenen Vorgabeframes werden 
	sukzessive die zwischen_ausgelassen Frames eingefuegt.
	sprich: zwischen_ausgelassen = 0,1,2,3,..

	s und I_vorgabe wurden nur einmal abgespeichert fuer zwischen_ausgelassen = 0
=#

armijo_bas			= 0.5
armijo_sig			= 0.0
armijo_maxtry		= 40

@everywhere const alpha	= 0.001
@everywhere const beta	= 0.001

# maxsteps 			= 10
maxsteps 			= 100000

save_every			= 0

time_regularization				= false  # geht nicht mit velocities_at interfaces
@everywhere interpolate_w_time	= false

# velocities_at		= "interfaces"
velocities_at		= "centers"

transport_parallel	= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll

# das Verfahren mit Zeitregularisierung parallelisiert 
# automatisch die Dimensionen, wenn mehr als ein Worker existiert

grad_parallel		= true # betrifft nur die Verfahren ohne Zeitregularisierung
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

@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1
@everywhere const T_vorgabe				= auswahl_vorgabe(auslassen, n_samples)[end]
@everywhere const dt	= 1/(T-1)
@everywhere const dx	= 1/(max(m,n) -1)

I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)
s			= I_vorgabe[:,:,auswahl_vorgabe(auslassen, n_samples)] 
velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= 0* ones( m-1, n, T-1 )
end

include("verfahren.jl") 
@everywhere rootdir = "../out/demo/exp_n_zwischen_quadrat/alpha_$(alpha)/$(zwischen_ausgelassen)/"
# make_output_dir(rootdir)

# @time I, u, v, p, L2_errs, H1_errs, Js, H1_J_ws, steps = verfahren_grad(s, u, v, 1, 1.0)
# save_endergebnis(rootdir)

using JLD
@load "$(rootdir)res.jld"

# vorgabe_fehler	= diff_vorgabe(I_vorgabe, I, auslassen, zwischen_ausgelassen)

# save_demo_rot_disc([(".png", 100),(".eps", 1200)])
# nzw_table("Test Anzahl Frames", "dummy")

save_displacement(rootdir, ".png", 100)
save_displacement(rootdir, ".eps", 1200)
pygui(true)
