armijo_bas			= 0.9
armijo_sig			= 0.0
armijo_maxtry		= 80

@everywhere const alpha	= 0.01
@everywhere const beta	= 0.01

# maxsteps 			= 2
maxsteps 			= 100000

save_every			= 0

@everywhere time_regularization	= false 
@everywhere velocities_at		= "centers"
# @everywhere velocities_at		= "interfaces"
@everywhere interpolate_w_time	= true

transport_parallel	= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll
grad_parallel		= true # betrifft nur die Verfahren ohne Zeitregularisierung
project_divfree		= false # betrifft nur velocities_at = "interfaces"
poisson_solver		= "lufact" #fur gegebene Probleme am besten. Eigentlich Cholesky-Faktorisierung fuer die interfaces und LU-Faktorisierung fuer center
stokes_solver		= "lufact"#fur gegebene Probleme am besten
timereg_solver		= "multig"#fur gegebene Probleme am besten

@everywhere const mg_tol = 1e-1 

@everywhere with_cfl_check = true

include("view.jl")

@everywhere const m					= 100
@everywhere const n					= 100

include("beispiele.jl")

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples				= 2

@everywhere const auslassen				= 10 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 5 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1
@everywhere const T_vorgabe				= auswahl_vorgabe(auslassen, n_samples)[end]
@everywhere const dt	= 1/(T-1)
@everywhere const dx	= 1/(max(m,n) -1)

I_vorgabe   	= init_vorgabe(rot_circle_ex, 2*m,2*n, T_vorgabe)[m+1:2*m, n+1:2*n, :]
# I_vorgabe_lat   = init_vorgabe(rot_circle_ex_alt, 2*m,2*n, T_vorgabe)[m+1:2*m, n+1:2*n, :]

# @everywhere rfac=0.2
# srand(1) #random seed setzen, damit fuer verschiedene Durchlaeufe vergleichbarere Ergebnisse 
# randerr	= randn(size(I_vorgabe))
# I_vorgabe+= rfac*randerr

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

@everywhere rootdir = "../out/demo/exp_rot_disc/$(velocities_at)/time_reg_$(time_regularization)/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)/"

# make_output_dir(rootdir)
# echo=_echolog
@time I, u, v, p, L2_errs, H1_errs, Js, H1_J_ws, steps = verfahren_grad(s, u, v, 1, 1.0)
# save_endergebnis(rootdir)

# using JLD
# @load "$(rootdir)res.jld"

# echo("==============")
# echo("Gradnorm", H1_J_ws[end])
# echo("L2err", L2_errs[end])
# echo("unweightened space reg error", H1_norm_w_noweight_space(u,v))
# echo("unweightened time reg error", H1_norm_w_noweight_time(u,v))
# 
# vorgabe_fehler	= diff_vorgabe(I_vorgabe, I, auslassen, zwischen_ausgelassen)
# echo("L2( I-I_vorgabe )", L2norm(vorgabe_fehler))
# echo("linf( I-I_vorgabe )", l_inf(vorgabe_fehler))
# for l=1:T_vorgabe
# 	echo("vorgabefehler", l, "psnr", psnr(vorgabe_fehler[:,:,l]), "L2", vorgabe_fehler[:,:,l][:]'*B*vorgabe_fehler[:,:,l][:], "Linf", l_inf(vorgabe_fehler[:,:,l]))
# end
# 
# demo_table("demoRotDisc", "demo_rot_disc")
# save_demo_rot_disc([(".png", 100),(".eps", 1200)])
# save_demo_rot_disc_rauschen([(".png", 100),(".eps", 1200)])
# save_displacement(rootdir, ".png", 100)
# save_displacement(rootdir, ".eps", 1200)

_="fertig"
