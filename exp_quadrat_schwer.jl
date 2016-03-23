# armijo_bas			= 0.5
# armijo_sig			= 0
# armijo_maxtry		= 40

armijo_bas			= 0.9
armijo_sig			= 0
armijo_maxtry		= 80
grad_bound			= 1e-8

@everywhere const alpha	= 0.0000001
@everywhere const beta	= 0.0000001

# maxsteps 			= 10
maxsteps 			= 100000

save_every			= 0

time_regularization	= false  # geht nicht mit velocities_at interfaces

# velocities_at		= "interfaces"
velocities_at		= "centers"

transport_parallel				= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll
@everywhere interpolate_w_time	= false

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

@everywhere with_cfl_check = true

# Zeitregularisierung funktioniert nur mit Flussdiskretisierung an Zellmittelpunkten
# diese Zeile ist zu Sicherheit, damit man nichts falsch einstellt
# velocities_at		= ~time_regularization ? velocities_at : "centers"

include("view.jl")

@everywhere const m					= 60 
@everywhere const n					= 60

include("beispiele.jl")

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples				= 2

@everywhere const auslassen				= 19 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 3 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.
# @everywhere const zwischen_ausgelassen	= 12 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

# die Anzahl zwischen den Referenzframes zu generierenden Frames. 
@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

# ...................... T, alle Frames/Zeitpunkte, also T-1 Zeitschritte von einem Frame auf den naechsten
@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1

@everywhere const T_vorgabe				= auswahl_vorgabe(auslassen, n_samples)[end]

@everywhere const dt	= 1/(T-1)
@everywhere const dx	= 1/(max(m,n) -1)

I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)
# I_vorgabe   = init_vorgabe(rot_circle_ex, 2*m,2*n, T_vorgabe)[m+1:2*m, n+1:2*n, :]
# I_vorgabe   = init_vorgabe(transl_circle, m, n, T_vorgabe)
# I_vorgabe   = init_vorgabe(deform_circle, m, n, T_vorgabe)
# I_vorgabe   = init_vorgabe(__rot_circle_ex, m,n, T_vorgabe)
# I_vorgabe   = readtaxi_alt()[3:192, 3:end-2, 1:T_vorgabe] # die dlm-dateien wurden am Rand mit Nullen aufgefuellt

@everywhere rfac=0.3
srand(1) #random seed setzen, damit fuer verschiedene Durchlaeufe vergleichbarere Ergebnisse 
randerr	= randn(size(I_vorgabe))
I_vorgabe+= rfac*randerr

s			= I_vorgabe[:,:,auswahl_vorgabe(auslassen, n_samples)] 

# Aufloesungsunabhaengig, ohne Vorgabe
# s[:,:,1] = zeros(m,n)
# s[:,:,2] = zeros(m,n)
# for i = 1:m
# 	for j = 1:n
# 		s[i,j,1]	= ((j-1)*dx >= 0.3) & ((j-1)*dx<=0.6) & (((i-1)*dx) >=0.3 )& (((i-1)*dx)<=0.6)
# 		s[i,j,2]	= ((j-1)*dx >= 0.5) & ((j-1)*dx<=0.8) & (((i-1)*dx) >=0.5 )& (((i-1)*dx)<=0.8)
# 	end
# end


velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= 0* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= 0* ones( m-1, n, T-1 )
end

include("verfahren.jl") 

@everywhere rootdir = "../out/rausch/$(rfac)/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)/"
make_output_dir(rootdir)

echo=_echolog
@time I, u, v, p, L2_errs, H1_errs, J, H1_J_ws, steps = verfahren_grad(s, u, v, 1, 1.0, grad_bound)
save_endergebnis(rootdir)

# # Differenz zur Vorgabe
vorgabe_fehler	= diff_vorgabe(I_vorgabe, I, auslassen, zwischen_ausgelassen)
echo("L2( I-I_vorgabe )", L2norm(vorgabe_fehler))
echo("linf( I-I_vorgabe )", l_inf(vorgabe_fehler))
echo("PNSR( I-I_vorgabe )", psnr(vorgabe_fehler))
echo("Gradnorm", H1_J_ws[end])
echo("L2( I-I_vorgabe )", L2norm(vorgabe_fehler))
echo("linf( I-I_vorgabe )", l_inf(vorgabe_fehler))
for l=1:T_vorgabe
	echo("vorgabefehler", l, "psnr", psnr(vorgabe_fehler[:,:,l]), "L2", vorgabe_fehler[:,:,l][:]'*B*vorgabe_fehler[:,:,l][:], "Linf", l_inf(vorgabe_fehler[:,:,l]))
end

demo_table("test", "test")
save_demo_taxi([(".png", 100),(".eps", 1200)])

# save_all()

_="fertig"
