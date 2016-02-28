# armijo_bas			= 0.5
# armijo_sig			= 0
# armijo_maxtry		= 40

armijo_bas			= 0.9
armijo_sig			= 0
armijo_maxtry		= 80

# @everywhere const alpha	= 0.2 
# @everywhere const beta	= 0.2

# @everywhere const alpha	= 0.005
# @everywhere const beta	= 0.005

@everywhere const alpha	= 0.001
@everywhere const beta	= 0.001

# @everywhere const alpha	= 0.0001006
# @everywhere const beta	= 0.0001006
 
# @everywhere const alpha	= 0.0001 #* 0.5556353538255331
# @everywhere const beta	= 0.0001

# @everywhere const alpha	= 0.00005 
# @everywhere const beta	= 0.00005 

# @everywhere const alpha	= 0.0001 * 0.5556353538255331

# @everywhere const alpha	= 0.00005
# @everywhere const beta	= 0.00005

# maxsteps 			= 3
maxsteps 			= 100000

save_every			= 0

time_regularization	= false  # geht nicht mit velocities_at interfaces

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

# @everywhere const m					= 60
# @everywhere const n					= 60
# 
@everywhere const m					= 100 
@everywhere const n					= 100

# @everywhere const m					= 140
# @everywhere const n					= 140

include("beispiele.jl")

# fuer die Konstruktion der Zeitregularisierungsmatrizen muss n_samples >=2 und n_zwischensamples >=3 sein!
@everywhere const n_samples				= 60

@everywhere const auslassen				= 0 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 1 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.
# @everywhere const zwischen_ausgelassen	= 12 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

# die Anzahl zwischen den Referenzframes zu generierenden Frames. 
@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

# ...................... T, alle Frames/Zeitpunkte, also T-1 Zeitschritte von einem Frame auf den naechsten
@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1

@everywhere const T_vorgabe				= auswahl_vorgabe(auslassen, n_samples)[end]

# @everywhere const dt	= 1/(T-1)
# @everywhere const dx	= 1/(max(m,n) -1)

@everywhere const dt	= 1.0
@everywhere const dx	= 1.0

@everywhere rootdir = "../out/transporttest/"


v_=0.5

# v_=0.01
# @everywhere const T						= 1000


I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)
s			= I_vorgabe[:,:,auswahl_vorgabe(auslassen, n_samples)] 
s0			= s[:,:,1]

velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= v_* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= v_* ones( m-1, n, T-1 )
end

include("transport.jl") 
I_center	= transport( s0, u, v, T-1 )


velocities_at		= "interfaces"
velocities_at == "centers" && begin
	u		= 0* ones( m, n, T-1 )
	v		= v_* ones( m, n, T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T-1 )
	v		= v_* ones( m-1, n, T-1 )
end

include("transport.jl") 
I_interf	= transport( s0, u, v, T-1 )

diff	= I_center - I_interf
I_center == I_interf

verr_center	= diff_vorgabe(I_vorgabe, I_center, auslassen, zwischen_ausgelassen)
verr_interf	= diff_vorgabe(I_vorgabe, I_interf, auslassen, zwischen_ausgelassen)

pygui(true)
