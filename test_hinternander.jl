@everywhere const m					= 60
@everywhere const n					= 60

# insgesamt. diese variable wird nicht zur Verwendung der globalen Matrizen verwendet.
@everywhere const _n_samples		= 10
# fuer die Hintereinanderausfuehrung, jeweils von einem Sample zum nächsten
@everywhere const n_samples			= 2

armijo_bas			= 0.5
armijo_sig			= 0.0

@everywhere const dt	= 0.01
@everywhere const dx	= 0.01

@everywhere const alpha	= 0.00009 #*0.001 #warum ist das nicht dasselbe, wie die norm noch mal mit alpha zu multiplizieren? siehe CostNormOp --> thr ahh, wegen der ruecksubstitution nach ellipt gleichung. aber warum funktioniert es so gut, wenn man die norm noch mal mit alpha multipliziert. teste das auch mal ohne zeitreg
@everywhere const beta	= 0.00009

# maxsteps 			= 1
maxsteps 			= 100000

save_every			= 0

time_regularization	= false  # geht nicht mit velocities_at interfaces

# velocities_at		= "interfaces"
velocities_at		= "centers"

transport_parallel	= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll

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
include("beispiele.jl")

@everywhere const auslassen				= 2 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 4 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.
# die Anzahl zwischen den Referenzframes zu generierenden Frames. 
@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen
# Zuordnung Samplenummer zu Zeitpunkt 
@everywhere const sample_times			= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples-1 ]
@everywhere const _sample_times			= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:_n_samples-1 ]

@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1
@everywhere const _T					= (_n_samples-1)*(n_zwischensamples+1) +1

@everywhere const vorgabe_used_indices	= (1:(auslassen+1):(auslassen+1)*_n_samples) 
@everywhere const T_vorgabe				= vorgabe_used_indices[end]
@everywhere const samples_to_vorgabe	= [(k, vorgabe_used_indices[k]) for k in 1:_n_samples]

@everywhere const vorgabe_frames		= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
# @everywhere const vorgabe_to_frames		= [(k,vorgabe_frames[k]) for k in 1:T_vorgabe] #wird nicht wirklich gebraucht

I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe) 	###aaaahhrg
s			= I_vorgabe[:,:,vorgabe_used_indices] 		###aaaahhrg
# s			= inits(quadrat)

velocities_at == "centers" && begin
	u		= 0* ones( m, n, _T-1 )
	v		= 0* ones( m, n, _T-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, _T-1 )
	v		= 0* ones( m-1, n, _T-1 )
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

I	= zeros(m,n,_T)


@time for i in 1:length(_sample_times)-1
	@show _sample_times[i]
	from_f	= _sample_times[i][2]
	till_f	= _sample_times[i+1][2]-1
	# u geht nur bis T-1
	u_		= copy(u[:,:,from_f:till_f])
	v_		= copy(v[:,:,from_f:till_f])
	@show from_f till_f
	@show size(u_)


	#samples
	from_s	= _sample_times[i][1]
	till_s	= from_s+1
	@show from_s till_s
	s_		= copy(s[:,:,from_s:till_s])
	
	# I geht bis T, der letzte Frame wird durch das entsprechende sample ersetzt, siehe oben
	@time I_, u_, v_, p_, L2_err_, H1_err_, J_, H1_J_w_, steps_ = verfahren_grad(s_, u_, v_, steps)

	# zusammenkopieren fuers endergebnis
	@show size(I_)
	@show size(I_[:,:,1:end-1])
	@show from_f:till_f
	I[:,:,from_f:till_f]	= I_[:,:,1:end-1]
	#thr das letzte sample muss noch nach I kopiert werden!!!
	#hae!! das ist doch falsch! moment
	I[:,:,till_f+1]			= s_[:,:,end]

	# !!!!! das optische Flussfeld stimmt dann so nicht !!!!!!
	u[:,:,from_f:till_f]		= copy(u_[:,:,1:end])
	v[:,:,from_f:till_f]		= copy(v_[:,:,1:end])

	println( '\n')
end


# Differenz zur Vorgabe
diff_vorgabe	= zeros( size(I_vorgabe) )
for t in 1:T_vorgabe
	@show j				= vorgabe_frames[t]
	diff_vorgabe[:,:,t]	= I_vorgabe[:,:,t] - I[:,:,j]
end

_="fertig"

# nochmal mit restarts


