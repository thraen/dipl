@everywhere const m					= 60
@everywhere const n					= 60

@everywhere const n_samples_gesamt	= 5

armijo_bas			= 0.5
armijo_sig			= 0.0

# @everywhere const dt	= 0.01
@everywhere const dt	= 0.005
@everywhere const dx	= 0.01

@everywhere const alpha	= 0.00009 #*0.001 #warum ist das nicht dasselbe, wie die norm noch mal mit alpha zu multiplizieren? siehe CostNormOp --> thr ahh, wegen der ruecksubstitution nach ellipt gleichung. aber warum funktioniert es so gut, wenn man die norm noch mal mit alpha multipliziert. teste das auch mal ohne zeitreg
@everywhere const beta	= 0.00009

# maxsteps 			= 1
maxsteps 			= 100000

save_every			= 0

time_regularization	= false  # geht nicht mit velocities_at interfaces

velocities_at		= "centers"

transport_parallel	= false # geht nicht gut, erst ab ca 500x500 Pixel sinnvoll

grad_parallel		= false # betrifft nur die Verfahren ohne Zeitregularisierung

project_divfree		= false # betrifft nur velocities_at = "interfaces"

poisson_solver		= "lufact" #fur gegebene Probleme am besten. Eigentlich Cholesky-Faktorisierung fuer die interfaces und LU-Faktorisierung fuer center

stokes_solver		= "lufact"#fur gegebene Probleme am besten

timereg_solver		= "multig"#fur gegebene Probleme am besten
@everywhere const mg_tol = 1e-1 

@everywhere with_cfl_check = false

include("view.jl")
include("beispiele.jl")

# @everywhere const auslassen				= 2 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
# @everywhere const zwischen_ausgelassen	= 4 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

@everywhere const auslassen				= 4 # die Referenzsamples werden so gewählt, dass aus der Vorgabe werden immer `auslassen` Frames weggelassen werden
@everywhere const zwischen_ausgelassen	= 8 # zwischen zwei ausgelassenen Frames sollen so viele Zwischenframes generiert werden.

@everywhere const n_zwischensamples		= auslassen + (auslassen+1) * zwischen_ausgelassen

@everywhere const n_samples				= 2
@everywhere		  samples_to_frames		= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples-1 ]
@everywhere		  smpls_to_frms_gesamt	= [ (k+1, k*(n_zwischensamples+1)+1) for k in 0:n_samples_gesamt-1 ]

@everywhere const T						= (n_samples-1)*(n_zwischensamples+1) +1
@everywhere const T_gesamt				= (n_samples_gesamt-1)*(n_zwischensamples+1) +1

@everywhere const vorgabe_used_indices	= (1:(auslassen+1):(auslassen+1)*n_samples_gesamt) 
@everywhere const T_vorgabe				= vorgabe_used_indices[end]
@everywhere		  samples_to_vorgabe	= [(k, vorgabe_used_indices[k]) for k in 1:n_samples_gesamt]

@everywhere const vorgabe_frames		= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 

I_vorgabe	= init_vorgabe(char_quadrat, m,n, T_vorgabe)
s			= I_vorgabe[:,:,vorgabe_used_indices]

velocities_at == "centers" && begin
	u		= 0* ones( m, n, T_gesamt-1 )
	v		= 0* ones( m, n, T_gesamt-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T_gesamt-1 )
	v		= 0* ones( m-1, n, T_gesamt-1 )
end

include("verfahren.jl") 

@everywhere rootdir = "../out/hinternander/$(m)_x_$(n)_$(n_samples_gesamt)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/separat/"
run(`mkdir -p $rootdir/src`)
run(`sh -c "cp *jl $rootdir/src"`)
run(`sh -c "git log -1 > $rootdir/this_git_commit"`) #thr

steps=1

I	= zeros(m,n,T_gesamt)

tic()
@time for i in 1:length(smpls_to_frms_gesamt)-1
	#frames
	from_f	= smpls_to_frms_gesamt[i][2]
	till_f	= smpls_to_frms_gesamt[i+1][2]-1
	# u geht nur bis T-1
	u_		= copy(u[:,:,from_f:till_f])
	v_		= copy(v[:,:,from_f:till_f])

	#samples
	from_s	= smpls_to_frms_gesamt[i][1]
	till_s	= from_s+1
	s_		= copy(s[:,:,from_s:till_s])
	
	# I geht bis T, der letzte Frame wird durch das entsprechende sample ersetzt, siehe oben
	I_, u_, v_, p_, L2_err_, H1_err_, J_, H1_J_w_, steps_ = verfahren_grad(s_, u_, v_, steps)

	# zusammenkopieren fuers endergebnis
	I[:,:,from_f:till_f]	= I_[:,:,1:end-1]
	I[:,:,till_f+1]			= s_[:,:,end]

	# !!!!! das optische Flussfeld stimmt dann natuerlich jeweils aufs naechste sample nicht !!!!!!
	u[:,:,from_f:till_f]		= copy(u_[:,:,1:end])
	v[:,:,from_f:till_f]		= copy(v_[:,:,1:end])
end
elapsed_t = toc()
echo("\n\nelapsed time: ", elapsed_t)

# Differenz zur Vorgabe
diff_vorgabe	= zeros( size(I_vorgabe) )
for t in 1:T_vorgabe
	@show j				= vorgabe_frames[t]
	diff_vorgabe[:,:,t]	= I_vorgabe[:,:,t] - I[:,:,j]
end
echo("\n\nL2 Differenz: ", L2norm(diff_vorgabe))
save_all()

#nochmal ohne restart

@everywhere const n_samples				= n_samples_gesamt
				  samples_to_frames		= smpls_to_frms_gesamt
@everywhere const T						= T_gesamt

velocities_at == "centers" && begin
	u		= 0* ones( m, n, T_gesamt-1 )
	v		= 0* ones( m, n, T_gesamt-1 )
end 
velocities_at == "interfaces" && begin
	u		= 0* ones( m, n-1, T_gesamt-1 )
	v		= 0* ones( m-1, n, T_gesamt-1 )
end

include("verfahren.jl") 

@everywhere rootdir = "../out/hinternander/$(m)_x_$(n)_$(n_samples_gesamt)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/gesamt/"

steps=1

tic()
I, u, v, p, L2_err, H1_err, J, H1_J_w, steps = verfahren_grad(s, u, v, steps)
elapsed_t = toc()
echo("\n\nelapsed time: ", elapsed_t)

diff_vorgabe	= zeros( size(I_vorgabe) )
for t in 1:T_vorgabe
	@show j				= vorgabe_frames[t]
	diff_vorgabe[:,:,t]	= I_vorgabe[:,:,t] - I[:,:,j]
end

echo("\n\nL2 Differenz: ", L2norm(diff_vorgabe))

save_all()

_="fertig"
