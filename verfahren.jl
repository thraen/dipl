@everywhere m					=  40
@everywhere n					=  40

@everywhere n_samples			=   5
@everywhere n_zwischensamples	=  25   # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
sample_times		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]

armijo_lim			= 0.1
armijo_bas			= 0.5
armijo_sig			= 0.5

const dx			= 1
const dt			= 0.95#(600/6-1/6) /100 #0.99833

@everywhere alpha	= 0.01

maxsteps 			= 1000

#@everywhere rootdir = "$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)/"
@everywhere rootdir = "/tmp/out/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)/"

include("marcel_matrizen.jl")
include("misc.jl")
include("transport.jl")
include("poisson.jl")
include("view.jl")
include("beispiele.jl")

# marcel Finite-Differenzen-Matrizen
_, Cx, Cy, Dx, Dy = generateMatrices3(n, 1)

function grad_J(I, p, u, v)
	println( "================calculate gradient $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	for t= 1:T-1
		# p[2:m-1,2:n-1] vielleicht lieber im Voraus fuer alle t berechnen, damit die p[:,:,t] nicht umkopiert werden muessen

		# marcel's Variante
		##pI_x_			= Cx*sI_h.values.*sp_h.values;
		#I_x_			= reshape(Cx*reshape(I[:,:,t] , n*m) , m, n )
		#I_y_			= reshape(Cy*reshape(I[:,:,t] , n*m) , m, n )

		#pI_x_			= I_x_.*p[:,:,t]
		#pI_y_			= I_y_.*p[:,:,t]

		pI_x_			= reshape(Cx*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		pI_y_			= reshape(Cy*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		phi_x_			= poissolv( LU, pI_x_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		phi_y_			= poissolv( LU, pI_y_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		grd_u_J[:,:,t]	= phi_x_+alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y_+alpha*v[:,:,t]

		# meine Variante
		#pI_x			= p[2:m-1,2:n-1,t].*central_diff_x( I[:,:,t] )
		#pI_y			= p[2:m-1,2:n-1,t].*central_diff_y( I[:,:,t] )
		# what? meiner meinung nach sollte hier -pI_x stehen und im naechsten +alpha*u
		#phi_x			= poissolv( LU, pI_x, zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		#phi_y			= poissolv( LU, pI_y, zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		# ..und hier +alpha*u
		#grd_u_J[:,:,t]	= phi_x+alpha*u[:,:,t] 
		#grd_v_J[:,:,t]	= phi_y+alpha*v[:,:,t]
	end
	return grd_u_J, grd_v_J
end

s		= inits(quadrat)
#s		= inits(rot_circle)
s0		= s[:,:,1]

norm_s	= L2norm(s)

echo("norm_s", norm_s)

u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

println("=============START $n x $m x Aufloesung $T ($n_samples samples $n_zwischensamples zwischsamples)")

# matrix fuer poissongleichung wird nur einmal vorbestimmt und spaeter wiederverwertet
L			= generate_laplace(m, n) 
LU			= factorize(L)

H1_err		= _H1_norm( u, v )

I			= transport(s0, u, v, T-1)

# what! hier sollte ,..,-u, -v uebergeben werden # hier ist es eh noch egal, wenn u0 = v0 = 0
p, L2_err	= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
L2_err, _	= sample_err_L2(I,s,norm_s)

echo("initial L2_err", L2_err)

grd_u_J, grd_v_J = grad_J( I, p, u, v )
#echo("max grd_J", maximum((grd_u_J)), maximum((grd_v_J)), maximum( max( (grd_u_J), (grd_v_J)) ) )
#echo("min grd_J", minimum((grd_u_J)), minimum((grd_v_J)), minimum( max( (grd_u_J), (grd_v_J)) ) )

J	= L2_err/2 + alpha*H1_err/2
J0	= J
H0	= H1_err
L0	= L2_err

# Armijo-Schrittweite
armijo_exp = 0
steps = 0
while steps < maxsteps
	while (armijo_exp < 40)
		t 					= armijo_bas^armijo_exp

		echo()
		echo("step", steps, armijo_exp,"test armijo step length ", t)
		echo("step", steps, armijo_exp,"test armijo step length ", t)
		echo()

		u_next				= u - t*grd_u_J
		v_next				= v - t*grd_v_J

		H1_err_next			= _H1_norm(u_next, v_next)
		H1_J_w				= _H1_norm(grd_u_J, grd_v_J)

		I_next				= transport( s0, u_next, v_next, T-1 )
		L2_err_next, _		= sample_err_L2(I_next,s,norm_s)

		J_next 				= L2_err_next/2 + alpha*H1_err_next/2

		#echo( "L2errors", L2_err, L2_err_next, tmp, t2 )

		#echo("max u\t", maximum(abs(u)), "max u_next", maximum(abs(u_next)))
		#echo("max v\t", maximum(abs(v)), "max v_next", maximum(abs(v_next)))
		#echo("max I\t", maximum(abs(I)), "max I_next", maximum(abs(I_next)))

		echo("L2errors",  L2_err, L2_err_next, L2_err-L2_err_next)
		echo("H1_errors", H1_err, H1_err_next, H1_err-H1_err_next)
		echo("alpha H1_errors", alpha*H1_err, alpha*H1_err_next, alpha*(H1_err-H1_err_next))
		echo("J        ", J, J_next,J-J_next)
		echo("H1_J_w", H1_J_w)
		echo()

		#wtf?
		if (J_next < J) 
		#if J_next < J-armijo_sig*t*H1_J_w
			I					= I_next
			u					= u_next
			v					= v_next

			H1_err				= H1_err_next
			L2_err				= L2_err_next

			p, _				= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )

			grd_u_J, grd_v_J	= grad_J( I, p, u, v )

			J					= L2_err/2 + alpha*H1_err/2

			armijo_exp = 0
			echo("\n****** NEW GRADIENT *****")
			echo("max grd_J", maximum((grd_u_J)), maximum((grd_v_J)), maximum( max( (grd_u_J), (grd_v_J)) ) )
			echo("min grd_J", minimum((grd_u_J)), minimum((grd_v_J)), minimum( max( (grd_u_J), (grd_v_J)) ) )
			break 
		end
		
		armijo_exp += 1
	end

	steps +=1
end

#save_images_(s, "s")


