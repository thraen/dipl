using HDF5, JLD
include("echo.jl")

include("transport.jl")

include("matrizen.jl")
const B		= generate_L2(m, n, dx)
# const B		= generateB(m, dx)

include("misc.jl")
L2norm		= function(s) return Xnorm(s, B) end
sample_err	= sample_err_L2

time_regularization == true		&& include("grad_time_reg.jl")  
~time_regularization && velocities_at == "centers"		&& include("grad_centers.jl") 	
~time_regularization && velocities_at == "interfaces"	&& include("grad_interfaces.jl")

function next_w!(I, p, u, v, alpha)
	for t= 1:T-1
		pI_x_		= reshape(Cx*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		pI_y_		= reshape(Cy*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		u[:,:,t]	= poissolv( -pI_x_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) ) /alpha
		v[:,:,t]	= poissolv( -pI_y_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) ) /alpha
	end
	return u, v
end

function verfahren_direkt(s, u, v)
	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dy=$dy, alpha=$alpha, beta=$beta")
	s0			= s[:,:,1]
	norm_s		= L2norm(s) #thr. gewicht?

	I			= zeros(m, n, T)
	p			= zeros(m, n, T)

	L2_err, _	= sample_err(I,s,norm_s)
	H1_err		= H1_norm_w( u, v )
	J			= L2_err/2 + alpha*H1_err/2

	steps 		= 1
	while steps < maxsteps
		I			= transport(s0, u, v, T-1)
		p			= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
		u, v		= next_w!(I, p, u, v, alpha)

		L2_err, _	= sample_err(I,s,norm_s)
		H1_err		= H1_norm_w( u, v )
		J			= L2_err/2 + alpha*H1_err/2

		echo(steps, "L2errors",  L2_err, "H1_errors", H1_err, 
			,"\nJ", J, "\n")
		steps+=1
	end

	return I, u, v, p, L2_err, H1_err, J, steps
end

function verfahren_grad(s, u, v, steps=1)
	s0			= s[:,:,1]
	norm_s		= L2norm(s)

	H1_err		= H1_norm_w( u, v )

	@time I			= transport(s0, u, v, T-1)
	@time p			= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)

	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dt=$dt, alpha=$alpha, beta=$beta",
		 "\nnorm_s", norm_s,
		 "\ninitial L2_err", L2_err)

	@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

	@show H1_J_w			= H1_norm_grd(grd_u_J, grd_v_J)

	@show J	= L2_err/2 + alpha*H1_err/2
	@show J0	= J
	@show H0	= H1_err
	@show L0	= L2_err

	# Armijo-Schrittweite
	armijo_exp	= 0

	while steps < maxsteps  &&  armijo_exp < 40  &&  H1_J_w > 1e-8 
		while (armijo_exp < 40)
			t 					= armijo_bas^armijo_exp

			u_next				= u - t*grd_u_J
			v_next				= v - t*grd_v_J

			H1_err_next			= H1_norm_w(u_next, v_next)
			# thr, das sollte besser beim gradientenupdate stehen
			H1_J_w				= H1_norm_grd(grd_u_J, grd_v_J)

			@time I_next		= transport( s0, u_next, v_next, T-1 )
			L2_err_next, _		= sample_err(I_next,s,norm_s)

			J_next 				= L2_err_next/2 + alpha*H1_err_next/2

			#echo("max u\t", maximum(abs(u)), "max u_next", maximum(abs(u_next)))
			#echo("max v\t", maximum(abs(v)), "max v_next", maximum(abs(v_next)))
			#echo("max I\t", maximum(abs(I)), "max I_next", maximum(abs(I_next)))

			echo("\nstep", steps, armijo_exp,"test armijo step length ", t, 
				 "\nL2errors",  		L2_err, L2_err_next, L2_err-L2_err_next, 
				 "\nalpha H1_errors", alpha*H1_err, alpha*H1_err_next, alpha*(H1_err-H1_err_next),
				 "\nJ        ", J, J_next,J-J_next,"\n")

			if J_next < J - armijo_sig * t *H1_J_w
				I					= I_next
				u					= u_next
				v					= v_next

				H1_err				= H1_err_next
				L2_err				= L2_err_next

				@time p					= ruecktransport(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)
				@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

				J					= L2_err/2 + alpha*H1_err/2 # thr! das stimmt doch nicht mit Zeitregularisierung

				armijo_exp = 0
				echo("\n****** NEW GRADIENT *****\n", 
					 "max abs grd_u_J", maximum(abs(grd_u_J)),
					 "max abs grd_v_J", maximum(abs(grd_v_J)), "\n",
					 "H1_J_w", H1_J_w, "\n")
				break 
			end
			
			armijo_exp += 1
		end

		if (save_every > 0) && (steps % save_every == 0)
			save_jld(steps, dx, dt, alpha, beta, s, I, p, u, v, grd_u_J, grd_v_J)
		end

		steps +=1
	end

	return I, u, v, p, L2_err, H1_err, J, H1_J_w, steps
end

function verfahren_grad_goldstein(s, u, v, steps=1)
	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dt=$dt, alpha=$alpha, beta=$beta")
	s0			= s[:,:,1]
	norm_s		= L2norm(s)
	echo("norm_s", norm_s)

	H1_err		= H1_norm_w( u, v )

	@time I		= transport(s0, u, v, T-1)
	@time p		= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)

	echo("initial L2_err", L2_err)

	@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)
	H1_J_w					= H1_norm_grd(grd_u_J, grd_v_J)

	J	= L2_err/2 + alpha*H1_err/2

	u_next = 0
	v_next = 0
	H1_err_next			= 0
	L2_err_next			= 0
	J_next 				= 0
	I_next				= 0

	function try_J(t)
		u_next				= u - t*grd_u_J
		v_next				= v - t*grd_v_J

		H1_err_next			= H1_norm_w(u_next, v_next)

		@time I_next		= transport( s0, u_next, v_next, T-1 )
		L2_err_next, _		= sample_err(I_next,s,norm_s)

		J_next 				= L2_err_next/2 + alpha*H1_err_next/2
		# echo( "try t", t, L2_err_next, H1_err_next, J_next)

		return J_next
	end
	
	function update()
		echo( "update", L2_err_next, H1_err_next, J_next)
		I					= I_next
		u					= u_next
		v					= v_next

		H1_err				= H1_err_next
		L2_err				= L2_err_next

		@time p					= ruecktransport(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)
		@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

		H1_J_w				= H1_norm_grd(grd_u_J, grd_v_J)
		J					= L2_err/2 + alpha*H1_err/2
		steps += 1
	end

	function upper_J(t)
		return J - sig * t * H1_J_w
	end

	function lower_J(t)
		return J - (1-sig) * t * H1_J_w
	end

	tu	= 0
	to	= 1
	while steps < maxsteps
		echo( "STEP ", steps, J, H1_J_w )
		@show try_J(to)  lower_J(to)
		while try_J(to) < lower_J(to)
			echo("calibrating window", tu, to)
			tu=to
			to*=2
		end

		if J_next <= upper_J(to)
			echo("passt")
			update()
			break
		end

		t=(tu+to)/2
		J_next = try_J(t)
		while ~( lower_J(t) <= J_next <= upper_J(t) )
			t=(tu+to)/2
			try_J(t)
			echo(tu, t, to)
			echo(lower_J(t), J_next, upper_J(t))
			if 		J_next < lower_J(t)
				tu	= t
			end
			if  	J_next > upper_J(t)
				to	= t
			end
		end
		update()
	end

	return I, u, v, p, L2_err, H1_err, J, H1_J_w, steps
end

function verfahren_grad_try_par(s, u, v, steps=1)
	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dt=$dt, alpha=$alpha, beta=$beta")
	norm_s		= L2norm(s)
	echo("norm_s", norm_s)

	H1_err		= H1_norm_w( u, v )

	@time I		= transport(s[:,:,1], u, v, T-1)
	@time p		= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)

	echo("initial L2_err", L2_err)

	@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)
	H1_J_w					= H1_norm_grd(grd_u_J, grd_v_J)

	J	= L2_err/2 + alpha*H1_err/2

	function try_J(u, v, t, grd_u_J, grd_v_J, s, norm_s)
		u_next				= u - t*grd_u_J
		v_next				= v - t*grd_v_J

		H1_err_next			= H1_norm_w(u_next, v_next)

		I_next				= transport( s[:,:,1], u_next, v_next, T-1 )
		L2_err_next, _		= sample_err(I_next,s,norm_s)

		J_next 				= L2_err_next/2 + alpha*H1_err_next/2
		# echo( "try t", t, L2_err_next, H1_err_next, J_next)

		return I_next, L2_err_next, H1_err_next, J_next
	end
	
	function update()
		echo( "update", L2_err_next, H1_err_next, J_next)
		I					= I_next
		u					= u_next
		v					= v_next

		H1_err				= H1_err_next
		L2_err				= L2_err_next

		@time p					= ruecktransport(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)
		@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

		H1_J_w				= H1_norm_grd(grd_u_J, grd_v_J)
		J					= L2_err/2 + alpha*H1_err/2
		steps += 1
	end

	armijo_exp = 0
	while steps < maxsteps
		#for proc = 1:nprocs
			#rc	jkkk
		#end
	

		steps +=1
	end


	return I, u, v, p, L2_err, H1_err, J, H1_J_w, steps
end


