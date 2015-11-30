using HDF5, JLD
include("echo.jl")
include("misc.jl")

@everywhere function solve_lin_multig(A,b)
	# call pyamg 
	return A[:solve](b, tol=1e-3)
end

# das gmres aus KrylowMethods ist leider schlechter (langsamer) als das von IterativeSolvers
# aber das Gauss-Seidel-Verfahren ist besser.
#@everywhere using IterativeSolvers
#@everywhere using KrylovMethods   #thr das spaeter noch mal probieren. das Paket ist suboptimal

@everywhere function solve_lin_gmres(A,b)
	#x, conv_hist	= gmres(A, reshape(b, (T-1)*n*m), restart=5)
	x, conv_hist	= gmres(A, b, restart=5)
	return x
end

@everywhere function solve_lin_elim(A,b)
	return A\b
end

#solverf = solve_lin_gmres
#solverf = solve_lin_multig
solverf = solve_lin_elim

function solve_stokes(grd_u_J, grd_v_J)
	ndofu	= m*(n-1)
	ndofv	= (m-1)*n

	rhs		= [ Lx*grd_u_J[:] ; Ly*grd_v_J[:] ; zeros(m*m) ]
	rhs[end]= 1
	
	res = SLU\rhs

	u_proj	= res[1              : ndofu]
	v_proj	= res[ndofu+1        : ndofu+ndofv]
	p		= res[ndofu+ndofv+1 : end]

	return reshape(u_proj, m, n-1), reshape(v_proj, m-1, n)
end

function grad_J_nobeta_interf(I, p, u, v)
	echo( "================Calculate gradient beta =0 $m x $n" )
	grd_u_J	= zeros( m, n-1, T-1 )
	grd_v_J	= zeros( m-1, n, T-1 )

	for t= 1:T-1
		# p interpolieren
		p_zgx			= P_zgx * reshape(p[:,:,t], m*n)
		p_zgy			= P_zgy * reshape(p[:,:,t], m*n)
	
		pI_x			= Cx_zg * reshape(I[:,:,t], m*n) .* p_zgx
		pI_y			= Cy_zg * reshape(I[:,:,t], m*n) .* p_zgy

		#phi_x			= reshape( solverf(Lx, -pI_x), m, n-1 )
		#phi_y			= reshape( solverf(Ly, -pI_y), m-1, n )
		phi_x			= reshape( solverf(LxLU, -pI_x), m, n-1 )
		phi_y			= reshape( solverf(LyLU, -pI_y), m-1, n )

		
		grd_u_J[:,:,t], grd_v_J[:,:,t]	 = solve_stokes( phi_x + alpha*u[:,:,t] , phi_y + alpha*v[:,:,t] )

		#grd_u_J[:,:,t]	= phi_x + alpha*u[:,:,t] 
		#grd_v_J[:,:,t]	= phi_y + alpha*v[:,:,t] 
	end
	return grd_u_J, grd_v_J
end

function grad_J_nobeta(I, p, u, v)
	echo( "================Calculate gradient beta =0 $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	for t= 1:T-1
		#thr hier auch *dt^2? schau noch mal nach! es geht besser mit. vielleicht ist hier auch ein dx^2 irgendwo zuviel
		pI_x			= Cx*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n) #*dx^2 #*dt^2 #wtf?
		pI_y			= Cy*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n) #*dx^2 #*dt^2
		phi_x			= poissolv_(pI_x, m, n)
		phi_y			= poissolv_(pI_y, m, n)

		grd_u_J[:,:,t]	= phi_x + alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y + alpha*v[:,:,t] 
	end
	return grd_u_J, grd_v_J
end

@everywhere function solve_ellip_beta(I, p, uv, Cxy, L, ellOp)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv			= L*  reshape(uv[:,:,t], n*m)
		pI_xy		= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)

		#thr, hier muss eigentlich eigentlich mit dt^2 multipliziert werden
		#thr!!
		rhs[:,:,t]	= ((beta-alpha)* Luv + pI_xy) * dt^2

		if (t==1) || (t==T-1)
			rhs[:,:,t] /= 2
		end
	end

	zuv = solverf( ellOp, reshape(rhs, (T-1)*n*m) )
	return zuv
end

function grad_J_beta(I, p, u, v) 
	echo( "================Calculate gradient $m x $n" )
	zu = solve_ellip_beta(I, p, u, Cx, L, ml)
	zv = solve_ellip_beta(I, p, v, Cy, L, ml)
	grd_u_J, grd_v_J	= reshape(zu, m, n, T-1)+beta*u, reshape(zv, m, n, T-1)+beta*v
	return grd_u_J, grd_v_J
end

function grad_J_beta_parallel(I, p, u, v)
	echo( "================Calculate gradient $m x $n" )
	zu = @spawn solve_ellip_beta(I, p, u, Cx, L, ellOp)
	zv = @spawn solve_ellip_beta(I, p, v, Cy, L, ellOp)
	return reshape(fetch(zu), m, n, T-1) +beta*u, reshape(fetch(zv), m, n, T-1)+beta*v
end

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

		echo()
		echo(steps, "L2errors",  L2_err, "H1_errors", H1_err)
		echo("J", J)
		echo()
		steps+=1
	end

	return I, u, v, p, L2_err, H1_err, J, steps
end

function verfahren_grad(s, u, v, steps=1)
	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dt=$dt, alpha=$alpha, beta=$beta")
	s0			= s[:,:,1]
	norm_s		= L2norm(s)
	echo("norm_s", norm_s)

	H1_err		= H1_norm_w( u, v )

	@time I			= transport(s0, u, v, T-1)
	@time p			= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)

	echo("initial L2_err", L2_err)

	@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

	@show H1_J_w			= H1_norm_grd(grd_u_J, grd_v_J)

	@show J	= L2_err/2 + alpha*H1_err/2
	@show J0	= J
	@show H0	= H1_err
	@show L0	= L2_err


	# Armijo-Schrittweite
	armijo_exp	= 0
	while steps < maxsteps  &&  armijo_exp < 80  &&  H1_J_w > 1e-8 
		while (armijo_exp < 80)
			t 					= armijo_bas^armijo_exp

			echo()
			echo("step", steps, armijo_exp,"test armijo step length ", t)
			echo("step", steps, armijo_exp,"test armijo step length ", t)
			echo()

			u_next				= u - t*grd_u_J
			v_next				= v - t*grd_v_J

			H1_err_next			= H1_norm_w(u_next, v_next)
			# thr, das sollte besser beim gradientenupdate stehen
			H1_J_w				= H1_norm_grd(grd_u_J, grd_v_J)

			@time I_next		= transport( s0, u_next, v_next, T-1 )
			L2_err_next, _		= sample_err(I_next,s,norm_s)

			J_next 				= L2_err_next/2 + alpha*H1_err_next/2

			#echo( "L2errors", L2_err, L2_err_next, tmp, t2 )

			#echo("max u\t", maximum(abs(u)), "max u_next", maximum(abs(u_next)))
			#echo("max v\t", maximum(abs(v)), "max v_next", maximum(abs(v_next)))
			#echo("max I\t", maximum(abs(I)), "max I_next", maximum(abs(I_next)))

			echo("L2errors",  		L2_err, L2_err_next, L2_err-L2_err_next)
			echo("alpha H1_errors", alpha*H1_err, alpha*H1_err_next, alpha*(H1_err-H1_err_next))
			echo("J        ", J, J_next,J-J_next)
			echo("H1_J_w", H1_J_w)
			echo()

			#if (J_next < J) 
			if J_next < J - armijo_sig * t *H1_J_w
				I					= I_next
				u					= u_next
				v					= v_next

				H1_err				= H1_err_next
				L2_err				= L2_err_next

				@time p					= ruecktransport(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)
				@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

				J					= L2_err/2 + alpha*H1_err/2

				armijo_exp = 0
				echo("\n****** NEW GRADIENT *****")
				#echo("max grd_J", maximum((grd_u_J)), maximum((grd_v_J)), maximum( max( (grd_u_J), (grd_v_J)) ) )
				#echo("min grd_J", minimum((grd_u_J)), minimum((grd_v_J)), minimum( min( (grd_u_J), (grd_v_J)) ) )
				break 
			end
			
			armijo_exp += 1
		end

		if (save_every > 0) && (steps % save_every == 0)
			try
				info("\n zwischenspeichern\n $steps")
				save("$(rootdir)zwischenergebnis_$steps.jld", 
					 	"dx", dx,
						"dt", dt,
					 	"alpha", alpha,
						"beta", beta,
					 	"s", s,
						"I", I, 
						"p", p,
						"u", u,
						"v", v, 
						"grd_u_J", grd_u_J, 
						"grd_v_J", grd_v_J)	
			catch e
				warn("ZWISCHENERGEBNIS KONNTE NICHT GESPEICHERT WERDEN!", e)
			end
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

