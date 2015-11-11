using HDF5, JLD
include("echo.jl")
include("misc.jl")
#include("transport.jl")
include("transport_neu.jl")
#include("transport_dxfalsch_divfalsch.jl")

function grad_J_nobeta(I, p, u, v)
	echo( "================Calculate gradient beta =0 $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	for t= 1:T-1
		#thr hier auch *dt^2? schau noch mal nach! es geht besser mit. vielleicht ist hier auch ein dx^2 irgendwo zuviel
		pI_x__			= Cx*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n) #*dx^2 #*dt^2 #wtf?
		pI_y__			= Cy*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n) #*dx^2 #*dt^2
		phi_x__			= poissolv_(pI_x__, m, n)
		phi_y__			= poissolv_(pI_y__, m, n)

		grd_u_J[:,:,t]	= phi_x__+alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y__+alpha*v[:,:,t] 
	end
	return grd_u_J, grd_v_J
end

@everywhere function parallel_dim(I, p, uv, Cxy, L, ellOp)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv			= L*  reshape(uv[:,:,t], n*m)
		pI_xy		= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)

		#thr, hier muss eigentlich eigentlich mit dt^2 multipliziert werden
		#thr!!
		rhs[:,:,t]	= ((beta-alpha)* Luv + pI_xy) #* dt^2

		if (t==1) || (t==T-1)
			rhs[:,:,t] /= 2
		end
	end

	zuv, conv_hist	= gmres(ellOp, reshape(rhs, (T-1)*n*m), restart=5)
	return zuv
end

@everywhere function parallel_dim_pymultig(I, p, uv, Cxy, L, ml)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv			= L*  reshape(uv[:,:,t], n*m)
		pI_xy		= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)

		#thr, hier muss eigentlich eigentlich mit dt^2 multipliziert werden
		#thr!!
		rhs[:,:,t]	= ((beta-alpha)* Luv + pI_xy) #* dt^2

		if (t==1) || (t==T-1)
			rhs[:,:,t] /= 2
		end
	end

	# call pyamg via python
	zuv			= ml[:solve](reshape(rhs, (T-1)*n*m), tol=1e-3)
	return zuv
end

#function grad_J_zellgrenzen_multig(I, p, u, v)
	#echo( "================Calculate gradient $m x $n" )
	#zu	= parallel_dim_pymultig(I, p, u, Cx, L
#end

function grad_J_beta_multig(I, p, u, v)
	echo( "================Calculate gradient $m x $n" )
	zu = parallel_dim_pymultig(I, p, u, Cx, L, ml)
	zv = parallel_dim_pymultig(I, p, v, Cy, L, ml)
	grd_u_J, grd_v_J	= reshape(zu, m, n, T-1)+beta*u, reshape(zv, m, n, T-1)+beta*v
	return grd_u_J, grd_v_J
end

function grad_J_beta_multig_parallel(I, p, u, v)
	echo( "================Calculate gradient $m x $n" )
	zu = @spawn parallel_dim_pymultig(I, p, u, Cx, L, ml)
	zv = @spawn parallel_dim_pymultig(I, p, v, Cy, L, ml)
	return reshape(fetch(zu), m, n, T-1) +beta*u, reshape(fetch(zv), m, n, T-1)+beta*v
end

function grad_J_beta_parallel(I, p, u, v)
	echo( "================Calculate gradient $m x $n" )
	zu = @spawn parallel_dim(I, p, u, Cx, L, ellOp)
	zv = @spawn parallel_dim(I, p, v, Cy, L, ellOp)
	return reshape(fetch(zu), m, n, T-1) +beta*u, reshape(fetch(zv), m, n, T-1)+beta*v
end

function grad_J_beta(I, p, u, v)
	echo( "================Calculate gradient $m x $n" )
	zu = parallel_dim(I, p, u, Cx, L, ellOp)
	zv = parallel_dim(I, p, v, Cy, L, ellOp)

	grd_u_J, grd_v_J	= reshape(zu, m, n, T-1)+beta*u, reshape(zv, m, n, T-1)+beta*v
	return grd_u_J, grd_v_J
end

function grad_J_alt(I, p, u, v, alpha)
	echo( "================Calculate gradient $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	for t= 1:T-1
		# p[2:m-1,2:n-1] vielleicht lieber im Voraus fuer alle t berechnen, damit die p[:,:,t] nicht umkopiert werden muessen

		pI_x_			= reshape(Cx*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		pI_y_			= reshape(Cy*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		phi_x_			= poissolv( pI_x_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		phi_y_			= poissolv( pI_y_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		grd_u_J[:,:,t]	= phi_x_+alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y_+alpha*v[:,:,t]
	end
	return grd_u_J, grd_v_J
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
	H1_J_w					= H1_norm_grd(grd_u_J, grd_v_J)

	J	= L2_err/2 + alpha*H1_err/2
	J0	= J
	H0	= H1_err
	L0	= L2_err

	# Armijo-Schrittweite
	armijo_exp = 0
	while steps < maxsteps
		while (armijo_exp < 40)
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
				echo("max grd_J", maximum((grd_u_J)), maximum((grd_v_J)), maximum( max( (grd_u_J), (grd_v_J)) ) )
				echo("min grd_J", minimum((grd_u_J)), minimum((grd_v_J)), minimum( min( (grd_u_J), (grd_v_J)) ) )
				break 
			end
			
			armijo_exp += 1
		end

		if (save_every > 0) && (steps % save_every == 0)
			try
				info("\n zwischenspeichern\n")
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

