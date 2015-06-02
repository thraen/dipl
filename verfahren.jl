include("misc.jl")
include("transport.jl")

function grad_J_nobeta(I, p, u, v)
	println( "================calculate gradient $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	#thr parallelisieren
	for t= 1:T-1
		pI_x__			= Cx*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n)
		pI_y__			= Cy*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n)
		phi_x__			= poissolv_(pI_x__, m, n)
		phi_y__			= poissolv_(pI_y__, m, n)

		grd_u_J[:,:,t]	= phi_x__+alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y__+alpha*v[:,:,t]
	end
	return grd_u_J, grd_v_J
end

@everywhere function parallel_dim(I, p, uv, Cxy, L, WaveOp, gmres)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv				= L*  reshape(uv[:,:,t], n*m)
		pI_xy			= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)
		rhs[:,:,t]	= (beta-alpha)* Luv + pI_xy 
		if (t==1) || (t==T-1)
			rhs[:,:,t] /= 2
		end
	end
	grd_uv_J, conv_hist	= gmres(WaveOp, reshape(rhs, (T-1)*n*m))
	return grd_uv_J
end

function grad_J_beta_parallel(I, p, u, v)
	println( "================calculate gradient $m x $n" )
	grd_u_J = @spawn parallel_dim(I, p, u, Cx, L, WaveOp, gmres)
	grd_v_J	= @spawn parallel_dim(I, p, v, Cy, L, WaveOp, gmres)

	return reshape(fetch(grd_u_J), m, n, T-1), reshape(fetch(grd_v_J), m, n, T-1)
end

function grad_J_beta(I, p, u, v)
	println( "================calculate gradient $m x $n" )
	rhs_x	= zeros( m, n, T-1 )
	rhs_y	= zeros( m, n, T-1 )

	for t= 1:T-1
		# thr STOP! reshape ordnet Spalten vor Zeilen
		Lu				= L* reshape(u[:,:,t], n*m)
		Lv				= L* reshape(v[:,:,t], n*m)
		pI_x			= Cx*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)
		pI_y			= Cy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)
		# Cx, Cy setzen schon 0-Randbedingungen auf dem Rand, 0-setzen der Raender nicht noetig
		rhs_x[:,:,t]	= (beta-alpha)* Lu + pI_x 
		rhs_y[:,:,t]	= (beta-alpha)* Lv + pI_y
		if (t==1) || (t==T-1)
			rhs_x[:,:,t] /= 2
			rhs_y[:,:,t] /= 2
		end
	end

	#grd_u_J	= WaveOpLU \ reshape(rhs_x, (T-1)*n*m)
	#grd_v_J	= WaveOpLU \ reshape(rhs_y, (T-1)*n*m)

	grd_u_J, conv_hist	= gmres(WaveOp, reshape(rhs_x, (T-1)*n*m))
	grd_v_J, conv_hist	= gmres(WaveOp, reshape(rhs_y, (T-1)*n*m))

	return reshape(grd_u_J, m, n, T-1), reshape(grd_v_J, m, n, T-1)
end

function grad_J_alt(I, p, u, v, alpha)
	println( "================calculate gradient $m x $n" )
	grd_u_J	= zeros( m, n, T-1 )
	grd_v_J	= zeros( m, n, T-1 )
	#thr parallelisieren
	for t= 1:T-1
		# p[2:m-1,2:n-1] vielleicht lieber im Voraus fuer alle t berechnen, damit die p[:,:,t] nicht umkopiert werden muessen

		# marcel's Variante
		pI_x_			= reshape(Cx*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		pI_y_			= reshape(Cy*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		phi_x_			= poissolv( pI_x_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		phi_y_			= poissolv( pI_y_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		grd_u_J[:,:,t]	= phi_x_+alpha*u[:,:,t] 
		grd_v_J[:,:,t]	= phi_y_+alpha*v[:,:,t]

		# meine Variante
		#pI_x			= p[2:m-1,2:n-1,t].*central_diff_x( I[:,:,t] )
		#pI_y			= p[2:m-1,2:n-1,t].*central_diff_y( I[:,:,t] )
		#phi_x			= poissolv( pI_x, zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		#phi_y			= poissolv( pI_y, zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) )
		#grd_u_J[:,:,t]	= phi_x+alpha*u[:,:,t] 
		#grd_v_J[:,:,t]	= phi_y+alpha*v[:,:,t]
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
	println("=============START $n x $m x Aufloesung $T ($n_samples samples $n_zwischensamples zwischsamples)")
	s0			= s[:,:,1]
	norm_s		= L2norm(s)

	I			= zeros(m, n, T)
	p			= zeros(m, n, T)

	L2_err, _	= sample_err(I,s,norm_s)
	H1_err		= H1_norm( u, v )
	J			= L2_err/2 + alpha*H1_err/2

	steps 		= 0
	while steps < maxsteps
		I			= transport(s0, u, v, T-1)
		p			= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
		u, v		= next_w!(I, p, u, v, alpha)

		L2_err, _	= sample_err(I,s,norm_s)
		H1_err		= H1_norm( u, v )
		J			= L2_err/2 + alpha*H1_err/2

		echo()
		echo(steps, "L2errors",  L2_err, "H1_errors", H1_err)
		echo("J", J)
		echo()
		steps+=1
	end

	return I, u, v, p, L2_err, H1_err, J, steps
end

function verfahren_grad(s, u, v)
	println("=============START $n x $m x Aufloesung $T ($n_samples samples $n_zwischensamples zwischsamples)")
	s0			= s[:,:,1]
	norm_s		= L2norm(s)
	echo("norm_s", norm_s)

	H1_err		= H1_norm( u, v )

	I			= transport(s0, u, v, T-1)

	p			= ruecktransport( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)

	echo("initial L2_err", L2_err)

	grd_u_J, grd_v_J	= grad_J(I, p, u, v)
	H1_J_w				= H1_norm(grd_u_J, grd_v_J)
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

			H1_err_next			= H1_norm(u_next, v_next)
			H1_J_w				= H1_norm(grd_u_J, grd_v_J)

			I_next				= transport( s0, u_next, v_next, T-1 )
			L2_err_next, _		= sample_err(I_next,s,norm_s)

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

				p					= ruecktransport(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)

				tic()
				grd_u_J, grd_v_J	= grad_J(I, p, u, v)
				toc()

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

	return I, u, v, p, L2_err, H1_err, J, H1_J_w, steps
end

