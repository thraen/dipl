@everywhere const ellOp, GradNormOp, CostNormOp, Sreg, Treg	= generate_ellip_beta(m, n, T, dt, dx, alpha, beta)
const L			= generate_laplace(m, n, dx)
const Cx, Cy	= generate_differentiation_central(m, n, dx) 

timereg_solver == "lufact" && begin
	info("Vorsicht: LU-Zerlegung fuer den Zeitregularisierungsoperator braucht sehr viel RAM!")
	info("LU-Zerlegung fuer Zeitregularisierung funktioniert nicht parallel")
	const ellOpLU	= factorize(ellOp)
	function solve_timereg_x(b)
		return ellOpLU\b
	end
	function solve_timereg_y(b)
		return ellOpLU\b
	end
end

timereg_solver == "gmres" && begin
	@everywhere using IterativeSolvers
	#@everywhere using KrylovMethods   #thr das spaeter noch mal probieren. das Paket ist suboptimal
	@everywhere function solve_timereg_x(b)
		x, conv_hist	= gmres(ellOp, b, restart=5) 
		return x
	end
	@everywhere function solve_timereg_y(b)
		x, conv_hist	= gmres(ellOp, b, restart=5) 
		return x
	end
end

timereg_solver == "multig" && begin
	include("pyamg.jl")
	# fuer die Parallelausfuehrung in den beiden Dimensionen ist es wichtig, dass die 
	# zwei verschiedene Multigrid-Objekte angelegt werden, weil der Zustand zur Laufzeit jeweils 
	# darin steckt. Wenn dann verschiedene Prozesse denselben MG-Solver benutzen geht es schief.

	@everywhere const ellOp_mg_x		= construct_mgsolv(ellOp)
	@everywhere const ellOp_mg_y		= construct_mgsolv(ellOp)
	@everywhere function solve_timereg_x(b)
		return ellOp_mg_x[:solve](b, tol=mg_tol, accel="cg")
	end
	@everywhere function solve_timereg_y(b)
		return ellOp_mg_y[:solve](b, tol=mg_tol, accel="cg")
	end
end

function H1_norm_w_timereg(u,v)
	u_	= reshape(u, m*n*(T-1))
	v_	= reshape(v, m*n*(T-1))
	# alpha und beta stecken in CostNormOp
	ret	= dx*dx* (u_'*CostNormOp*u_ + v_'*CostNormOp*v_)	
	return  ret[1] 
end

function H1_norm_grd_timereg(u,v)
	u_	= reshape(u, m*n*(T-1))
	v_	= reshape(v, m*n*(T-1))
	ret	= dx*dx* (u_'*GradNormOp*u_ + v_'*GradNormOp*v_)
	return  ret[1]
end

function H1_norm_w_noweight_space_centers(u,v)
	return dx*dx* (u[:]'*Sreg*u[:] + v[:]'*Sreg*v[:])
end

function H1_norm_w_noweight_time_centers(u,v)
	return dx*dx* (u[:]'*Treg*u[:] + v[:]'*Treg*v[:])
end

H1_norm_w   = H1_norm_w_timereg
H1_norm_grd = H1_norm_grd_timereg

H1_norm_w_noweight_space	= H1_norm_w_noweight_space_centers
H1_norm_w_noweight_time		= H1_norm_w_noweight_time_centers

@everywhere function constr_rhs_beta(I, p, uv, Cxy, L)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv			= L*  reshape(uv[:,:,t], n*m)
		pI_xy		= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)
		#thr, warum hier muss hier eigentlich eigentlich mit dt^2 multipliziert werden?
		rhs[:,:,t]	= ((beta-alpha)* Luv + pI_xy) * dt^2
		if (t==1) || (t==T-1)
			#@show t
			#@show maximum(abs(I[:,:,t]))
			#@show maximum(abs(p[:,:,t]))
			#@show maximum(abs(pI_xy))
			rhs[:,:,t] /= 2
		end
	end

	return reshape(rhs, (T-1)*n*m)
end

# 	pygui(true)
@everywhere function grad_J_beta_dim_x(I, p, uv, Cxy, L)
	rhs	= constr_rhs_beta(I, p, uv, Cxy, L)
# 	@show(reshape(rhs, m, n, T-1)[:,:,5])
# 	surf(reshape(rhs, m, n, T-1)[:,:,5], cstride=1, rstride=1)
	#imshow(reshape(rhs, m, n, T-1)[:,:,5])
	#savefig("~/tr_cent.png")
	zuv = solve_timereg_x( rhs )
	return zuv
end

@everywhere function grad_J_beta_dim_y(I, p, uv, Cxy, L)
	rhs	= constr_rhs_beta(I, p, uv, Cxy, L)
	zuv = solve_timereg_y( rhs )
	return zuv
end

function grad_J(I, p, u, v)
	echo( "================Calculate gradient with time regularization $m x $n parallel=$grad_parallel" )
	zu = @spawn grad_J_beta_dim_x( I, p, u, Cx, L )
	zv = @spawn grad_J_beta_dim_y( I, p, v, Cy, L )
	return reshape(fetch(zu), m, n, T-1) +beta*u, reshape(fetch(zv), m, n, T-1)+beta*v
end

