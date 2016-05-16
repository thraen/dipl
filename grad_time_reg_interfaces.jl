include("matrizen_zellgrenzen.jl")
@everywhere const ellOp_x, GradNormOp_x, CostNormOp_x, Sreg_x, Treg_x	= generate_ellip_beta_u(m, n, T, dt, dx, alpha, beta)
@everywhere const ellOp_y, GradNormOp_y, CostNormOp_y, Sreg_y, Treg_y	= generate_ellip_beta_v(m, n, T, dt, dx, alpha, beta)
@everywhere const Lx	= -generateLu(m,n,dx)
@everywhere const Ly	= -generateLv(m,n,dx)

const Cx_zg	= generate_differentiation_interfx(m,n,dx)
const Cy_zg	= generate_differentiation_interfy(m,n,dx)

const P_zgx = abs(Cx_zg)*dx/2
const P_zgy = abs(Cy_zg)*dx/2

# timereg_solver == "lufact" && begin
# 	info("Vorsicht: LU-Zerlegung fuer den Zeitregularisierungsoperator braucht sehr viel RAM!")
# 	info("LU-Zerlegung fuer Zeitregularisierung funktioniert nicht parallel")
# 	const ellOpLU	= factorize(ellOp)
# 	function solve_timereg_x(b)
# 		return ellOpLU\b
# 	end
# 	function solve_timereg_y(b)
# 		return ellOpLU\b
# 	end
# end

timereg_solver == "gmres" && begin
	@everywhere using IterativeSolvers
	#@everywhere using KrylovMethods   #thr das spaeter noch mal probieren. das Paket ist suboptimal
	@everywhere function solve_timereg_x(b)
		x, conv_hist	= gmres(ellOp_x, b, restart=5) 
		return x
	end
	@everywhere function solve_timereg_y(b)
		x, conv_hist	= gmres(ellOp_y, b, restart=5) 
		return x
	end
end

timereg_solver == "multig" && begin
	include("pyamg.jl")
	# fuer die Parallelausfuehrung in den beiden Dimensionen ist es wichtig, dass die 
	# zwei verschiedene Multigrid-Objekte angelegt werden, weil der Zustand zur Laufzeit jeweils 
	# darin steckt. Wenn dann verschiedene Prozesse denselben MG-Solver benutzen geht es schief.

	@everywhere const ellOp_mg_x		= construct_mgsolv(ellOp_x)
	@everywhere const ellOp_mg_y		= construct_mgsolv(ellOp_y)
	@everywhere function solve_timereg_x(b)
		# return ellOp_mg_x[:solve](b, tol=mg_tol)
		return ellOp_mg_x[:solve](b, tol=mg_tol, accel="cg")
	end
	@everywhere function solve_timereg_y(b)
		#return ellOp_mg_y[:solve](b, tol=mg_tol)
		return ellOp_mg_y[:solve](b, tol=mg_tol, accel="cg")
	end
end

function H1_norm_w_timereg(u,v)
	u_	= reshape(u, m*(n-1)*(T-1))
	v_	= reshape(v, (m-1)*n*(T-1))
	# alpha und beta stecken in CostNormOp
	ret	= dx*dx* (u_'*CostNormOp_x*u_ + v_'*CostNormOp_y*v_)	
	return  ret[1] 
end

function H1_norm_grd_timereg(u,v)
	u_	= reshape(u, m*(n-1)*(T-1))
	v_	= reshape(v, (m-1)*n*(T-1))
	ret	= dx*dx* (u_'*GradNormOp_x*u_ + v_'*GradNormOp_y*v_)
	return  ret[1]
end

#thr die hier koennen nicht bei mehreren durchlaeufen mit versch. parametern umdefiniert werden
# musst du noch was machen!
function H1_norm_w_noweight_space_interf(u,v)
	return dx*dx* (u[:]'*Sreg_x*u[:] + v[:]'*Sreg_y*v[:])
end

function H1_norm_w_noweight_time_interf(u,v)
	return dx*dx* (u[:]'*Treg_x*u[:] + v[:]'*Treg_y*v[:])
end

H1_norm_w   = H1_norm_w_timereg
H1_norm_grd = H1_norm_grd_timereg

H1_norm_w_noweight_space	= H1_norm_w_noweight_space_interf
H1_norm_w_noweight_time		= H1_norm_w_noweight_time_interf

@everywhere function constr_rhs_beta_x(I, p, u, Cx_zg, Lx, P_zgx)
	rhs_x	= zeros( m, n-1, T-1 )
	for t= 1:T-1
		Lu			= Lx*  reshape(u[:,:,t], (n-1)*m)
		p_zgx		= P_zgx * reshape(p[:,:,t], m*n)

		pI_x		= Cx_zg * reshape(I[:,:,t], m*n) .* p_zgx

		rhs_x[:,:,t]	= ((beta-alpha)* Lu + pI_x) * dt^2
		if (t==1) || (t==T-1)
			rhs_x[:,:,t] /= 2
		end
	end

	return reshape(rhs_x, (T-1)*(n-1)*m)
end

@everywhere function constr_rhs_beta_y(I, p, v, Cy_zg, Ly, P_zgy)
	rhs_y	= zeros( m-1, n, T-1 )
	for t= 1:T-1
		Lv			= Ly*  reshape(v[:,:,t], n*(m-1))
		p_zgy		= P_zgy * reshape(p[:,:,t], m*n)

		pI_y		= Cy_zg * reshape(I[:,:,t], m*n) .* p_zgy

		rhs_y[:,:,t]	= ((beta-alpha)* Lv + pI_y) * dt^2
		if (t==1) || (t==T-1)
			rhs_y[:,:,t] /= 2
		end
	end

	return reshape(rhs_y, (T-1)*n*(m-1))
end

@everywhere function grad_J_beta_dim_x(I, p, u, Cx_zg, Lx, P_zgx)
	rhs_x	= constr_rhs_beta_x(I, p, u, Cx_zg, Lx, P_zgx)
	zv		= solve_timereg_x( rhs_x )
	return zv
end

@everywhere function grad_J_beta_dim_y(I, p, v, Cy_zg, Ly, P_zgy)
	rhs_y	= constr_rhs_beta_y(I, p, v, Cy_zg, Ly, P_zgy)
	zu		= solve_timereg_y( rhs_y )
	return zu
end


function grad_J(I, p, u, v)
	echo( "================Calculate gradient with time regularization $m x $n parallel=$grad_parallel" )
	zu = @spawn grad_J_beta_dim_x(I, p, u, Cx_zg, Lx, P_zgx)
	zv = @spawn grad_J_beta_dim_y(I, p, v, Cy_zg, Ly, P_zgy)
	return reshape(fetch(zu), m, n-1, T-1) +beta*u, reshape(fetch(zv), m-1, n, T-1)+beta*v
end

