#const ellOp, GradNormOp, CostNormOp	= generate_ellip_beta(n, T, dt, dx, alpha, beta)
#const ellOpLU	= factorize(ellOp)

#include("pyamg.jl")
#const ellOp_ml		= construct_mgsolv(ellOp)

#@everywhere const ellOp, GradNormOp, CostNormOp	= generate_ellip_beta(n, T, dt, dx, alpha, beta)
#@everywhere const ellOp_ml1		= construct_mgsolv(ellOp)
#@everywhere const ellOp_ml2		= construct_mgsolv(ellOp)


function H1_norm_beta_w(u,v)
	u_	= reshape(u, m*n*(T-1))
	v_	= reshape(v, m*n*(T-1))
	ret	= dx*dx* (u_'*CostNormOp*u_ + v_'*CostNormOp*v_)	
	return  ret[1]
end

function H1_norm_beta_grd(u,v)
	u_	= reshape(u, m*n*(T-1))
	v_	= reshape(v, m*n*(T-1))
	ret	= dx*dx* (u_'*GradNormOp*u_ + v_'*GradNormOp*v_)
	return  ret[1]
end

@everywhere function solve_lin_multig(A,b)
	# call pyamg 
	return A[:solve](b, tol=mg_tol, accel="cg")
end

function solve_ellip_beta(b)
	#return ellOp \ b	
	#x, conv_hist	= gmres(ellOp, b, restart=5)
	return ellOp_ml[:solve](b, tol=mg_tol, accel="cg")
end

# das ist leider notwendig, da die python-handles anders nicht auf die worker-prozesse kopiert werden koennen
@everywhere function solve_ellip_beta1(b)
	return ellOp_ml1[:solve](b, tol=mg_tol, accel="cg")
end

@everywhere function solve_ellip_beta2(b)
	return ellOp_ml2[:solve](b, tol=mg_tol, accel="cg")
end

@everywhere function constr_rhs_beta(I, p, uv, Cxy, L)
	rhs	= zeros( m, n, T-1 )
	for t= 1:T-1
		Luv			= L*  reshape(uv[:,:,t], n*m)
		pI_xy		= Cxy*reshape(I[:,:,t], n*m) .* reshape(p[:,:,t], n*m)
		#thr, warum hier muss hier eigentlich eigentlich mit dt^2 multipliziert werden?
		rhs[:,:,t]	= ((beta-alpha)* Luv + pI_xy) * dt^2
		if (t==1) || (t==T-1)
			rhs[:,:,t] /= 2
		end
	end

	return reshape(rhs, (T-1)*n*m)
end

@everywhere function _grad_J_beta_dim(I, p, uv, Cxy, L)
	@time rhs = constr_rhs_beta(I, p, uv, Cxy, L)
	@time zuv = solve_ellip_beta( rhs )
	return zuv
end

@everywhere function grad_J_beta_dim(I, p, uv, Cxy, L, ellOp)
	@time rhs = constr_rhs_beta(I, p, uv, Cxy, L)

	#thr closure testen.
	#@time zuv = solve_ellip_beta( rhs )
	@time zuv = solverf( ellOp,  reshape(rhs, (T-1)*n*m)  )
	return zuv
end

function grad_J_beta(I, p, u, v) 
	echo( "================Calculate gradient serial $m x $n" )
	#@time zu = _grad_J_beta_dim(I, p, u, Cx, L)
	#@time zv = _grad_J_beta_dim(I, p, v, Cy, L)

	@time zu = grad_J_beta_dim(I, p, u, Cx, L, ellOp_ml)
	@time zv = grad_J_beta_dim(I, p, v, Cy, L, ellOp_ml)

	grd_u_J = reshape(zu, m, n, T-1)+beta*u
	grd_v_J	= reshape(zv, m, n, T-1)+beta*v
	return grd_u_J, grd_v_J
end

@everywhere function grad_J_beta_dim1(I, p, uv, Cxy, L)
	rhs	= constr_rhs_beta(I, p, uv, Cxy, L)
	zuv = solve_ellip_beta1( rhs )
	return zuv
end

@everywhere function grad_J_beta_dim2(I, p, uv, Cxy, L)
	rhs	= constr_rhs_beta(I, p, uv, Cxy, L)
	zuv = solve_ellip_beta2( rhs )
	return zuv
end

function grad_J_beta_parallel(I, p, u, v)
	echo( "================Calculate gradient parallel $m x $n" )
	#zu = @spawn grad_J_beta_dim(I, p, u, Cx, L, solve_ellip_beta1)
	#zv = @spawn grad_J_beta_dim(I, p, v, Cy, L, solve_ellip_beta2)
	
	zu = @spawn grad_J_beta_dim1( I, p, u, Cx, L )
	zv = @spawn grad_J_beta_dim2( I, p, v, Cy, L )

	#rhsx	= constr_rhs_beta(I, p, u, Cx, L)
	#rhsy	= constr_rhs_beta(I, p, u, Cy, L)
	#zu		= @spawn solve_ellip_beta1( reshape(rhsx, (T-1)*n*m)  )
	#zv		= @spawn solve_ellip_beta2( reshape(rhsy, (T-1)*n*m)  )

	#zu = remotecall(2, grad_J_beta_dim1, I, p, u, Cx, L)
	#zv = remotecall(3, grad_J_beta_dim2, I, p, v, Cy, L)

	#@time rhsx	= reshape( constr_rhs_beta(I, p, u, Cx, L), (T-1)*n*m )
	#@time rhsy	= reshape( constr_rhs_beta(I, p, u, Cy, L), (T-1)*n*m )

	#@time zu		= remotecall(2, solve_ellip_beta1, rhsx )
	#@time zv		= remotecall(3, solve_ellip_beta2, rhsy )

	#@time resu	= fetch(zu)
	#@time resv	= fetch(zv)

	#return reshape(resu, m, n, T-1) +beta*u, reshape(resv, m, n, T-1)+beta*v

	return reshape(fetch(zu), m, n, T-1) +beta*u, reshape(fetch(zv), m, n, T-1)+beta*v
end

