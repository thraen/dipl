const L			= generate_laplace(m, n, dx)
const Cx, Cy	= generate_differentiation_central(m, n, dx) 

poisson_solver == "lufact" && begin
	const LU		= factorize(L)
	@everywhere function solve_poisson(LU, b)
		return LU\b
	end
end

poisson_solver == "gmres" && begin
	#@everywhere using IterativeSolvers
	#@everywhere using KrylovMethods   #thr das spaeter noch mal probieren. das Paket ist suboptimal
	@everywhere function solve_poisson(b)
		#x, conv_hist	= gmres(L, reshape(b, (T-1)*n*m), restart=5)
		x, conv_hist	= gmres(L, b, restart=5)
		return x
	end
end

poisson_solver == "multig" && begin
	grad_parallel && error( "Verfahren ohne Zeitregularisierung und Multigrid geht nicht parallel" ) #thr vielleicht doch
	include("pyamg.jl")
	const L_mg		= construct_mgsolv(L)
	@everywhere function solve_poisson(b)
		return L[:solve](b, tol=mg_tol, accel="cg")
	end
end

# Integral genaehert durch Integral einer Stueckweise linearen Interpolation
function H1_norm(u, v)
	ret = 0
	for t=1:T-1
		u_ = reshape(u[:,:,t], n*m)
		v_ = reshape(v[:,:,t], n*m)
		ret_  = (u_'*L*u_) 
		ret_ += (v_'*L*v_) 
		if t==1 || t==T
			ret_ /= 2
		end
		ret += ret_
	end
	# thr! dx steckt auch in L? stimmt das so?
	return dt* dx*dx*ret[1]
end

# Integration durch Integral einer stueckweise konstanten Naeherung
#function H1_norm(u, v)
	#return dx*dx* dt* ( sum( central_diff_x( u ).^2 + sum( central_diff_y( u ).^2 ) ) + sum( central_diff_x( v ).^2 + sum( central_diff_y( v ).^2 ) ) )
#end

H1_norm_w	= function (u,v) return alpha* H1_norm(u,v) end
H1_norm_grd	= H1_norm

@everywhere @inline function grad_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx, Cy, LU, t)
	# die 0-Randbedingung steckt in der Multiplikation mit den Differentiationsmatrizen. Die Matrix setzt den schon Rand 0!
	pI_x			= Cx*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n)
	pI_y			= Cy*reshape(I[:,:,t], n*m).* reshape(p[:,:,t], m*n)

	phi_x			= solve_poisson(LU, pI_x)
	phi_y			= solve_poisson(LU, pI_y)

# 	if t==1
# 		imshow(reshape( pI_y, m, n ))
# 		imshow(reshape( pI_x, m, n ))
# 		imshow(reshape( phi_y, m, n ))
# 		imshow(reshape( phi_x, m, n ))
# 	end

	grd_u_J[:,:,t]	= reshape(phi_x, m, n) + alpha*u[:,:,t] 
	grd_v_J[:,:,t]	= reshape(phi_y, m, n) + alpha*v[:,:,t] 
end

function grad_J(I, p, u, v)
	echo( "================Calculate gradient, no time reg $m x $n  parallel=$(grad_parallel)")
	grd_u_J	= @init_grad(m,n,T)
	grd_v_J	= @init_grad(m,n,T)

	#convert(SharedArray{Float64}, I)
	#convert(SharedArray{Float64}, p)
	#convert(SharedArray{Float64}, u)
	#convert(SharedArray{Float64}, v)
	#=
		@show macroexpand(:(
		@do_par_when_defined for t= 1:T-1
			grad_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx, Cy, LU, t)
		end
		))
	=#

	@do_par_when_defined for t= 1:T-1
		grad_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx, Cy, LU, t)
	end

	#clf()
	#quiver(grd_u_J[:,:,1], grd_v_J[:,:,1])
	#sleep(1)

	return grd_u_J, grd_v_J
end 

