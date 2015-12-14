include("matrizen_zellgrenzen.jl")
@everywhere const Lx	= generateLu(m,n,dx)
@everywhere const Ly	= generateLv(m,n,dx)

const Cx_zg	= generate_differentiation_interfx(m,n,dx)
const Cy_zg	= generate_differentiation_interfy(m,n,dx)

# Interpolationsmatrizen
const P_zgx = abs(Cx_zg)*dx/2
const P_zgy = abs(Cy_zg)*dx/2

poisson_solver == "lufact" && begin
	# die Cholesky-Zerlegung kann nicht zu anderen Prozessen kopiert werden. 
	# daher muss man die Zerlegungen auf den anderen Prozessen generieren
	@everywhere const LxFact, LyFact	= factorize(Lx), factorize(Ly)
	@everywhere function solve_poisson_x(b) 
		return LxFact\b
	end
	@everywhere function solve_poisson_y(b) 
		return LyFact\b
	end
end

poisson_solver == "gmres" && begin
	#@everywhere using KrylovMethods   #thr das spaeter noch mal probieren. das Paket ist suboptimal
	@everywhere using IterativeSolvers
	@everywhere function solve_poisson_x(b)
		x, conv_hist	= gmres(L, b, restart=5)
		return x
	end
	@everywhere function solve_poisson_y(b)
		x, conv_hist	= gmres(Ly, b, restart=5)
		return x
	end
end

poisson_solver == "multig" && begin
	error("Zwischenzellgeschwindigkeiten koennen nich mit Multigridverfahren benutzt werden")
end

function H1_norm(u, v)
	ret = 0
	@show m,n
	@show size(u)
	@show size(v)
	for t=1:T-1
		u_ = reshape(u[:,:,t], (n-1)*m)
		v_ = reshape(v[:,:,t], n*(m-1))
		ret_  = (u_'*Lx*u_) 
		ret_ += (v_'*Ly*v_) 
		if t==1 || t==T
			ret_ /= 2
		end
		ret += ret_
	end
	# thr! dx steckt auch in L? stimmt das so?
	return -ret[1] *dx*dx *dt
end

H1_norm_w	= H1_norm
H1_norm_grd	= H1_norm

project_divfree	== true && const S		= generate_stokes(m, n, dx, Cx_zg, Cy_zg, Lx, Ly)
project_divfree	== true && const SLU	= factorize(S)

# die Multigridverfahren funktionieren sehr schlecht mit der Stokesmatrix
# include("pyamg.jl")
# const S_ml			= construct_mgsolv(S)

@everywhere function solve_stokes(grd_u_J, grd_v_J)
	ndofu	= m*(n-1)
	ndofv	= (m-1)*n

	rhs		= [ Lx*grd_u_J[:] ; Ly*grd_v_J[:] ; zeros(m*m) ]
	rhs[end]= 1
	
	res = SLU\rhs

	# geht leider nicht so schnell, wie erhofft. LU ist fuer Taxi noch bedeutend schneller
	#res	= solve_lin_multig( S_ml, rhs )

	u_proj	= res[1              : ndofu]
	v_proj	= res[ndofu+1        : ndofu+ndofv]
	p		= res[ndofu+ndofv+1 : end]

	return reshape(u_proj, m, n-1), reshape(v_proj, m-1, n)
end

@everywhere macro do_stokes_when_defined()
	return :(grd_u_J[:,:,t], grd_v_J[:,:,t]	 = solve_stokes( phi_x + alpha*u[:,:,t] , phi_y + alpha*v[:,:,t] ))
end

@everywhere @inline function grad_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx_zg, Cy_zg, P_zgx, P_zgy, t)
	# p interpolieren
	p_zgx			= P_zgx * reshape(p[:,:,t], m*n)
	p_zgy			= P_zgy * reshape(p[:,:,t], m*n)

	pI_x			= Cx_zg * reshape(I[:,:,t], m*n) .* p_zgx
	pI_y			= Cy_zg * reshape(I[:,:,t], m*n) .* p_zgy

	phi_x			= reshape( solve_poisson_x(-pI_x), m, n-1 )
	phi_y			= reshape( solve_poisson_y(-pI_y), m-1, n )

	#grd_u_J[:,:,t], grd_v_J[:,:,t]	 = solve_stokes( phi_x + alpha*u[:,:,t] , phi_y + alpha*v[:,:,t] )
	#@do_stokes_when_defined

	grd_u_J[:,:,t]	= phi_x + alpha*u[:,:,t] 
	grd_v_J[:,:,t]	= phi_y + alpha*v[:,:,t] 
end

function grad_J(I, p, u, v)
	echo( "================Calculate gradient no time reg $m x $n  parallel=$(grad_parallel)" )
	grd_u_J	= @init_grad( m, n-1, T-1 )
	grd_v_J	= @init_grad( m-1, n, T-1 )

	@show typeof(grd_u_J)

	@do_par_when_defined for t= 1:T-1
		grad_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx_zg, Cy_zg, P_zgx, P_zgy, t)
	end
	return grd_u_J, grd_v_J
end

