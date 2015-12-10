include("matrizen_zellgrenzen.jl")
const Lx	= generateLu(m,n,dx)
const Ly	= generateLv(m,n,dx)
const LxLU, LyLU	= factorize(Lx), factorize(Ly)

const Cx_zg	= generate_differentiation_interfx(m,n,dx)
const Cy_zg	= generate_differentiation_interfy(m,n,dx)

# Interpolationsmatrizen
const P_zgx = abs(Cx_zg)*dx/2
const P_zgy = abs(Cy_zg)*dx/2

function H1_norm_nobeta_interf(u, v)
	ret = 0
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

project_divfree	== true && const S		= generate_stokes(m, n, dx, Cx_zg, Cy_zg, Lx, Ly)
project_divfree	== true && const SLU	= factorize(S)

# die Multigridverfahren funktionieren sehr schlecht mit der Stokesmatrix
# include("pyamg.jl")
# const S_ml			= construct_mgsolv(S)

function solve_stokes(grd_u_J, grd_v_J)
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

@everywhere @inline function grad_nobeta_interf_time_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx_zg, Cy_zg, P_zgx, P_zgy, LxLU, LyLU, t)
	# p interpolieren
	p_zgx			= P_zgx * reshape(p[:,:,t], m*n)
	p_zgy			= P_zgy * reshape(p[:,:,t], m*n)

	pI_x			= Cx_zg * reshape(I[:,:,t], m*n) .* p_zgx
	pI_y			= Cy_zg * reshape(I[:,:,t], m*n) .* p_zgy

	phi_x			= reshape( solverf(LxLU, -pI_x), m, n-1 )
	phi_y			= reshape( solverf(LyLU, -pI_y), m-1, n )

	# thr makro einfuegen
	#grd_u_J[:,:,t], grd_v_J[:,:,t]	 = solve_stokes( phi_x + alpha*u[:,:,t] , phi_y + alpha*v[:,:,t] )

	grd_u_J[:,:,t]	= phi_x + alpha*u[:,:,t] 
	grd_v_J[:,:,t]	= phi_y + alpha*v[:,:,t] 
end

function grad_J_nobeta_interf_ser(I, p, u, v)
	echo( "================Calculate gradient no time reg serial $m x $n" )
	grd_u_J	= zeros( m, n-1, T-1 )
	grd_v_J	= zeros( m-1, n, T-1 )
	for t= 1:T-1
		grad_nobeta_interf_time_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx_zg, Cy_zg, P_zgx, P_zgy, LxLU, LyLU, t)
	end
	return grd_u_J, grd_v_J
end

function grad_J_nobeta_interf_par(I, p, u, v)
	echo( "================Calculate gradient no time reg parallel $m x $n" )
	grd_u_J	= SharedArray(Float64, (m, n-1, T-1), init= S -> S[localindexes(S)] = 0.0)
	grd_v_J	= SharedArray(Float64, (m-1, n, T-1), init= S -> S[localindexes(S)] = 0.0)
	@sync @parallel for t= 1:T-1
		grad_nobeta_interf_time_slice!(grd_u_J, grd_v_J, I, p, u, v, Cx_zg, Cy_zg, P_zgx, P_zgy, LxLU, LyLU, t)
	end
	return grd_u_J, grd_v_J
end

