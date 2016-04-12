
@everywhere function generate_laplace_inner(m, n)
	_m, _n	= m-2, n-2
	Dxx = spdiagm( (ones(_n-1),-2*ones(_n),ones(_n-1)), (-1,0,1), _n, _n) #1D discrete Laplacian in the x-direction ;
	Dyy = spdiagm( (ones(_m-1),-2*ones(_m),ones(_m-1)), (-1,0,1), _m, _m) #1D discrete Laplacian in the y-direction ;
	return -(kron(Dyy, speye(_n)) + kron(speye(_m), Dxx))
end

function generate_differentiation_central(m, n, dx)
	println("generate central differences Matrices for cell centers $m x $n")
	ndiag_x_m	= -[ repmat( [ 0; ones(m-2); 0], n-2) ; zeros(m) ] /(2*dx)
	ndiag_x_p	=  [ zeros(m) ; repmat( [ 0; ones(m-2); 0], n-2) ] /(2*dx)
	Cx	= spdiagm( (ndiag_x_m, ndiag_x_p), (-m, m) )

	ndiag_y_m	= -[ zeros(m-1); repmat( [0; ones(m-2); 0], n-2); zeros(m) ] /(2*dx)
	ndiag_y_p	=  [ zeros(m); repmat( [0; ones(m-2); 0], n-2); zeros(m-1) ] /(2*dx)
	Cy	= spdiagm( (ndiag_y_m, ndiag_y_p), (-1, 1) )
	return Cx, Cy
end

function L2diag(m, n, dx)
	return [2; 6*ones(m-2); 4;   repmat([6; 12*ones(m-2); 6], n-2);     4; 6*ones(m-2); 2]* (dx^2/24)
end

function L2diagM(m, n, dx)
	return repmat([1; 2*ones(m-2); 1], n-1)* (dx^2/24)
end

function L2diag1(m, n, dx)
	return [ ones(m-1); 0;  repmat([2*ones(m-1); 0], n-2); ones(m-1) ]* (dx^2/24)
end

function L2diagpMm1(m, n, dx)
  	return [ repmat([0; 2*ones(m-1)], n-1); 0 ]* (dx^2/24)
end

function generate_L2(m, n, dx)
	diags = (L2diagpMm1(m, n, dx), L2diagM(m, n, dx), L2diag1(m, n, dx), L2diag(m, n, dx), L2diag1(m,n,dx), L2diagM(m,n,dx), L2diagpMm1(m,n,dx))
	return spdiagm( diags, ( -(m-1), -m, -1, 0, 1, m, (m-1) ) )
end

@everywhere function laplace_diags(m,n)
	diag_seg	= [ 1; 4*ones(m-2) ; 1]
	diag		= [ ones(m); repmat(diag_seg, n-2); ones(m) ]

	ndiag1_seg	= [ 0; -ones( (m-2) ); 0]
	ndiagl1		= [ zeros(m-1); repmat(ndiag1_seg, n-2); zeros(m) ]
	ndiagr1		= [ zeros(m  ); repmat(ndiag1_seg, n-2); zeros(m-1) ]

	ndiagl2		= [ repmat(ndiag1_seg, n-2); zeros(m) ]
	ndiagr2		= [ zeros(m); repmat(ndiag1_seg, n-2) ]

	return ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2
end

@everywhere function generate_laplace(m,n,dx)
	println("generate Laplace Matrix")
	return spdiagm( laplace_diags(m,n), (-m, -1, 0, 1, m) ) / (dx*dx)
end

#thr! teste, ob das wirklich stimmt!
@everywhere function generate_block_laplace(m, n, T, dt, dx) 
	ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2 = laplace_diags(m,n)
	block_diag			= [ diag/2; repmat(diag, T-3); diag/2 ] 
	block_ndiagl1		= [ ndiagl1/2; repmat([0;ndiagl1], T-3); 0; ndiagl1/2 ]
	block_ndiagr1		= [ ndiagr1/2; repmat([0;ndiagr1], T-3); 0; ndiagr1/2 ]

	block_ndiagl2		= [ ndiagl2/2; repmat([zeros(m) ; ndiagl2], T-3) ; zeros(m) ; ndiagl2/2 ]
	block_ndiagr2		= [ ndiagr2/2; repmat([zeros(m) ; ndiagr2], T-3) ; zeros(m) ; ndiagr2/2 ]

	#thr *dt^2 ?
	return spdiagm( (block_ndiagl2, block_ndiagl1, block_diag, block_ndiagr1, block_ndiagr2), (-m, -1, 0, 1, m) ) * dt^2 / (dx*dx)
end

@everywhere function generate_ellip_beta(m, n, T, dt, dx, alpha, beta)
	println("generate elliptic operator")
	LT		= generate_block_laplace(m,n,T,dt, dx)

	# thr! hier wird das globale m verwendet
	R_diag	= [ones(m*n); 2*ones(m*n*(T-3)); ones(m*n)]
	R_ndiag	= -ones(m*n*(T-2))

	R		= spdiagm( (R_ndiag, R_diag, R_ndiag), (-(m*n), 0, m*n) )

    ellOp = LT + R

    GradNormOp = (LT + R )/dt
    CostNormOp = (alpha * LT + beta * R)/dt

	return ellOp, GradNormOp, CostNormOp, LT/dt, R/dt
end

@everywhere function ellop_inner_hilfsmatrizen(m,n,k)
	_m, _n	= m-2, n-2
	dh0 	= [[1] ; 2*ones(k-3); [1]]
	dhpm1 	= -ones(k-2)
	ht		= spdiagm( (dhpm1, dh0, dhpm1), (-1,0,1) )
	hl		= spdiagm( [[0.5] ; ones(k-3) ; [0.5] ] , 0)
	return ht, hl
end

@everywhere function _generate_ellop_beta(m,n,T, dt, dx, alpha, beta)
	ht, hl		= ellop_inner_hilfsmatrizen(m,n,T)
	L			= generate_laplace(m,n,dx)
	L2			= kron(hl,L) *dt^2

	# das ist T in der niederschrift
	R			= kron(ht, spdiagm(ones(m*n), 0)) 
	ellOp		= L2 + R
	GradNormOp	= (L2 + R) /dt
	CostNormOp	= (alpha * L2 + beta * R)/dt

	return ellOp, GradNormOp, CostNormOp, L2/dt, R/dt
end

# testn	= 100
# testm	= 100
# testT	= 100
# 
# testdx	= 0.5
# testdt	= 0.4
# 
# testalph	= 0.2
# testbet	= 0.3
# 
# @time el , gr , cst , lt , r =generate_ellip_beta(testm, testn, testT, testdt, testdx, testalph, testbet)
# @time el_, gr_, cst_, lt_, r_ =_generate_ellop_beta(testm, testn, testT, testdt, testdx, testalph, testbet)
# 
# println(el == el )
# println(gr == gr )
# println(cst == cst )
# println(lt == lt )
# println(r  ==r)
