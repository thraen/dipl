
function generate_differentiation_central(n, dx)
	println("generate central differences Matrices for cell centers")
	nDOF= n^2
	cy	= sparse(diagm(vec([-ones(n-2,1);0]),-1)+diagm(vec([0;ones(n-2,1)]),1))
	cx	= sparse(diagm(vec([0; ones(n-2,1);0])))

	tCx	= copy(cx)
	
	for k=3:n-1
		tCx = blkdiag(tCx, cx)
	end

	z	= spzeros(n,n)
	Cy	= copy(z)
	for k=2:n-1
		Cy = blkdiag(Cy, cy)
	end

	Cy	= blkdiag(Cy, z)

	Cx	= spzeros(nDOF,nDOF)

	Cx[n+1:end-n, 1:end-2*n] = -tCx
	Cx[n+1:end-n, 2*n+1:end] = Cx[n+1:end-n, 2*n+1:end] + tCx

	Cy	= Cy / 2 / dx
	Cx	= Cx / 2 / dx

	return  Cx, Cy
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
	diag_seg	= [ 1; repmat([4], n-2) ; 1]
	diag		= [ ones(n); repmat(diag_seg, m-2); ones(n) ]

	ndiag1_seg	= [ 0; -ones( (n-2) ); 0]
	ndiagl1		= [ zeros(n-1); repmat(ndiag1_seg, m-2); zeros(n) ]
	ndiagr1		= [ zeros(n  ); repmat(ndiag1_seg, m-2); zeros(n-1) ]

	ndiagl2		= [ repmat(ndiag1_seg, m-2); zeros(n) ]
	ndiagr2		= [ zeros(n); repmat(ndiag1_seg, m-2) ]

	return ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2
end

function generate_laplace(m,n,dx)
	println("generate Laplace Matrix")
	return spdiagm( laplace_diags(m,n), (-n, -1, 0, 1, n) )/ (dx*dx)
end

@everywhere function generate_block_laplace(m, n, T, dt, dx)
	ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2 = laplace_diags(m,n)
	block_diag			= [ diag/2; repmat(diag, T-3); diag/2 ] 
	block_ndiagl1		= [ ndiagl1/2; repmat([0;ndiagl1], T-3); 0; ndiagl1/2 ]
	block_ndiagr1		= [ ndiagr1/2; repmat([0;ndiagr1], T-3); 0; ndiagr1/2 ]
	block_ndiagl2		= [ ndiagl2/2; repmat([zeros(n) ; ndiagl2], T-3) ; zeros(n) ; ndiagl2/2 ]
	block_ndiagr2		= [ ndiagr2/2; repmat([zeros(n) ; ndiagr2], T-3) ; zeros(n) ; ndiagr2/2 ]
	return spdiagm( (block_ndiagl2, block_ndiagl1, block_diag, block_ndiagr1, block_ndiagr2), (-n, -1, 0, 1, n) ) * dt^2 / (dx*dx)
end

@everywhere function generate_ellip_beta(n, T, dt, dx, alpha, beta)
	println("generate elliptic operator")
	LT		= generate_block_laplace(m,n,T,dt, dx)

	R_diag	= [ones(m*n); 2*ones(m*n*(T-3)); ones(m*n)]
	R_ndiag	= -ones(m*n*(T-2))

	R		= spdiagm( (R_ndiag, R_diag, R_ndiag), (-(m*n), 0, m*n) )

    ellOp = LT + R

    GradNormOp = (LT + R )/dt
    CostNormOp = (alpha * LT + beta * R)/dt

	return ellOp, GradNormOp, CostNormOp
end
