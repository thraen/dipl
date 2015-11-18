
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

function generateB(n, h)
	println("generate B")
	nDOF		= n^2
	B			= spzeros(nDOF, nDOF)

	SW			= spzeros(1, n+1)
	SW[1,1]		= 2
	SW[1,2]		= 1
	SW[1, n+1]	= 1

	W			= spzeros(1, n+2)
	W[1,1]		= 1
	W[1,2]		= 6
	W[1,3]		= 1
	W[1, n+1]	= 2
	W[1, n+2]	= 2

	NW			= spzeros(1, 2*n)
	NW[1, n-1]	= 1
	NW[1, n]	= 4
	NW[1, 2*n-1]= 2
	NW[1, 2*n]	= 1

	S			= spzeros(1, 2*n+1)
	S[1,1]		= 1
	S[1,2]		= 2
	S[1, n+1]	= 6
	S[1, n+2]	= 2
	S[1, 2*n+1]	= 1

	C			= spzeros(1, 2*n+2)
	C[1, 2]		= 2
	C[1, 3]		= 2
	C[1, n+1]	= 2
	C[1, n+2]	= 12
	C[1, n+3]	= 2
	C[1, 2*n+1]	= 2
	C[1, 2*n+2]	= 2

	N			= spzeros(1, 3*n)
	N[1, n]		= 1
	N[1, 2*n-1]	= 2
	N[1, 2*n]	= 6
	N[1, 3*n-1]	= 2
	N[1, 3*n]	= 1

	SE			= spzeros(1, 2*n)
	SE[1,1] = 1
	SE[1,2] = 2
	SE[1, n+1] = 4
	SE[1, n+2] = 1

	E			= spzeros(1, 2*n)
	E[1, end-n-1] = 2
	E[1, end-n] = 2
	E[1, end-2] = 1
	E[1, end-1] = 6
	E[1, end] = 1

	NE			= spzeros(1, 2*n)
	NE[1, end-n] = 1
	NE[1, end-1] = 1
	NE[1, end] = 2

	Links		= spzeros(n, 2*n)

	Links[1, 1:n+1] = SW
	for k = 1:n-2
		Links[k+1, k:1+n+k] = W
	end
	Links[n, 1:2*n] = NW

	Mitte		= spzeros(n, 3*n)

	Mitte[1, 1:2*n+1] = S
	for k = 1:n-2
		Mitte[k+1, k:1+2*n+k] = C
	end
	Mitte[n, 1:3*n] = N

	Rechts = spzeros(n, 3*n)

	Rechts[end, end-2*n+1:end] = NE
	for k = 1:n-2
		Rechts[end-k, end-2*n-k+2:end-k+1] = E 
	end
	Rechts[end-n+1, end-2*n+1:end] = SE 

	B[1:n, 1:2*n] = Links
	for k = 1:n-2
	   B[k*n+1:(k+1)*n, (k-1)*n+1:(k-1)*n+3*n] =  Mitte
	end
	B[end-n+1:end, end-3*n+1:end] = Rechts

	return B * h^2/24
end

function laplace_diags(m,n)
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

function generate_block_laplace(m, n, T, dt, dx)
	ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2 = laplace_diags(m,n)
	block_diag			= [ diag/2; repmat(diag, T-3); diag/2 ] 
	block_ndiagl1		= [ ndiagl1/2; repmat([0;ndiagl1], T-3); 0; ndiagl1/2 ]
	block_ndiagr1		= [ ndiagr1/2; repmat([0;ndiagr1], T-3); 0; ndiagr1/2 ]
	block_ndiagl2		= [ ndiagl2/2; repmat([zeros(n) ; ndiagl2], T-3) ; zeros(n) ; ndiagl2/2 ]
	block_ndiagr2		= [ ndiagr2/2; repmat([zeros(n) ; ndiagr2], T-3) ; zeros(n) ; ndiagr2/2 ]
	return spdiagm( (block_ndiagl2, block_ndiagl1, block_diag, block_ndiagr1, block_ndiagr2), (-n, -1, 0, 1, n) ) * dt^2 / (dx*dx)
end

function generate_ellip_beta(n, T, dt, dx, alpha, beta)
	println("elliptischer Operatormatrix")
	LT		= generate_block_laplace(m,n,T,dt, dx)

	R_diag	= [ones(m*n); 2*ones(m*n*(T-3)); ones(m*n)]
	R_ndiag	= -ones(m*n*(T-2))

	R		= spdiagm( (R_ndiag, R_diag, R_ndiag), (-(m*n), 0, m*n) )

    ellOp = LT + R
	# thr, was hats damit auf sich
    #GradNormOp = ellOp / sTime.dt;
    #CostNormOp = ellOp / sTime.dt * sRegParam.alpha;

    GradNormOp = (LT + R )/dt
    CostNormOp = (alpha * LT + beta * R)/dt

	println("factorize")
	#ellOPLU	= factorize(ellOp)
	#ellOPLU	= lufact(ellOp)
	ellOPLU	= ellOp
	println("factorized")

	return ellOp, ellOPLU, GradNormOp, CostNormOp
end
#(isdefined(:L)		&& (m*n==size(L,1))) 			|| (const L 	= generate_laplace(m, n, dx))
#(isdefined(:LU)		&& (m*n==size(L,1)))			|| (const LU	= factorize(L))
#(isdefined(:B)		&& (m*n==size(B,1))) 			|| (const B		= generateB(m, dx))
#(isdefined(:Cx)		&& (m*n==size(Cx,1))) 			|| (const Cx, Cy = generate_differentiation_central(n, dx) )
#(isdefined(:ellOp)	&& (m*n*(T-1)==size(ellOp,1)))	|| (const ellOp, ellOPLU, GradNormOp, CostNormOp	= generate_ellip_beta(n, T, dt, dx, alpha, beta))


const L 	= generate_laplace(m, n, dx)
const LU	= factorize(L)
const B		= generateB(m, dx)
const Cx, Cy = generate_differentiation_central(n, dx) 
const ellOp, ellOPLU, GradNormOp, CostNormOp	= generate_ellip_beta(n, T, dt, dx, alpha, beta)

include("pyamg.jl")
const ml			= construct_mgsolv(ellOp)

include("matrizen_zellgrenzen.jl")
const Lx	= generateLu(m,n,dx)
const Ly	= generateLv(m,n,dx)
const Cx_zg	= generate_differentiation_interfx(m,n,dx)
const Cy_zg	= generate_differentiation_interfy(m,n,dx)

# Interpolationsmatrizen
const P_zgx = abs(Cx_zg)*dx/2
const P_zgy = abs(Cy_zg)*dx/2

