# Cx, Cy zentrale finite Differenzen
# Dx, Dy ?

function need_redefine()
	return !(isdefined(:L) && size(L) == (m*n, m*n) && L[1,1]==dx*dx)
	#return !(isdefined(:WaveOp) && size(WaveOp) == (m*n*T)^2 && L[1,1]==dx*dx)
end

function  generateMatrices3(n, h)
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
	T	= copy(Cx)

	Cx[n+1:end-n, 1:end-2*n] = -tCx
	Cx[n+1:end-n, 2*n+1:end] = Cx[n+1:end-n, 2*n+1:end] + tCx

	T[n+1:end-n, n+1:end-n] = -2*tCx

	Dx	= T + abs(Cx)
	Dy	= T + abs(Cy)

	Cy	= Cy / 2 / h
	Cx	= Cx / 2 / h
	Dy	= Dy / 2 / h
	Dx	= Dx / 2 / h


	#println("Cx",  Cx)
	#println("Cy ", Cy)
	#println("Dx ", Dx)
	#println("Dy ", Dy)

	return  Cx, Cy, Dx, Dy
end

function generateB(n, h)
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

	## bis hier alle gecheckt

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

function _generate_laplace(m,n,dx) 
	ind_i_diag	= 1:m*n
	ind_j_diag 	= 1:m*n

	diag_seg	= [ 1; repmat([4], n-2) / (dx*dx); 1]
	diag		= [ ones(n); repmat(diag_seg, m-2); ones(n) ]

	ind_j_ndiag = zeros(1:(n-2)*(m-2)*4) # das hier ist nur zur Deklaration der Variablen
	ind_i_ndiag = zeros(1:(n-2)*(m-2)*4) # als Index (Integerarray) in der gewuenschten Groesse
	lap_ndiag_ind!(ind_i_ndiag, ind_j_ndiag, m, n)

	ndiag	 	= -ones( (n-2)*(m-2)*4 )/ (dx*dx)

	ind_i		= [ind_i_diag; ind_i_ndiag]
	ind_j		= [ind_j_diag; ind_j_ndiag]
	val			= [diag; ndiag]

	return sparse(ind_i, ind_j, val, m*n, m*n) 
end

# die Eintraege fuer den Rand stimmen nicht, die duerfen nicht durch dx^2 geteilt werden.
function __generate_laplace(m,n,dx)
	ind_i_diag	= 1:m*n
	ind_j_diag 	= 1:m*n

	diag_seg	= [ 1; repmat([4], n-2) ; 1]
	diag		= [ ones(n); repmat(diag_seg, m-2); ones(n) ]

	ind_j_ndiag = zeros(1:(n-2)*(m-2)*4) # das hier ist nur zur Deklaration der Variablen
	ind_i_ndiag = zeros(1:(n-2)*(m-2)*4) # als Index (Integerarray) in der gewuenschten Groesse
	lap_ndiag_ind!(ind_i_ndiag, ind_j_ndiag, m, n)

	ndiag	 	= -ones( (n-2)*(m-2)*4 )

	ind_i		= [ind_i_diag; ind_i_ndiag]
	ind_j		= [ind_j_diag; ind_j_ndiag]
	val			= [diag; ndiag]

	return sparse(ind_i, ind_j, val, m*n, m*n) / (dx*dx)
end

function lap_ndiag_ind!(ind_i_ndiag, ind_j_ndiag, m, n)
	for i = 2:m-1
		for j = 2:n-1
			off	= ((i-2)*(n-2)+j-2)*4
			diag= (i-1)*n+j
			#println( ((i-2)*(n-2)+j-2)*4,' ', (i-1)*n+j )

			#ind_j_ndiag[ off+1:off+4  ] = [(i-1)*n+j-n, (i-1)*n+j-1, (i-1)*n+j+1,(i-1)*n+j+n]
			#ind_i_ndiag[ off+1:off+4  ] = [(i-1)*n+j,   (i-1)*n+j,   (i-1)*n+j,  (i-1)*n+j  ]
			
			ind_j_ndiag[ off+1 ] = diag -n
			ind_j_ndiag[ off+2 ] = diag -1
			ind_j_ndiag[ off+3 ] = diag +1
			ind_j_ndiag[ off+4 ] = diag +n

			ind_i_ndiag[ off+1 ] = diag
			ind_i_ndiag[ off+2 ] = diag
			ind_i_ndiag[ off+3 ] = diag
			ind_i_ndiag[ off+4 ] = diag
		end
	end
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
	return spdiagm( laplace_diags(m,n), (-n, -1, 0, 1, n) )/ (dx*dx)
end

function generate_block_laplace(m,n,T,dx)
	ndiagl2, ndiagl1, diag, ndiagr1, ndiagr2 = laplace_diags(m,n)
	block_diag			= [ diag/2; repmat(diag, T-3); diag/2 ] 
	block_ndiagl1		= [ ndiagl1/2; 0; repmat(ndiagl1, T-3) ; 0; ndiagl1/2 ]
	block_ndiagr1		= [ ndiagr1/2; 0; repmat(ndiagr1, T-3) ; 0; ndiagr1/2 ]
	block_ndiagl2		= [ ndiagl2/2; zeros(4) ; repmat(ndiagl2, T-3) ; zeros(4) ; ndiagl2/2 ]
	block_ndiagr2		= [ ndiagr2/2; zeros(4) ; repmat(ndiagr2, T-3) ; zeros(4) ; ndiagr2/2 ]
	return spdiagm( (block_ndiagl2, block_ndiagl1, block_diag, block_ndiagr1, block_ndiagr2), (-n, -1, 0, 1, n) ) * dt^2 / (dx*dx)
end

#need_redefine() && println("matrizen neu machen")
#need_redefine() && const L	= generate_laplace(m, n, dx) 
#need_redefine() && const LU	= factorize(L)
#need_redefine() && const B	= generateB(m, dx)
#need_redefine() && const Cx, Cy, Dx, Dy	= generateMatrices3(n, dx) #thr

const L					= generate_laplace(m, n, dx) 
#const L_				= __generate_laplace(m, n, dx) 

const LU				= factorize(L)
const B					= generateB(m, dx)
const Cx, Cy, Dx, Dy	= generateMatrices3(n, dx) #thr

function generate_wave_op(n, T, dt, alpha, beta)
	LT		= generate_block_laplace(m,n,T,dx)

	R_diag	= [ones(m*n); 2*ones(m*n*(T-3)); ones(m*n)]
	R_ndiag	= -ones(m*n*(T-2))
	R		= spdiagm( (R_ndiag, R_diag, R_ndiag), (-(m*n), 0, m*n) )

    WaveOp = LT + R
    #%GradNormOp = WaveOp / sTime.dt;
    #%CostNormOp = WaveOp / sTime.dt * sRegParam.alpha;

    GradNormOp = (LT + R )/dt
    CostNormOp = (alpha * LT + beta * R)/dt

	println("factorize")
	#WaveOpLU	= factorize(WaveOp)
	#WaveOpLU	= lufact(WaveOp)
	WaveOpLU	= WaveOp
	println("factorized")

	return WaveOp, WaveOpLU, GradNormOp, CostNormOp
end

function __generate_wave_op(n, T, dt, alpha, beta)
	nDOF		= n^2

	Id	= speye(nDOF)

	# thr langsam. mit repmat geht das wohl besser
    LT	= copy(L)
	R	= 2*Id
	R2	= -Id
	# thr
    for k = 3:T-2
	print(k)
        LT	= blkdiag(LT, L)
        R	= blkdiag(R, 2*Id)
        R2	= blkdiag(R2, -Id)
    end

    LT	= blkdiag(L/2,LT,L/2)*dt^2

    To	= blkdiag(-Id, R2)
    Tu	= blkdiag(R2, -Id)
    R	= blkdiag(Id, R, Id)
       
    R[nDOF+1:end, 1:end-nDOF] = R[nDOF+1:end, 1:end-nDOF] + Tu
    R[1:end-nDOF, nDOF+1:end] = R[1:end-nDOF, nDOF+1:end] + To

    # hae?
	#R[nDOF+1:end-nDOF,:] = R[nDOF+1:end-nDOF,:]

    WaveOp = LT + R
    #%GradNormOp = WaveOp / sTime.dt;
    #%CostNormOp = WaveOp / sTime.dt * sRegParam.alpha;
    GradNormOp = (LT + R )/dt
    CostNormOp = (alpha * LT + beta * R)/dt

	println("factorize")
	#WaveOpLU	= factorize(WaveOp)
	#WaveOpLU	= lufact(WaveOp)
	WaveOpLU	= WaveOp
	println("factorized")

	return WaveOp, WaveOpLU, GradNormOp, CostNormOp
end

#need_redefine() && 
@time const WaveOp, WaveOpLU, GradNormOp, CostNormOp	= generate_wave_op(n, T, dt, alpha, beta)
#@time const WaveOp_, WaveOpLU_, GradNormOp_, CostNormOp_	= __generate_wave_op(n, T, dt, alpha, beta)

#echo(WaveOp ==WaveOp_)
#echo(WaveOpLU ==WaveOpLU_)
#echo(GradNormOp ==GradNormOp_)
#echo(CostNormOp ==CostNormOp_	)

1
