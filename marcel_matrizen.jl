# Cx, Cy zentrale finite Differenzen
# Dx, Dy ?

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

function generate_laplace(m,n,dx)
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

const L					= generate_laplace(m, n, dx) 
const LU				= factorize(L)
const B					= generateB(m, dx)
const Cx, Cy, Dx, Dy	= generateMatrices3(n, dx) #thr

function generate_wave_op(n, schritte, dt, alpha, beta)
	nDOF		= n^2

	#B	= spdiagm(vec(ones(nDOF,1)),0);
	Id	= speye(nDOF)

    L2	= copy(L)
	R	= 2*Id
	R2	= -Id
    for k = 2:schritte-1
        L2	= blkdiag(L2, L)
        R	= blkdiag(R, 2*Id)
        R2	= blkdiag(R2, -Id)
    end

    L2	= blkdiag(L/2,L2,L/2)*dt^2

    To	= blkdiag(-Id, R2)
    Tu	= blkdiag(R2, -Id)
    R	= blkdiag(Id, R, Id)
       
    R[nDOF+1:end, 1:end-nDOF] = R[nDOF+1:end, 1:end-nDOF] + Tu
    R[1:end-nDOF, nDOF+1:end] = R[1:end-nDOF, nDOF+1:end] + To

    # hae?
	#R[nDOF+1:end-nDOF,:] = R[nDOF+1:end-nDOF,:]

    WaveOp = L2 + R
    #%GradNormOp = WaveOp / sTime.dt;
    #%CostNormOp = WaveOp / sTime.dt * sRegParam.alpha;
    GradNormOp = (L2 + R )/dt
    CostNormOp = (alpha * L2 + beta * R)/dt
    
    #setup.type = 'nofill';
    #[WaveL, WaveU] = ilu(WaveOp, setup);
	return L2, R, WaveOp, GradNormOp, CostNormOp
end

#L2, R, WaveOp, GradNormOp, CostNormOp	= generate_wave_op(n, T, dt, alpha, beta)

