# Cx, Cy zentrale finite Differenzen
# Dx, Dy ?
# B      identitaet

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

	#B	= spdiagm(ones(nDOF,1),0,nDOF,nDOF)
	B	= spdiagm(vec(ones(nDOF,1)),0);

	#println("Cx",  Cx)
	#println("Cy ", Cy)
	#println("Dx ", Dx)
	#println("Dy ", Dy)
	#println("B ", B)

	#return  B, transpose(Cx), transpose(Cy), Dx, Dy
	return  B, Cx, Cy, Dx, Dy
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
