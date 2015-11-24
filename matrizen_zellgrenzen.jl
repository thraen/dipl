function generate_differentiation_interfx(m,n,dx)
	println("differentiationmatrix x")
	#Knotenpunkte in x und y-Richtung
	nx = n-1
	ny = m

	# Hauptdiagonale
	d	= -ones(ny*nx)

	# ny. Nebendiagonale
	dny = ones(ny*nx, 1)

	Cx = spdiagm((d, dny), (0, ny), nx*ny, m*n)
	return Cx /dx
end

function generate_differentiation_interfy(m,n,dx)
	println("differentiationmatrix y")
	#Knotenpunkte in x und y-Richtung
	nx = n
	ny = m-1

	# Blockmatrix
	d = ones(ny, 1)
	D = spdiagm((-d, d), (0, 1), ny, m)

	b = speye(nx)
	Cy = kron(b, D) 
	return Cy /dx
end

function generateLu(m,n,dx)
	println("Laplace u")
	#Knotenpunkte in x und y-Richtung
	nx = n-1
	ny = m

	# Hauptdiagonale
	d = [-5; -4*ones(ny-2); -5]
	d = repmat(d, nx)

	# 1. Nebendiagonale. 
	# bei matlab werden die negativen Nebendiagonalen unten abgeschnitten.
	# die positiven Nebendiagonalen werden oben abgeschnitten

	# thr! ist die 0 da drin Absicht? oder sollte die nur nur einmal am anfang eingefuegt werden
	d1 = [0; ones(ny-1)]
	d1 = repmat(d1, nx)[2:end]

	# ny. Nebendiagonale
	dny = ones(ny, 1)
	dny = repmat(dny, nx-1)

	# -1. Nebendiagonale
	dm1 = [ones(ny-1); 0]
	dm1 = repmat(dm1, nx)[1:end-1]

	# ny. Nebendiagonale
	dmny = ones(ny)
	dmny = repmat(dmny, nx-1)

	# Laplace Matrix
	Lu = spdiagm((dmny, dm1, d, d1, dny), (-ny, -1, 0, 1, ny), nx*ny, nx*ny)
	return Lu / dx^2
end

function generateLv(m,n,dx)
	println("Laplace v")
	nx = n
	ny = m-1

	# Hauptdiagonale
	d = -4*ones(ny)
	d = [-5*ones(ny); repmat(d, nx-2); -5*ones(ny)]

	# 1. Nebendiagonale
	d1 = [0; ones(ny-1)]
	d1 = repmat(d1, nx)[2:end]

	# ny. Nebendiagonale
	dN_y = ones(ny)
	dN_y = repmat(dN_y, nx-1)

	# -1. Nebendiagonale
	dm1 = [ones(ny-1); 0]
	dm1 = repmat(dm1, nx)[1:end-1]

	# ny. Nebendiagonale
	dmN_y = ones(ny)
	dmN_y = repmat(dmN_y, nx-1)

	# Laplace Matrix
	Lv = spdiagm((dmN_y, dm1, d, d1, dN_y), (-ny, -1, 0, 1, ny), nx*ny, nx*ny)
	return Lv / dx^2
end

function generateS(m, n, dx, Cx_zg, Cy_zg, Lx, Ly)
	ndofu	= prod(m*(n-1))
	ndofv	= prod((m-1)*n)

	S 		= [	Lx						spzeros(ndofu, ndofv)	Cx_zg;
				spzeros(ndofv, ndofu)	Ly						Cy_zg;
				Cx_zg'					Cy_zg'					spzeros(m*n, m*n) ]

	S[end, :]	= 0
	S[end, end] = 1

	return S
end
