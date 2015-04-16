# we solve the 2-dimensional poisson equation
#
#  laplace(u) = -f
#
# with boundary condition 
#  u(x,y) = g on the border of the rectangular area [1..m]x[1..n]
#
# where f is given as ((m-2) x (n-2))-sample and 
# where g is given as 2 m-2 vectors gl, gr and 2 n vectors gu, gd
#
# 		=======gu======		_
#       |			  |
#		gl     f      gr	m-2
#		|			  |		_
# 		=======gd======		
#
#		|	  n  	  |
#

@everywhere function poissolv(f, gu, gd, gl, gr)
	m = size(f, 1) +2
	n = size(f, 2) +2

	# boundary condition
	g			= [gu ; [gl f gr] ; gd]

	b			= reshape(g, m*n)
	A			= generate_laplace(m,n)

	u 			= A\b

	return reshape(u, m, n)
end

# where LU is either precalculated LU-Decomposition of A or A itself, works both
@everywhere function poissolv(LU, f, gu, gd, gl, gr)
	m = size(f, 1) +2
	n = size(f, 2) +2

	# boundary condition
	g			= [gu ; [gl f gr] ; gd]
	b			= reshape(g, m*n)

	u 			= LU\b 

	return reshape(u, m, n)
end

@everywhere function generate_laplace(m,n)
	ind_i_diag	= 1:m*n
	ind_j_diag 	= 1:m*n

	diag_seg	= [ 1; repmat([4], n-2) ; 1]
	diag		= [ ones(n); repmat(diag_seg, m-2); ones(n) ]

	ind_j_ndiag = zeros(1:(n-2)*(m-2)*4) # das hier ist nur zur Deklaration der Variablen
	ind_i_ndiag = zeros(1:(n-2)*(m-2)*4) # als Index (Integerarray) in der gewuenschten Groesse
	ndiag_ind!(ind_i_ndiag, ind_j_ndiag, m, n)

	ndiag	 	= -ones( (n-2)*(m-2)*4 )

	ind_i		= [ind_i_diag; ind_i_ndiag]
	ind_j		= [ind_j_diag; ind_j_ndiag]
	val			= [diag; ndiag]

	return sparse(ind_i, ind_j, val, m*n, m*n) 
end

@everywhere function ndiag_ind!(ind_i_ndiag, ind_j_ndiag, m, n)
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

