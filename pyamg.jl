@everywhere using PyCall
@everywhere @pyimport pyamg
@everywhere @pyimport scipy.sparse as scipy_sparse

@everywhere function pycsr(A::SparseMatrixCSC)
	pya	= scipy_sparse.csc_matrix( (A.nzval, A.rowval-1, A.colptr-1) )
	return pya[:tocsr]()
end

# ist nur fuer symmetrische Matrizen richtig? ja
function _pycsr(A::SparseMatrixCSC)
	return	scipy_sparse.csr_matrix( (A.nzval, A.rowval-1, A.colptr-1) )
end

@everywhere function construct_mgsolv(A)
	print("load into Python")
	@time pyA	= pycsr(A)
	@show size(A)
	print("construct multig solver")

	# geht nicht fuer die Stokesmatrix
	@time ml	= pyamg.ruge_stuben_solver(pyA)

	#@time ml = pyamg.ruge_stuben_solver(pyA, max_levels=10)
	#@time ml = pyamg.ruge_stuben_solver(pyA, max_coarse=10)
	#@time ml = pyamg.ruge_stuben_solver(pyA)

	#geht fuer Stokesmatrix, ist aber langsam
	#@time ml = pyamg.smoothed_aggregation_solver(pyA)

	return ml
end

# generate 2D laplacian
#N = 1000
#L1 = spdiagm((-ones(N-1), 2*ones(N), -ones(N-1)), (-1,0,1), N, N) * N^2
#pyB = kron(speye(N), L1) + kron(L1, speye(N))

#pyB	= WaveOp

# solve with random RHS
#b = rand(size(B,1))

#print("pyamg")
#@time x_py		= ml[:solve](b, tol=1e-4)

#print("gmres")
#@time x_gm,h	= gmres(B, b, restart=5, tol=1e-8, maxiter=100)
#@time x_gm,h	= gmres(B, b, restart=5, tol=1e-8, maxiter=1)
#@time x_gm,h	= gmres(B, b, tol=1e-10)
#@time x_gm,h	= gmres(B, b, restart=5)

#print("exact")
#@time x			= B\ b

# check result
#println("|x - x_py| = ", norm(x_py - x	, Inf))
#println("|x - x_gm| = ", norm(x_gm - x	, Inf))

