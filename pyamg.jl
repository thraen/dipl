@everywhere using PyCall
@everywhere @pyimport pyamg
@everywhere @pyimport scipy.sparse as scipy_sparse

@everywhere function pycsr(A::SparseMatrixCSC)
	pya	= scipy_sparse.csc_matrix( (A.nzval, A.rowval-1, A.colptr-1) )
	return pya[:tocsr]()
end

@everywhere function construct_mgsolv(A)
	#print("load into Python")
	pyA	= pycsr(A)
	#print("construct multig solver")

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
#@everywhere N = 1000
#@everywhere Ltest = spdiagm((-ones(N-1), 2*ones(N), -ones(N-1)), (-1,0,1), N, N) * N^2
#@everywhere pyB = kron(speye(N), Ltest) + kron(Ltest, speye(N))

#@everywhere m1	= construct_mgsolv(pyB)
#@everywhere m2	= construct_mgsolv(pyB)

#@everywhere function solve_lin_test(A,b)
	#return A[:solve](b, tol=1e-3)
#end

@everywhere function solve_lin_test2(b)
	return m2[:solve](b, tol=1e-3)
end

@everywhere function solve_lin_test1(b)
	return m1[:solve](b, tol=1e-3)
end

function do_par(b)
	#	die varianten gehen alle nicht
	#cannot serialize pointer
	#e1 = @spawn solve_lin_test(m1, b)
	#e2 = @spawn solve_lin_test(m2, b)
	#cannot serialize pointer
	#remotecall(2, solve_lin_test, m1, b)
	#remotecall(3, solve_lin_test, m2, b)
	#cannot serialize pointer
	#@spawn m1[:solve](b)
	#@spawn m2[:solve](b)


	@time e1 = remotecall(2, solve_lin_test1, b)
	@time e2 = remotecall(3, solve_lin_test2, b)

	@time res1 = fetch(e1)
	@time res2 = fetch(e2)


	return res1, res2

	#return fetch(e1), fetch(e2)

	#return solve_lin_test1(b), solve_lin_test2(b)
end

function do_ser(b)
	@time res1 = solve_lin_test1(b)
	@time res2 = solve_lin_test2(b)

	return res1, res2
end

function test_mg_ser()
	for i=1:3
		b = rand(size(pyB,1))
		e1, e2 = do_ser(b)
		@show e1 == e1
	end
end

function test_mg_par()
	for i=1:3
		b = rand(size(pyB,1))
		e1, e2 = do_par(b)
		@show e1 == e1
	end
end

#println("\nparallel execution")
#@time test_mg_par()
#println("\nserial execution")
#@time test_mg_ser()


#pyB	= WaveOp

# solve with random RHS

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

