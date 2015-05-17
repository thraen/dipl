include("marcel_matrizen.jl")

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

	# where LU is either precalculated LU-Decomposition of A or A itself, both work
	u 			= LU\b 

	return reshape(u, m, n)
end

