include("matrizen.jl")

function save_jld(steps, dx, dt, alpha, beta, s, I, p, u, v, grd_u_J, grd_v_J)
	try
		info("\n zwischenspeichern\n $steps")
		save("$(rootdir)zwischenergebnis_$steps.jld", 
				"dx", dx,
				"dt", dt,
				"alpha", alpha,
				"beta", beta,
				"s", s,
				"I", I, 
				"p", p,
				"u", u,
				"v", v, 
				"grd_u_J", grd_u_J, 
				"grd_v_J", grd_v_J)	
	catch e
		warn("ZWISCHENERGEBNIS KONNTE NICHT GESPEICHERT WERDEN!", e)
	end
end

macro init_grad(m,n,T)
	if grad_parallel 
		echo("***** init grad parallel *****")
		return :( SharedArray(Float64, ($m, $n, T-1), init= S -> S[localindexes(S)] = 0.0) )
	else
		echo("***** init grad serial *****")
		return :( zeros( $m, $n, T-1 ) )
	end
end

macro do_par_when_defined(forloop)
	#print(forloop, "\n")
	return grad_parallel ? :( @sync @parallel $(forloop) ) : forloop
end

function central_diff_y( a )
	m, n = size(a)
	d = zeros(m-2, n-2)
	for j= 2:n-2
		for i= 2:m-2
			d[i-1,j-1] = (a[i+1,j] - a[i-1,j])/2
		end
	end
	return d
end

function central_diff_x( a )
	m, n = size(a)
	d = zeros(m-2, n-2)
	for j= 2:n-2
		for i= 2:m-2
			d[i-1,j-1] = (a[i,j+1] - a[i,j-1])/2
		end
	end
	return d
end

function l2norm(s)
	#l2_s	= sum([ sum(s[:,:,k].^2) for k=1:n_samples ])
	# thr! hier darf kein dt rein?
	return dx*dx* sum(s.^2)
end

function Xnorm(s, X) 
	ret	= 0
	T	= size(s,3)
	for k = 1:T
		s_	 = reshape(s[:,:,k] , n*m)
		ret	+= s_'*X*s_
	end
	return ret[1]
end

function sample_err_L2(I, s, norm_s)
	err		= zeros( size(s) )
	L2err	= 0
	for (k,t) in samples_to_frames
		err[:,:,k]	= I[:,:,t] - s[:,:,k]
	end
	L2err 	= L2norm(err)/norm_s
	return L2err, err
end

function sample_err_l2(I, s, norm_s)
	err		= zeros( size(s) )
	l2err	= 0
	for (k,t) in samples_to_frames
		err[:,:,k]	= I[:,:,t] - s[:,:,k]
	end
	# dasselbe wie, aber schneller als l2err=Xnorm(err, speye)
	l2err = dx*dx* sum(err.*err)/norm_s 
	return l2err, err
end

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
function poissolv(f, gu, gd, gl, gr)
	m = size(f, 1) +2
	n = size(f, 2) +2

	# boundary condition
	g			= [gu ; [gl f gr] ; gd]
	b			= reshape(g, m*n)

	# where LU is either precalculated LU-Decomposition of A or A itself, both work
	u 			= LU\b 

	return reshape(u, m, n)
end

# hier wird b zusammengesetzt aus f und g uebergeben
function poissolv_(b,m, n)
	# where LU is either precalculated LU-Decomposition of A or A itself, both work
	u 			= LU\b 
	return reshape(u, m, n)
end

# thr
# L*f ~= ∆f !
# am Rand stimmt's nicht! in unseren Beispielen egal, da dort und in der Naehe eh alles 0 ist
# das passiert in grad_J_beta, und grad_J_beta_parallel und h1norm
function laplace(f)
	# Lf		= L*f # wrong on boundary 
	# corrected = extrapolate inner values to boundary # geht allerdings von stetigem ∆f
	# return 	corrected
end

function checkdiv(u, v)
	mu, nu, T = size(u)
	mv, nv, T = size(v)
	#divx	= 0
	#divy	= 0

	mdu, ndu	= size( central_diff_x( u[:,:,1] ) )
	mdv, ndv	= size( central_diff_y( v[:,:,1] ) )
	divx = zeros( mdu, ndu, T)
	divy = zeros( mdv, ndv, T)
	for t=1:T
		divx[:,:,t] = central_diff_x( u[:,:,t] )
		divy[:,:,t] = central_diff_y( v[:,:,t] )
	end
	
	return divx , divy
end


