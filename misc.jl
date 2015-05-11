#Base.Multimedia.writemime(stream,::MIME"text/plain",x::Float64)=@printf("%1.2f",x)

function echo(args...)
	for x in args print( x, '\t') end 
	print('\n')
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

# teste H1_norm mit l2-Norm und L2-Norm
function __H1_norm(u, v)
	return dx*dx* dt* ( sum( central_diff_x( u ).^2 + sum( central_diff_y( u ).^2 ) ) + sum( central_diff_x( v ).^2 + sum( central_diff_y( v ).^2 ) ) )
end

function ___H1_norm(u, v, L)
	ret = 0
	for t=1:T-1
		u_ = reshape(u[:,:,t] , n*m)
		v_ = reshape(v[:,:,t] , n*m)
		ret_  = (u_'*L*u_) 
		ret_ += (v_'*L*v_) 
		if t==1 || t==T
			ret_ /= 2
		end
		ret += ret_
	end
	return dx*dx* dt* ret[1]
end

function l2norm(s)
	#l2_s	= sum([ sum(s[:,:,k].^2) for k=1:n_samples ])
	return dx*dx* dt* sum(s.^2)
end

# berechnet sum( L_X(s[:,:,k ) , k in I )
# s[:,:,k], k in 0:T Bildsequenz
# X quadratische Matrix
# I Indexmenge
function Xnorm(s, X, Idx) 
	ret	= 0
	for k in Idx
		s_	 = reshape(s[:,:,k] , n*m)
		ret	+= s_'*X*s_
	end
	return ret[1]
end

function Xnorm(s, X) 
	ret	= 0
	for k = 1:n_samples
		s_	 = reshape(s[:,:,k] , n*m)
		ret	+= s_'*X*s_
	end
	return ret[1]
end

B		= generateB(m, dx)
Beye	= speye(m*n)

function sample_err_L2(I, s, norm_s)
	err		= zeros( size(s) )
	L2err	= 0
	for (k,t) in sample_times
		err[:,:,k]	= I[:,:,t] - s[:,:,k]
	end
	L2err 	= dx*dx* L2norm(err)/norm_s
	return L2err, err
end

function sample_err_l2(I, s, norm_s)
	err		= zeros( size(s) )
	l2err	= 0
	for (k,t) in sample_times
		err[:,:,k]	= I[:,:,t] - s[:,:,k]
	end
	# dasselbe wie, aber schneller als l2err=Xnorm(err, speye, ..)
	l2err = dx*dx* sum(err.*err)/norm_s 
	return l2err, err
end

function L2norm(s)
	return Xnorm(s,B)
	#return Xnorm(s,Beye)
end

function sample_err(I, s, norm_s)
	return sample_err_L2(I, s, norm_s)
	#return sample_err_l2(I, s, norm_s)
end

function H1_norm(u,v)
	___H1_norm(u,v,L)
end
#H1_norm = __H1_norm



