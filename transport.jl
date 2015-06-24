@everywhere const r			= dt/(dx*dx) # thr marcel hat hier nur /dx
#@everywhere const r			= dt/dx # thr marcel hat hier nur /dx

#const infts = [999999.0, 999999.0, 999999.0, 999999.0]

@everywhere const infts = 999999
@everywhere const oness = 1
@everywhere const twos	= 2
@everywhere const zeross= 0

@everywhere function fluss_lxw1( a, u, v )
	return 0.5* ( a*(u+v) + r*a*a*(u-v) )
end

@everywhere function fluss_upw1( a, u, v )
	# wenn u bzw v <0 wird von ausserhalb -u bzw -v uebergeben, daher ist die Unterscheidung hier unnoetig
	#return (a>=0) ? (a*u) : (a*v)
	return a*u
end

@everywhere function theta( um, u, up )
	if ( (up-u) == 0 )
		return infts
	else
		return (u - um) / (up - u)
	end
end

@everywhere function sbee( thet )
	return max( zeross, max( min(oness, 2*thet), min(thet, twos) ) );
end

@everywhere function fluss_lim1( a, um, u, up )
	thet = theta(um, u, up)
	ret = (1- sbee( thet ) )* fluss_upw1(a, u, up ) +     sbee( thet )  * fluss_lxw1(a, u, up )
	return 	  ret
end

#
# approximiert u_t = w*grad(u)
#

function transport(I0, u, v, schritte)
	m, n	= size(I0)
	I		= zeros(m, n, schritte+1)
	I[:,:,1]= I0
	println("==============Transport=================$n x $m x $T")

	for t = 1:schritte
		for j = 3:n-2  # zuerst die spalten. ist etwas schneller
			for i = 3:m-2
				diffx = u[i,j+1,t] - u[i,j-1,t]
				diffy = v[i+1,j,t] - v[i-1,j,t]

				#echo(diffx, diffy)

				anteilx = (u[i,j,t] >= 0) ? fluss_lim1(  u[i,j,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u[i,j,t], I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) : fluss_lim1( -u[i,j,t], I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u[i,j,t], I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
				anteily = (v[i,j,t] >= 0) ? fluss_lim1(  v[i,j,t], I[i-1,j,t], I[i,j,t], I[i+1,j,t] ) - fluss_lim1(  v[i,j,t], I[i-2,j,t], I[i-1,j,t], I[i,j,t] ) : fluss_lim1( -v[i,j,t], I[i+1,j,t], I[i,j,t], I[i-1,j,t] ) - fluss_lim1( -v[i,j,t], I[i+2,j,t], I[i+1,j,t], I[i,j,t] ) 

				#I[i,j,t+1] = I[i,j,t] - r* ( anteilx + anteily ) 
				I[i,j,t+1] = I[i,j,t] - r* ( anteilx + anteily )  + I[i,j,t]*dt/2*(diffx+diffy) ## ausklammern

				#if abs(I[i,j,t]) > 1
					#echo(t,i,j, anteilx, anteily, u[i,j,t], v[i,j,t])
					#echo(I[i,j,t], I[i,j,t+1])
					#bla()
				#end
			end
		end
	end

	return I
end

function ruecktransport(s, I, u, v, n_samp, n_zsamp, norm_s)
	m, n, T = size(I)
	p		= zeros(m,n,T)
	#p		= SharedArray(typeof(I[1,1,1]), size(I), init=false)
	println("==============Ruecktransport============$n x $m x $n_samples $n_zwischensamples, $T")

	# thr!!! rechter oder linker grenzwert zum diracterm?
	sk		= 0

	for t = T:-1:2
		# thr!!! rechter oder linker grenzwert zum diracterm?
		if mod(t-1, n_zsamp) == 0 then
			#echo("am sample        ", t, "->", t-1, n_samp-sk)
			err			= I[:, :, t] - s[:, :, n_samp-sk] 
			p[:,:,t] 	= p[:,:,t] - err/norm_s
			sk 			+= 1
		#else
			#echo("zwischen samples ", t, "->", t-1)
		end

		for j = 3:n-2  # zuerst die spalten. ist etwas schneller
			for i = 3:m-2
				a = u[i,j,t-1]
				b = v[i,j,t-1]
				anteilx = (a >= 0) ? fluss_lim1(  a, p[i,j-1,t], p[i,j,t], p[i,j+1,t] ) - fluss_lim1(  a, p[i,j-2,t], p[i,j-1,t], p[i,j,t] ) : fluss_lim1( -a, p[i,j+1,t], p[i,j,t], p[i,j-1,t] ) - fluss_lim1( -a, p[i,j+2,t], p[i,j+1,t], p[i,j,t] ) 
				anteily = (b >= 0) ? fluss_lim1(  b, p[i-1,j,t], p[i,j,t], p[i+1,j,t] ) - fluss_lim1(  b, p[i-2,j,t], p[i-1,j,t], p[i,j,t] ) : fluss_lim1( -b, p[i+1,j,t], p[i,j,t], p[i-1,j,t] ) - fluss_lim1( -b, p[i+2,j,t], p[i+1,j,t], p[i,j,t] ) 

				p[i,j,t-1] = p[i,j,t] - r* (anteilx + anteily)
			end
		end
	end

	return p
end
