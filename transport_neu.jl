@everywhere const r			= dt/dx # thr marcel hat hier nur /dx

@everywhere function fluss_lxw1( a, u, v )
	return 0.5* ( a*(u+v) + r*a*a*(u-v) )
end

@everywhere function fluss_upw1( a, u, v )
	# wenn u bzw v <0 wird von ausserhalb -u bzw -v uebergeben, daher ist die Unterscheidung hier unnoetig
	return (a>=0) ? (a*u) : (a*v)
	#return a*u
end

@everywhere function theta( um, u, up )
	return (u - um +eps()) / (up - u +eps())
end

@everywhere function sbee( thet )
	return max( 0.0, max( min(1.0, 2*thet), min(thet, 2.0) ) );
end

@everywhere function fluss_lim1( a, um, u, up )
	thet = theta(um, u, up)
	ret = (1- sbee( thet ) )* fluss_upw1(a, u, up ) +     sbee( thet )  * fluss_lxw1(a, u, up )
	return 	  ret
end

@everywhere function fluss_lim_kons( a, umm, um, u, up )
	# donor cell flux as lower order flux
	anteil_low = fluss_upw1(a, um, u)
	# ratio of gradient of the current cell face and gradient of the next face 
	# in opposite upwind direction
	if a>=0 
		thet = (um - umm +eps()) / (u - um +eps())
	else
		thet = (up - u   +eps()) / (u - um +eps())
	end
	# this is the difference between lax-wendroff and donor cell fluxes, 
	# weightened by flux limiter dependent on gradient ratio
	anteil_antidiff	= sbee(thet) * 0.5*abs(a) * (1-r*abs(a)) * (u-um)
	return 	  anteil_low + anteil_antidiff
end

#
# approximiert u_t = w*grad(u)
#
function transport(I0, u, v, schritte)
	m, n	= size(I0)
	I		= zeros(m, n, schritte+1)
	I[:,:,1]= I0
	Ih		= zeros(m, n) #half step buffer for dimension splitting
	println("==============Transport=================$n x $m x $T")

	for t = 1:schritte
		for j = 3:n-2
			for i = 3:m-2
				#anteilx = (u[i,j,t] >= 0) 	? fluss_lim1(  u[i,j,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u[i,j,t], I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) : fluss_lim1( -u[i,j,t], I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u[i,j,t], I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
				if (u[i,j,t] >= 0) 	
					anteilx = fluss_lim1(  u[i,j,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u[i,j,t], I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) 
				else
					anteilx = fluss_lim1( -u[i,j,t], I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u[i,j,t], I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
				end
				Ih[i,j] = I[i,j,t] - r* (anteilx)
			end
		end

		for j = 3:n-2 
			for i = 3:m-2
				anteily = (v[i,j,t] >= 0) 	? fluss_lim1(  v[i,j,t], Ih[i-1,j], Ih[i,j], Ih[i+1,j] ) - fluss_lim1(  v[i,j,t], Ih[i-2,j], Ih[i-1,j], Ih[i,j] ) : fluss_lim1( -v[i,j,t], Ih[i+1,j], Ih[i,j], Ih[i-1,j] ) - fluss_lim1( -v[i,j,t], Ih[i+2,j], Ih[i+1,j], Ih[i,j] ) 
				I[i,j,t+1] = Ih[i,j] - r* (anteily)
			end
		end
	end

	return I
end

function ruecktransport(s, I, u, v, n_samp, n_zsamp, norm_s)
	m, n, T = size(I)
	p		= zeros(m,n,T)
	ph		= zeros(m,n)
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
				uxph	= (u[i,j,t-1] + u[i,j+1,t-1])/2
				uxmh	= (u[i,j,t-1] + u[i,j-1,t-1])/2
				anteilx = fluss_lim_kons( uxph, p[i,j-1,t], p[i,j,t], p[i,j+1,t], p[i,j+2,t]) - fluss_lim_kons( uxmh, p[i,j-2,t], p[i,j-1,t], p[i,j,t], p[i,j+1,t])

				ph[i,j] = p[i,j,t] - r* anteilx
			end
		end

		for j = 3:n-2  # zuerst die spalten. ist etwas schneller
			for i = 3:m-2
				vyph	= (v[i+1,j,t-1] + v[i,j,t-1])/2
				vymh	= (v[i-1,j,t-1] + v[i,j,t-1])/2
				anteily = fluss_lim_kons( vyph, ph[i-1,j], ph[i,j], ph[i+1,j], ph[i+2,j]) - fluss_lim_kons( vymh, ph[i-2,j], ph[i-1,j], ph[i,j], ph[i+1,j])

				p[i,j,t-1] = ph[i,j] - r* anteily
			end
		end
	end

	return p
end
