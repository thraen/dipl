#@everywhere const r			= dt/(dx*dx) # thr marcel hat hier nur /dx
@everywhere const r			= dt/dx # thr marcel hat hier nur /dx

#const infts = [999999.0, 999999.0, 999999.0, 999999.0]

@everywhere const infts = 999999.0
@everywhere const oness = 1.0
@everywhere const twos	= 2.0
@everywhere const zeross= 0.0

@everywhere function fluss_lxw1( a, u, v )
	return 0.5* ( a*(u+v) + r*a*a*(u-v) )
end

@everywhere function fluss_upw1( a, u, v )
	# wenn u bzw v <0 wird von ausserhalb -u bzw -v uebergeben, daher ist die Unterscheidung hier unnoetig
	return (a>=0) ? (a*u) : (a*v)
	#return a*u
end

#@everywhere function theta( um, u, up )
	#if ( (up-u) == 0 )
		#return infts
	#else
		#return (u - um) / (up - u)
	#end
#end

@everywhere function theta( um, u, up )
	return (u - um +eps()) / (up - u +eps())
end

@everywhere function sbee( thet )
	return max( zeross, max( min(oness, 2*thet), min(thet, twos) ) );
end

@everywhere function fluss_lim1( a, um, u, up )
	thet = theta(um, u, up)
	ret = (1- sbee( thet ) )* fluss_upw1(a, u, up ) +     sbee( thet )  * fluss_lxw1(a, u, up )
	return 	  ret
end

@everywhere function limited_hot(a, umm, um, u, up)
	if a>=0 
		thet = (um - umm +eps()) / (u - um +eps())
	else
		thet = (up - u   +eps()) / (u - um +eps())
	end
	anteil_hig	= 0.5*abs(a) * (1-r*abs(a)) * sbee(thet) * (u-um)

	return 	  anteil_hig
end

@everywhere function fluss_lim_kons( a, umm, um, u, up )
	anteil_low	= fluss_upw1(a, um, u)
	anteil_hig	= limited_hot(a, umm, um, u, up)
	return 	  anteil_low + anteil_hig
end

function transport(I0, u, v, schritte)
	m, n	= size(I0)
	I		= zeros(m, n, schritte+1)
	I[:,:,1]= I0
	Ih		= zeros(m, n) #half step buffer for dimension splitting
	println("==============Transport=================$n x $m x $T")
	for t = 1:schritte
		for j = 3:n-2
			for i = 3:m-2
				#uxph	= (u[i,j,t] + u[i,j+1,t])/2
				#uxmh	= (u[i,j,t] + u[i,j-1,t])/2

				uxph	= u[i,j,t] 
				uxmh	= u[i,j-1,t]
				
				uxph_m	= min( uxph, 0)
				uxmh_p	= max( uxmh, 0)

				Wmh		= I[i,j,t]   - I[i,j-1,t]
				Wph		= I[i,j+1,t] - I[i,j,t]

				anteil_higmh	= limited_hot( uxmh, I[i,j-2,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] )
				anteil_higph	= limited_hot( uxph, I[i,j-1,t], I[i,j,t], I[i,j+1,t], I[i,j+2,t] )

				anteil_low	= uxph_m*Wph + uxmh_p*Wmh

				Ih[i,j] = I[i,j,t] - r* anteil_low - r*(anteil_higph - anteil_higmh)
			end
		end
		for j = 3:n-2 
			for i = 3:m-2
				#vymh	= (v[i-1,j,t] + v[i,j,t])/2
				#vyph	= (v[i+1,j,t] + v[i,j,t])/2

				vyph	= v[i,j,t]
				vymh	= v[i-1,j,t]

				vymh_p	= max( vymh, 0)
				vyph_m	= min( vyph, 0)

				Wmh		= I[i,j,t]   - I[i-1,j,t]
				Wph		= I[i+1,j,t] - I[i,j,t]

				anteil_low		= vyph_m*Wph + vymh_p*Wmh

				anteil_higmh	= limited_hot( vymh, I[i-2,j,t], I[i-1,j,t], I[i,j,t], I[i+1,j,t] )
				anteil_higph	= limited_hot( vyph, I[i-1,j,t], I[i,j,t], I[i+1,j,t], I[i+2,j,t] )

				I[i,j,t+1] = Ih[i,j] - r* anteil_low - r*(anteil_higph - anteil_higmh)
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
				#uxph	= (u[i,j,t-1] + u[i,j+1,t-1])/2
				#uxmh	= (u[i,j,t-1] + u[i,j-1,t-1])/2

				uxph	= u[i,j,t-1]
				uxmh	= u[i,j-1,t-1]

				anteilx = fluss_lim_kons( uxph, p[i,j-1,t], p[i,j,t], p[i,j+1,t], p[i,j+2,t]) - fluss_lim_kons( uxmh, p[i,j-2,t], p[i,j-1,t], p[i,j,t], p[i,j+1,t])

				ph[i,j] = p[i,j,t] - r* anteilx
			end
		end

		for j = 3:n-2  # zuerst die spalten. ist etwas schneller
			for i = 3:m-2
				#vyph	= (v[i+1,j,t-1] + v[i,j,t-1])/2
				#vymh	= (v[i-1,j,t-1] + v[i,j,t-1])/2

				vyph	= v[i,  j,t-1] 
				vymh	= v[i-1,j,t-1] 

				anteily = fluss_lim_kons( vyph, p[i-1,j,t], p[i,j,t], p[i+1,j,t], p[i+2,j,t]) - fluss_lim_kons( vymh, p[i-2,j, t], p[i-1,j,t], p[i,j,t], p[i+1,j,t])

				p[i,j,t-1] = ph[i,j] - r* anteily
			end
		end
	end

	return p
end

