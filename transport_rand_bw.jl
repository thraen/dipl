@everywhere @inline function procchunk_x_bw_center_innerer_rand_LR!(p, ph, u, t, irange, jrange)
	for j = jrange
		for i = irange
			# wir machen hier nur Upwind. Man koennte wieder in Upwind-Richtung und nicht Upwind-Richtung unterscheiden, aber das ist zu viel Aufwand fuer Zero-Nutzen
			uxph	= (u[i,j,t-1] + u[i,j+1,t-1])/2
			uxmh	= (u[i,j,t-1] + u[i,j-1,t-1])/2
			anteilx = fluss_upw1( uxph, p[i,j,t], p[i,j+1,t] ) - fluss_upw1( uxmh, p[i,j-1,t], p[i,j,t] )
			ph[i,j] = p[i,j,t] - r* anteilx
		end
	end
end

@everywhere @inline function procchunk_y_bw_center_innerer_rand_OU!(p, ph, v, t, irange, jrange)
	for j = jrange
		for i = irange
			vyph		= (v[i+1,j,t-1] + v[i,j,t-1])/2
			vymh		= (v[i-1,j,t-1] + v[i,j,t-1])/2
			anteily		= fluss_upw1( vyph, ph[i,j], ph[i+1,j]) - fluss_upw1( vymh, ph[i-1,j], ph[i,j])
			p[i,j,t-1]	= ph[i,j] - r* anteily
		end
	end
end

@everywhere @inline function procchunk_x_bw_interf_innerer_rand_LR!(p, ph, u, t, irange, jrange)
	for j = jrange
		for i = irange
			@inbounds uxph	= u[i,j,t-1]
			@inbounds uxmh	= u[i,j-1,t-1]
			@inbounds anteilx = fluss_upw1( uxph, p[i,j,t], p[i,j+1,t]) - fluss_upw1( uxmh, p[i,j-1,t], p[i,j,t])

			@inbounds ph[i,j] = p[i,j,t] - r* anteilx
		end
	end
end

@everywhere @inline function procchunk_y_bw_interf_innerer_rand_OU!(p, ph, v, t, irange, jrange)
	for j = jrange
		for i = irange
			@inbounds vyph	= v[i,  j,t-1] 
			@inbounds vymh	= v[i-1,j,t-1] 
			@inbounds anteily = fluss_upw1( vyph, ph[i,j], ph[i+1,j]) - fluss_upw1( vymh, ph[i-1,j], ph[i,j])

			@inbounds p[i,j,t-1] = ph[i,j] - r* anteily
		end
	end
end
