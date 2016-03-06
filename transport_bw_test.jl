include("transport_rand_bw.jl")

@everywhere @inline function procchunk_x_bw_center!(p, ph, u, t, irange, jrange)
	for j = jrange
		for i = irange
# 			@inbounds uxph	= (u[i,j,t-1] + u[i,j+1,t-1])/2
# 			@inbounds uxmh	= (u[i,j,t-1] + u[i,j-1,t-1])/2

			@inbounds uxph	= u[i,j+0.5,t-1]
			@inbounds uxmh	= u[i,j-0.5,t-1] 
			@inbounds anteilx = fluss_lim_kons( uxph, p[i,j-1,t], p[i,j,t], p[i,j+1,t], p[i,j+2,t]) - fluss_lim_kons( uxmh, p[i,j-2,t], p[i,j-1,t], p[i,j,t], p[i,j+1,t])
			@inbounds ph[i,j] = p[i,j,t] - r* anteilx
		end
	end
end

@everywhere @inline function procchunk_y_bw_center!(p, ph, v, t, irange, jrange)
	for j = jrange
		for i = irange
# 			@inbounds vyph		= (v[i+1,j,t-1] + v[i,j,t-1])/2
# 			@inbounds vymh		= (v[i-1,j,t-1] + v[i,j,t-1])/2

			@inbounds vyph		= v[i+0.5,j,t-1]
			@inbounds vymh		= v[i-0.5,j,t-1]

			@inbounds anteily		= fluss_lim_kons( vyph, ph[i-1,j], ph[i,j], ph[i+1,j], ph[i+2,j]) - fluss_lim_kons( vymh, ph[i-2,j], ph[i-1,j], ph[i,j], ph[i+1,j])
			@inbounds p[i,j,t-1]	= ph[i,j] - r* anteily
		end
	end
end

@everywhere @inline function procchunk_x_bw_interf!(p, ph, u, t, irange, jrange)
	@inbounds begin
	for j = jrange
		for i = irange
			@inbounds uxph	= u[i,j,t-1]
			@inbounds uxmh	= u[i,j-1,t-1]
			@inbounds anteilx = fluss_lim_kons( uxph, p[i,j-1,t], p[i,j,t], p[i,j+1,t], p[i,j+2,t]) - fluss_lim_kons( uxmh, p[i,j-2,t], p[i,j-1,t], p[i,j,t], p[i,j+1,t])

			@inbounds ph[i,j] = p[i,j,t] - r* anteilx
		end
	end
	end
end

@everywhere @inline function procchunk_y_bw_interf!(p, ph, v, t, irange, jrange)
	#@inbounds begin
	for j = jrange
		for i = irange
			@inbounds vyph	= v[i,  j,t-1] 
			@inbounds vymh	= v[i-1,j,t-1] 
			@inbounds anteily = fluss_lim_kons( vyph, ph[i-1,j], ph[i,j], ph[i+1,j], ph[i+2,j]) - fluss_lim_kons( vymh, ph[i-2,j], ph[i-1,j], ph[i,j], ph[i+1,j])

			@inbounds p[i,j,t-1] = ph[i,j] - r* anteily
		end
	end
	#end
end

function ruecktransport_ser!(s, I, u, v, n_samp, n_zsamp, norm_s)
	m, n, T = size(I)
	p		= zeros(m,n,T)
	ph		= zeros(m,n)

	uip = interpolate(u, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear())), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus
	vip = interpolate(v, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear())), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus

# 	uip = interpolate(u, (BSpline(Cubic(Flat())), BSpline(Cubic(Flat())), BSpline(Cubic(Flat()))), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus
# 	vip = interpolate(v, (BSpline(Cubic(Flat())), BSpline(Cubic(Flat())), BSpline(Cubic(Flat()))), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus

	println("==============Ruecktransport============$n x $m x $n_samples $n_zwischensamples, $T")
	sk		= 0
	for t = T:-1:2
		# thr!!! rechter oder linker grenzwert zum diracterm?
		if mod(t-1, n_zsamp+1) == 0 then
# 			echo("am sample        ", t, "->", t-1, n_samp-sk)
			err			= I[:, :, t] - s[:, :, n_samp-sk] 
			p[:,:,t] 	= p[:,:,t] - err/norm_s
			ph			= p[:,:,t]
			sk 			+= 1
		else
# 			echo("zwischensample        ", t, "->", t-1, n_samp-sk)
		end
		procchunk_x_bw!(p, ph, uip, t, 3:m-2, 3:n-2 )
		procchunk_y_bw!(p, ph, vip, t, 3:m-2, 3:n-2 )
	end
	return p
end
#thr!
ruecktransport_par! = ruecktransport_ser!
