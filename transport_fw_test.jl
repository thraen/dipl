@everywhere macro cfl_verletzt(v)
	cflcheck =  quote
		if abs($v)*r > 1.0
			info("cfl verletzt, $i, $j, $t, " * string(v[i,j,t]))
		end
	end
	return with_cfl_check ? cflcheck : nothing
end

include("transport_rand_fw.jl")

@everywhere @inline function divergence(um, up, h)
	return (up-um)/h
end

#@everywhere function procchunk_x_fw_interf!(I::SharedArray, Ih::SharedArray, u::SharedArray, t, irange, jrange)
@everywhere @inline function procchunk_x_fw_interf!(I, Ih, u, t, irange, jrange)
	#@inbounds begin
	for j = jrange
		for i = irange
			#@inbounds uxph	= u[i,j,t] 
			#@inbounds uxmh	= u[i,j-1,t]

			#das einzige, was ein bisschen was bring mit cubic interpolation
			uxph	= u[i,j-0.5,t]
			uxmh	= u[i,j-0.5,t]
			
			@inbounds uxph_m	= min( uxph, 0)
			@inbounds uxmh_p	= max( uxmh, 0)

			@inbounds Wmh		= I[i,j,t]   - I[i,j-1,t]
			@inbounds Wph		= I[i,j+1,t] - I[i,j,t]

			@inbounds anteil_higmh	= limited_hot( uxmh, I[i,j-2,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] )
			@inbounds anteil_higph	= limited_hot( uxph, I[i,j-1,t], I[i,j,t], I[i,j+1,t], I[i,j+2,t] )

			@inbounds anteil_low	= uxph_m*Wph + uxmh_p*Wmh


			@inbounds Ih[i,j] = I[i,j,t] - r* anteil_low - r*(anteil_higph - anteil_higmh)

# 			u_	= u[i,j-0.5,t]
# 			@inbounds anteilx = (u_ >= 0) 	? fluss_lim1( u_, I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u_, I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) : fluss_lim1( -u_, I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u_, I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
# 			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
	#end
end

@everywhere @inline function procchunk_y_fw_interf!(I, Ih, v, t, irange, jrange)
	#@inbounds begin
	for j = jrange
		for i = irange
			#@inbounds vyph	= v[i,j,t]
			#@inbounds vymh	= v[i-1,j,t]

			@inbounds vyph	= v[i-0.5,j,t]
			@inbounds vymh	= v[i-0.5,j,t]

			@inbounds vymh_p	= max( vymh, 0)
			@inbounds vyph_m	= min( vyph, 0)

			@inbounds Wmh		= Ih[i,j]   - Ih[i-1,j]
			@inbounds Wph		= Ih[i+1,j] - Ih[i,j]

			@inbounds anteil_low		= vyph_m*Wph + vymh_p*Wmh

			@inbounds anteil_higmh	= limited_hot( vymh, Ih[i-2,j], Ih[i-1,j], Ih[i,j], Ih[i+1,j] )
			@inbounds anteil_higph	= limited_hot( vyph, Ih[i-1,j], Ih[i,j], Ih[i+1,j], Ih[i+2,j] )

			@inbounds I[i,j,t+1] = Ih[i,j] - r* anteil_low - r*(anteil_higph - anteil_higmh)

# 			v_	= v[i-0.5,j,t]
# 			@inbounds anteily = (v_ >= 0) 	? fluss_lim1(  v_, Ih[i-1,j], Ih[i,j], Ih[i+1,j] ) - fluss_lim1(  v_, Ih[i-2,j], Ih[i-1,j], Ih[i,j] ) : fluss_lim1( -v_, Ih[i+1,j], Ih[i,j], Ih[i-1,j] ) - fluss_lim1( -v_, Ih[i+2,j], Ih[i+1,j], Ih[i,j] ) 
# 			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
	#end
end

@everywhere @inline function procchunk_x_fw_center!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
# 			u_	= u[i,j,t]
			u_	= u[i,j,t+0.5]
			@cfl_verletzt u_
			@inbounds anteilx = (u_ >= 0) 	? fluss_lim1( u_, I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u_, I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) : fluss_lim1( -u_, I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u_, I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end

@everywhere @inline function procchunk_y_fw_center!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
# 			v_	= v[i,j,t]
			v_	= v[i,j,t+0.5]
			@cfl_verletzt v_
			@inbounds anteily = (v_ >= 0) 	? fluss_lim1(  v_, Ih[i-1,j], Ih[i,j], Ih[i+1,j] ) - fluss_lim1(  v_, Ih[i-2,j], Ih[i-1,j], Ih[i,j] ) : fluss_lim1( -v_, Ih[i+1,j], Ih[i,j], Ih[i-1,j] ) - fluss_lim1( -v_, Ih[i+2,j], Ih[i+1,j], Ih[i,j] ) 
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

@everywhere proc_chunks_shared_x!(I, Ih, u, t) = procchunk_x_fw!(I, Ih, u, t, range_part_x(Ih)...)
@everywhere proc_chunks_shared_y!(I, Ih, v, t) = procchunk_y_fw!(I, Ih, v, t, range_part_y(Ih)...)

#instabil
# Rand ist hier nicht implementiert
function transport_par(I0, u, v, schritte)
	@time begin
	m, n	= size(I0)
	I		= SharedArray(Float64, (m,n,schritte+1))
	I[:,:,1]= I0
	Ih		= SharedArray(Float64, (m,n))

	u = convert(SharedArray{Float64}, u)
	v = convert(SharedArray{Float64}, v)
	end
	println("==============Transport=================$n x $m x $T")
	for t = 1:schritte
		@sync begin
			for p in procs(Ih)
				@async remotecall_wait( p, proc_chunks_shared_x!, I, Ih, u, t )
			end
		end

		@sync begin
			for p in procs(Ih)
				@async remotecall_wait( p, proc_chunks_shared_y!, I, Ih, v, t )
			end
		end
	end
	return I
end

using Interpolations

function transport_ser!(I, u, v, schritte)
	m, n, T	= size(I)
	Ih		= I[:,:,1] #half step buffer for dimension splitting
	println("==============Transport=================$n x $m x $T parallel=$(transport_parallel)")

# 	uip = interpolate(u, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear())), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus
# 	vip = interpolate(v, (BSpline(Linear()), BSpline(Linear()), BSpline(Linear())), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus

	#Flat, Line (alternatively, Natural), Free, Periodic and Reflect
# 	uip = interpolate(u, (BSpline(Quadratic(Flat())), BSpline(Quadratic(Flat())), BSpline(Quadratic(Flat()))), OnGrid()) #thr ongrid? oncell?
# 	vip = interpolate(v, (BSpline(Quadratic(Flat())), BSpline(Quadratic(Flat())), BSpline(Quadratic(Flat()))), OnGrid()) #thr ongrid? oncell?

	uip = interpolate(u, (BSpline(Cubic(Flat())), BSpline(Cubic(Flat())), BSpline(Cubic(Flat()))), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus
	vip = interpolate(v, (BSpline(Cubic(Flat())), BSpline(Cubic(Flat())), BSpline(Cubic(Flat()))), OnGrid()) #thr ongrid? oncell? kommt das gleiche raus

	println("FUCKFUCK")
	println("FUCKFUCK")
	println("vorsicht testversion")

	for t = 1:schritte
		procchunk_x_fw!(I, Ih, uip, t, 3:m-2, 3:n-2 )
		procchunk_y_fw!(I, Ih, vip, t, 3:m-2, 3:n-2 )
	end
	return I
end

