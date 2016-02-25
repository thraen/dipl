#@everywhere function procchunk_x_fw_interf!(I::SharedArray, Ih::SharedArray, u::SharedArray, t, irange, jrange)
@everywhere @inline function procchunk_x_fw_interf!(I, Ih, u, t, irange, jrange)
	#@inbounds begin
	for j = jrange
		for i = irange
			@inbounds uxph	= u[i,j,t] 
			# thr! was passiert, wenn ich hier auch u[i,j,t] nehme? oder den Durchschnitt und dann wie center?
			@inbounds uxmh	= u[i,j-1,t]
			
			@inbounds uxph_m	= min( uxph, 0)
			@inbounds uxmh_p	= max( uxmh, 0)

			@inbounds Wmh		= I[i,j,t]   - I[i,j-1,t]
			@inbounds Wph		= I[i,j+1,t] - I[i,j,t]

			@inbounds anteil_higmh	= limited_hot( uxmh, I[i,j-2,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] )
			@inbounds anteil_higph	= limited_hot( uxph, I[i,j-1,t], I[i,j,t], I[i,j+1,t], I[i,j+2,t] )

			@inbounds anteil_low	= uxph_m*Wph + uxmh_p*Wmh

			@inbounds Ih[i,j] = I[i,j,t] - r* anteil_low - r*(anteil_higph - anteil_higmh)
		end
	end
	#end
end

@everywhere @inline function procchunk_y_fw_interf!(I, Ih, v, t, irange, jrange)
	#@inbounds begin
	for j = jrange
		for i = irange
			@inbounds vyph	= v[i,j,t]
			# thr! was passiert, wenn ich hier auch v[i,j,t] nehme?
			@inbounds vymh	= v[i-1,j,t]

			@inbounds vymh_p	= max( vymh, 0)
			@inbounds vyph_m	= min( vyph, 0)

# 			if vyph*r > 0.05 || vymh < -.05
# 				error("cfl verletzt, $i, $j")
# 			end

			@inbounds Wmh		= Ih[i,j]   - Ih[i-1,j]
			@inbounds Wph		= Ih[i+1,j] - Ih[i,j]

			@inbounds anteil_low		= vyph_m*Wph + vymh_p*Wmh

			@inbounds anteil_higmh	= limited_hot( vymh, Ih[i-2,j], Ih[i-1,j], Ih[i,j], Ih[i+1,j] )
			@inbounds anteil_higph	= limited_hot( vyph, Ih[i-1,j], Ih[i,j], Ih[i+1,j], Ih[i+2,j] )

			@inbounds I[i,j,t+1] = Ih[i,j] - r* anteil_low - r*(anteil_higph - anteil_higmh)
		end
	end
	#end
end

@everywhere @inline function procchunk_x_fw_center!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			@inbounds anteilx = (u[i,j,t] >= 0) 	? fluss_lim1(  u[i,j,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u[i,j,t], I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) : fluss_lim1( -u[i,j,t], I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u[i,j,t], I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end

@everywhere macro cfl_verletzt(v)
	cflcheck =  quote
		if abs($v)*r > 1.0
			info("cfl verletzt, $i, $j, $t, " * string(v[i,j,t]))
		end
	end
	return with_cfl_check ? cflcheck : nothing
end

@everywhere @inline function procchunk_y_fw_center!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			@cfl_verletzt v[i,j,t]
			@inbounds anteily = (v[i,j,t] >= 0) 	? fluss_lim1(  v[i,j,t], Ih[i-1,j], Ih[i,j], Ih[i+1,j] ) - fluss_lim1(  v[i,j,t], Ih[i-2,j], Ih[i-1,j], Ih[i,j] ) : fluss_lim1( -v[i,j,t], Ih[i+1,j], Ih[i,j], Ih[i-1,j] ) - fluss_lim1( -v[i,j,t], Ih[i+2,j], Ih[i+1,j], Ih[i,j] ) 
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

@everywhere proc_chunks_shared_x!(I, Ih, u, t) = procchunk_x_fw!(I, Ih, u, t, range_part_x(Ih)...)
@everywhere proc_chunks_shared_y!(I, Ih, v, t) = procchunk_y_fw!(I, Ih, v, t, range_part_y(Ih)...)

#instabil
function transport_par(I0, u, v, schritte)
	@time begin
	m, n	= size(I0)
	I		= SharedArray(Float64, (m,n,schritte+1))
	I[:,:,1]= I0
	Ih		= SharedArray(Float64, (m,n))

	u = convert(SharedArray{Float64}, u)
	v = convert(SharedArray{Float64}, v)
	end

# 	@show typeof(u)

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

function transport_ser(I0, u, v, schritte)
	m, n	= size(I0)
	I		= zeros(m, n, schritte+1)
	I[:,:,1]= I0
	Ih		= zeros(m, n) #half step buffer for dimension splitting
	println("==============Transport=================$n x $m x $T parallel=$(transport_parallel)")
	for t = 1:schritte
		procchunk_x_fw!(I, Ih, u, t, 3:m-2, 3:n-2 )
		procchunk_y_fw!(I, Ih, v, t, 3:m-2, 3:n-2 )
	end
	return I
end

