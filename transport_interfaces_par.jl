#@everywhere function procchunk_x!(I::SharedArray, Ih::SharedArray, u::SharedArray, t, irange, jrange)
#@everywhere function procchunk_x!(I, Ih::SharedArray, u, t, irange, jrange)
@everywhere function procchunk_x!(I, Ih, u, t, irange, jrange)
	#println(irange, " ",jrange)
	for j = jrange
		for i = irange
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
end

@everywhere function procchunk_y!(I, Ih, v, t, irange, jrange)
	#println(irange, " ",jrange)
	for j = jrange
		for i = irange
			vyph	= v[i,j,t]
			vymh	= v[i-1,j,t]

			vymh_p	= max( vymh, 0)
			vyph_m	= min( vyph, 0)

			Wmh		= Ih[i,j]   - Ih[i-1,j]
			Wph		= Ih[i+1,j] - Ih[i,j]

			anteil_low		= vyph_m*Wph + vymh_p*Wmh

			anteil_higmh	= limited_hot( vymh, Ih[i-2,j], Ih[i-1,j], Ih[i,j], Ih[i+1,j] )
			anteil_higph	= limited_hot( vyph, Ih[i-1,j], Ih[i,j], Ih[i+1,j], Ih[i+2,j] )

			I[i,j,t+1] = Ih[i,j] - r* anteil_low - r*(anteil_higph - anteil_higmh)
		end
	end
end

@everywhere function range_part_x(I::SharedArray)
    idx = indexpids(I)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(I))
    splits = [round(Int, s) for s in linspace(2,size(I,2)-2,nchunks+1)]
    return 1:size(I,1), splits[idx]+1:splits[idx+1]
end

@everywhere function range_part_y(I::SharedArray)
    idx = indexpids(I)
    if idx == 0
        # This worker is not assigned a piece
        return 1:0, 1:0
    end
    nchunks = length(procs(I))
	# thr die Randbedinungen sind hier anders!
    splits = [round(Int, s) for s in linspace(0,size(I,2),nchunks+1)]
    return 3:size(I,1)-2, splits[idx]+1:splits[idx+1]
end

#@everywhere proc_chunks_shared_x!(I::SharedArray, Ih::SharedArray, u::SharedArray, t) = procchunk!(I, Ih, u, t, range_part_x(Ih)...)
@everywhere proc_chunks_shared_x!(I, Ih, u, t) = procchunk_x!(I, Ih, u, t, range_part_x(Ih)...)

@everywhere proc_chunks_shared_y!(I, Ih, v, t) = procchunk_y!(I, Ih, v, t, range_part_y(Ih)...)

function transport_par(I0, u, v, schritte)
	m, n	= size(I0)
	#I		= zeros(m, n, schritte+1)
	I		= SharedArray(Float64, (m,n,schritte+1))
	#Ih		= zeros(m, n) #half step buffer for dimension splitting
	I[:,:,1]= I0
	
	Ih		= SharedArray(Float64, (m,n))

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

