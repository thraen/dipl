
@everywhere @inline function procchunk_x_fw_center_innerer_rand_LR!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			u_	= @interpolate_w_time_center_fw(u)
			# nur Upwind
			if (u_ >=0)
				anteilx = u_*( I[i,j,t] - I[i,j-1,t] )
			else
				anteilx = u_*( I[i,j+1,t] - I[i,j,t] ) # u<0
			end
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end

@everywhere @inline function procchunk_y_fw_center_innerer_rand_OU!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			#nur Upwind
			v_	= @interpolate_w_time_center_fw(v)
			if (v_ >= 0) 	
				anteily = v_*(Ih[i,j] - Ih[i-1,j])
			else
				anteily = v_*(Ih[i+1,j] - Ih[i,j])
			end
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

@everywhere @inline function procchunk_x_fw_interf_innerer_rand_LR!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			@inbounds uxph	= u[i,j-0.5,t]
			@inbounds uxmh	= u[i,j-0.5,t]
			
			@inbounds uxph_m	= min( uxph, 0)
			@inbounds uxmh_p	= max( uxmh, 0)

			@inbounds Wmh		= I[i,j,t]   - I[i,j-1,t]
			@inbounds Wph		= I[i,j+1,t] - I[i,j,t]
			# nur upwind
			@inbounds Ih[i,j] = I[i,j,t] - r*( uxph_m*Wph + uxmh_p*Wmh )
		end
	end
end

@everywhere @inline function procchunk_y_fw_interf_innerer_rand_OU!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			@inbounds vyph	= v[i-0.5,j,t]
			@inbounds vymh	= v[i-0.5,j,t]

			@inbounds vymh_p	= max( vymh, 0)
			@inbounds vyph_m	= min( vyph, 0)

			@inbounds Wmh		= Ih[i,j]   - Ih[i-1,j]
			@inbounds Wph		= Ih[i+1,j] - Ih[i,j]

			@inbounds I[i,j,t+1] = Ih[i,j] - r*( vyph_m*Wph + vymh_p*Wmh )
		end
	end
end
