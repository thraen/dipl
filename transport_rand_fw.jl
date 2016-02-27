@everywhere @inline function procchunk_x_fw_center_innerer_rand_LR!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			# nur Upwind
			if (u[i,j,t] >=0)
				anteilx = u[i,j,t]*( I[i,j,t] - I[i,j-1,t] )
			else
				anteilx = u[i,j,t]*( I[i,j+1,t] - I[i,j,t] ) # u<0
			end
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end

@everywhere @inline function procchunk_x_fw_center_innerer_rand_R!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			if (u[i,j,t] >=0)
				# Bei Ausfluss, also u zeigt Richtung Rand kann man das normale Verfahren anwenden, alle Informationen liegen vor
				anteilx = fluss_lim1(  u[i,j,t], I[i,j-1,t], I[i,j,t], I[i,j+1,t] ) - fluss_lim1(  u[i,j,t], I[i,j-2,t], I[i,j-1,t], I[i,j,t] ) 
			else
				# Bei Einfluss wuerde das Laxwendroff-Verfahren einen Punkt außerhalb des Rand brauchen. Wir benutzten
				# an solchen Stelle daher ein Upwind-Verfahren. (Aequivalent dazu waere, den außerhalb liegenden Punkt auf den Randpunkt zu setzten => thet=0)
				anteilx = u[i,j,t]*( I[i,j+1,t] - I[i,j,t] ) # u<0
			end
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end
@everywhere @inline function procchunk_x_fw_center_innerer_rand_L!(I, Ih, u, t, irange, jrange)
	for j = jrange
		for i = irange
			if (u[i,j,t] >=0)
				# Bei Einfluss wuerde das Laxwendroff-Verfahren einen Punkt außerhalb des Rand brauchen. Wir benutzten
				# an solchen Stelle daher ein Upwind-Verfahren. (Aequivalent dazu waere, den außerhalb liegenden Punkt auf den Randpunkt zu setzten => thet=0)
				anteilx = u[i,j,t]*( I[i,j,t] - I[i,j-1,t] ) #u>0
			else
				# Bei Ausfluss, also u zeigt Richtung Rand kann man das normale Verfahren anwenden, alle Informationen liegen vor
				anteilx = fluss_lim1( -u[i,j,t], I[i,j+1,t], I[i,j,t], I[i,j-1,t] ) - fluss_lim1( -u[i,j,t], I[i,j+2,t], I[i,j+1,t], I[i,j,t] ) 
			end
			@inbounds Ih[i,j] = I[i,j,t] - r* (anteilx)
		end
	end
end

@everywhere @inline function procchunk_y_fw_center_innerer_rand_OU!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			#nur Upwind
			if (v[i,j,t] >= 0) 	
				anteily = v[i,j,t]*(Ih[i,j] - Ih[i-1,j])
			else
				anteily = v[i,j,t]*(Ih[i+1,j] - Ih[i,j])
			end
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

@everywhere @inline function procchunk_y_fw_center_innerer_rand_O!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			if (v[i,j,t] >= 0) 	
				# v zeigt nach aussen, alles gut
				anteily = fluss_lim1(  v[i,j,t], Ih[i-1,j], Ih[i,j], Ih[i+1,j] ) - fluss_lim1(  v[i,j,t], Ih[i-2,j], Ih[i-1,j], Ih[i,j] ) 
			else
				# v zeigt nach innen, upwind
				anteily = v[i,j,t]*(Ih[i+1,j] - Ih[i,j])
			end
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

@everywhere @inline function procchunk_y_fw_center_innerer_rand_U!(I, Ih, v, t, irange, jrange)
	for j = jrange
		for i = irange
			if (v[i,j,t] >= 0) 	
				# v zeigt nach innen, upwind
				anteily = v[i,j,t]*(Ih[i,j] - Ih[i-1,j])
			else
				# v zeigt nach aussen, alles gut
				anteily = fluss_lim1( -v[i,j,t], Ih[i+1,j], Ih[i,j], Ih[i-1,j] ) - fluss_lim1( -v[i,j,t], Ih[i+2,j], Ih[i+1,j], Ih[i,j] ) 
			end
			@inbounds I[i,j,t+1] = Ih[i,j] - r* (anteily)
		end
	end
end

# @everywhere @inline function procchunk_x_fw_center_rand!(I, Ih, u, t, irange, jrange)
# 	#hier ist u=0. man muss nichts machen
# 	for j = jrange
# 		for i = irange
# 			u[i,j,t] == 0 || println("FUCK u ist auf dem Rand nicht 0")
# 		end
# 	end
# end

# @everywhere @inline function procchunk_x_fw_center_rand_R!(I, Ih, u, t, irange, jrange)
# 	#hier ist u=0. man muss nichts machen
# 	for j = jrange
# 		for i = irange
# 			u[i,j,t] == 0 || println("FUCK u ist auf dem Rand nicht 0")
# 		end
# 	end
# end
# @everywhere @inline function procchunk_x_fw_center_rand_L!(I, Ih, u, t, irange, jrange)
# 	for j = jrange
# 		for i = irange
# 			u[i,j,t] == 0 || println("FUCK u ist auf dem Rand nicht 0")
# 		end
# 	end
# end
# @everywhere @inline function procchunk_y_fw_center_rand_O!(I, Ih, u, t, irange, jrange)
# 	#hier ist u=0. man muss nichts machen
# 	for j = jrange
# 		for i = irange
# 			u[i,j,t] == 0 || println("FUCK u ist auf dem Rand nicht 0")
# 		end
# 	end
# end
# @everywhere @inline function procchunk_y_fw_center_rand_U!(I, Ih, u, t, irange, jrange)
# 	#hier ist u=0. man muss nichts machen
# 	for j = jrange
# 		for i = irange
# 			u[i,j,t] == 0 || println("FUCK u ist auf dem Rand nicht 0")
# 		end
# 	end
# end
