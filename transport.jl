using Interpolations
@everywhere const r			= dt/dx # thr marcel hat hier nur /dx

@everywhere @inline function fluss_lxw1( a, u, v )
	return 0.5* ( a*(u+v) + r*a*a*(u-v) )
end

@everywhere @inline function fluss_upw1( a, u, v )
	return (a>=0) ? (a*u) : (a*v)
end

@everywhere @inline function theta( um, u, up )
	return (u - um +eps()) / (up - u +eps())
end

@everywhere @inline function sbee( thet )
	return max( 0.0, max( min(1.0, 2*thet), min(thet, 2.0) ) );
end

@everywhere @inline function fluss_lim1( a, um, u, up )
	thet = theta(um, u, up)
	ret = (1- sbee( thet ) )* fluss_upw1(a, u, up ) +     sbee( thet )  * fluss_lxw1(a, u, up )
	return 	  ret
end

@everywhere @inline function fluss_lim_kons( a, umm, um, u, up )
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

@everywhere @inline function limited_hot(a, umm, um, u, up)
	if a>=0 
		thet = (um - umm +eps()) / (u - um +eps())
	else
		thet = (up - u   +eps()) / (u - um +eps())
	end
	anteil_hig	= 0.5*abs(a) * (1-r*abs(a)) * sbee(thet) * (u-um)
	return 	  anteil_hig
end

@everywhere macro interpolate_w_time_center_fw(uv)
	ohne =  quote
		$(uv)[i,j,t]
	end
	mit =  quote
		$(uv)[i,j,t+0.5]
	end
	return interpolate_w_time ? mit : ohne
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

# include("transport_fw_test.jl")  #THR! MACH DAS WIEDER WEG!
# include("transport_bw_test.jl") 

include("transport_fw.jl")
include("transport_bw.jl")

velocities_at == "centers" && begin
	@everywhere procchunk_x_fw!	= procchunk_x_fw_center!
	@everywhere procchunk_y_fw!	= procchunk_y_fw_center!
	# hier wird nur Upwind gemacht
	@everywhere procchunk_x_fw_innerer_rand_LR!	= procchunk_x_fw_center_innerer_rand_LR!
	@everywhere procchunk_y_fw_innerer_rand_OU!	= procchunk_y_fw_center_innerer_rand_OU!

	@everywhere procchunk_x_bw!	= procchunk_x_bw_center!
	@everywhere procchunk_y_bw!	= procchunk_y_bw_center!
	@everywhere procchunk_x_bw_innerer_rand_LR!	= procchunk_x_bw_center_innerer_rand_LR!
	@everywhere procchunk_y_bw_innerer_rand_OU!	= procchunk_y_bw_center_innerer_rand_OU!
end

velocities_at == "interfaces" && begin
	@everywhere procchunk_x_fw!	= procchunk_x_fw_interf!
	@everywhere procchunk_y_fw!	= procchunk_y_fw_interf!
	@everywhere procchunk_x_fw_innerer_rand_LR!	= procchunk_x_fw_interf_innerer_rand_LR!
	@everywhere procchunk_y_fw_innerer_rand_OU!	= procchunk_y_fw_interf_innerer_rand_OU!

	@everywhere procchunk_x_bw!	= procchunk_x_bw_interf!
	@everywhere procchunk_y_bw!	= procchunk_y_bw_interf!
	@everywhere procchunk_x_bw_innerer_rand_LR!	= procchunk_x_bw_interf_innerer_rand_LR!
	@everywhere procchunk_y_bw_innerer_rand_OU!	= procchunk_y_bw_interf_innerer_rand_OU!
end

transport_parallel || begin
	transport!		= transport_ser!
	ruecktransport!	= ruecktransport_ser!
end

transport_parallel && begin
	transport!		= transport_par!
	ruecktransport!	= ruecktransport_par!
end

