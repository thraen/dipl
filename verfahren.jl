using HDF5, JLD
include("echo.jl")

include("transport.jl")

include("matrizen.jl")
const B		= generate_L2(m, n, dx)

include("misc.jl")
L2norm		= function(s) return Xnorm(s, B) end # thr das sollte besser L2sum heissen
# L2norm = l2norm #riemannsummen auf den zellen

sample_err	= sample_err_L2

time_regularization && velocities_at =="centers" 		&& include("grad_time_reg.jl")  
time_regularization && velocities_at =="interfaces" 	&& include("grad_time_reg_interfaces.jl")  

~time_regularization && velocities_at == "centers"		&& include("grad_centers.jl") 	
# ~time_regularization && velocities_at == "centers"		&& include("grad_centers_inner.jl") 	

~time_regularization && velocities_at == "interfaces"	&& include("grad_interfaces.jl")

function plot_grad_section(lastJ, u, v, grd_u_J, grd_v_J, H1_J_w, s0, norm_s, k, steps)
	Js		= zeros(k)
	stepws	= zeros(k)
	for exp in 0:k-1
		t 				= armijo_bas^exp
		u_next			= u - t*grd_u_J
		v_next			= v - t*grd_v_J
		H1_err_next		= H1_norm_w(u_next, v_next)
		I_next			= transport( s0, u_next, v_next, T-1 )
		L2_err_next, _	= sample_err(I_next,s,norm_s)
		stepws[exp+1]	= t
		Js[exp+1]		= (L2_err_next + H1_err_next)/2 
	end
	clf()
	rng=[1:k]
	plot(stepws[rng], [lastJ for exp in 1:k][rng], color="black" )
	plot(stepws[rng], Js[rng], "bx")
	plot(stepws[rng], [lastJ-t*H1_J_w for t in stepws][rng], color="green")
	savefig(rootdir * "/gradientenrichtung" * lpad(steps, 8,"0") * isuff, dpi=dpi)
# 	plot(stepws[rng], [lastJ-0.3*H1_J_w*t for t in stepws][rng])
# 	plot(stepws[rng], [lastJ-armijo_sig*H1_J_w*t for t in stepws][rng])
end

function next_w!(I, p, u, v, alpha)
	for t= 1:T-1
		pI_x_		= reshape(Cx*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		pI_y_		= reshape(Cy*reshape(I[:,:,t], n*m) , m, n).*p[:,:,t]
		u[:,:,t]	= poissolv( -pI_x_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) ) /alpha
		v[:,:,t]	= poissolv( -pI_y_[2:m-1,2:n-1], zeros(1,n), zeros(1,n), zeros(m-2), zeros(m-2) ) /alpha
	end
	return u, v
end

function verfahren_grad(s, u, v, steps=1, normierung=1.0, grad_bound=1e-8)
	s0			= s[:,:,1]
	norm_s		= normierung

	H1_err		= H1_norm_w( u, v )
	H1_errs		= [H1_err]

	I	= zeros(m,n,T)
	for t = 1:T I[:,:,t] = s0 end

	@time I		= transport!(I, u, v, T-1)
	@time p		= ruecktransport!( s, I, -u, -v, n_samples, n_zwischensamples, norm_s )
	L2_err, _	= sample_err(I,s,norm_s)
	L2_errs		= [L2_err]

	echo("START $n x $m x $T ($n_samples samples x $n_zwischensamples zwischsamples), dx = $dx, dt=$dt, alpha=$alpha, beta=$beta",
		 "\ninitial L2_err", L2_err)

	@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

	@show H1_J_w= H1_norm_grd(grd_u_J, grd_v_J)
	H1_J_ws	= [H1_J_w]

	J		= (L2_err + H1_err)/2
	Js		= [J]

	# Armijo-Schrittweite
	armijo_exp	= 0

	while steps < maxsteps  &&  armijo_exp < armijo_maxtry  &&  H1_J_w > grad_bound 
# 		if steps % 5 == 0
# 			plot_grad_section(J, u, v, grd_u_J, grd_v_J, H1_J_w, s0, norm_s, 140, steps)
# 		end
		
		while (armijo_exp < armijo_maxtry)
			t 					= armijo_bas^armijo_exp

			u_next				= u - t*grd_u_J
			v_next				= v - t*grd_v_J

			H1_err_next			= H1_norm_w(u_next, v_next)

			@time I_next		= transport!( I, u_next, v_next, T-1 )
			L2_err_next, _		= sample_err(I_next,s,norm_s)

			J_next				= (L2_err_next + H1_err_next)/2 

			echo("\nstep", steps, armijo_exp,"test armijo step length ", t, 
				 "\nL2errors ",   L2_err, L2_err_next, L2_err-L2_err_next, 
				 "\nH1_errors",   H1_err, H1_err_next, H1_err-H1_err_next,
				 "\nJ        ",   J, J_next,J-J_next,"\n")

			if J_next < (J - armijo_sig * t *H1_J_w)
				I					= I_next
				u					= u_next
				v					= v_next

				H1_err				= H1_err_next
				push!(H1_errs, H1_err)
				L2_err				= L2_err_next
				push!(L2_errs, L2_err)

				@time p					= ruecktransport!(s, I, -u, -v, n_samples, n_zwischensamples, norm_s)
				@time grd_u_J, grd_v_J	= grad_J(I, p, u, v)

				H1_J_w					= H1_norm_grd(grd_u_J, grd_v_J)
				push!(H1_J_ws, H1_J_w)

				J					= (L2_err + H1_err)/2
				push!(Js, J)

				armijo_exp = 0
				echo("\n****** NEW GRADIENT *****\n", 
					 "max abs grd_u_J", maximum(abs(grd_u_J)),
					 "max abs grd_v_J", maximum(abs(grd_v_J)), "\n",
					 "H1_J_w", H1_J_w, "\n")
				break 
			end
			
			armijo_exp += 1
		end

		if (save_every > 0) && (steps % save_every == 0)
			save_jld(steps, dx, dt, alpha, beta, s, I, p, u, v, grd_u_J, grd_v_J)
		end

		steps%50 == 0 && @everywhere gc()

		steps +=1
	end

	return I, u, v, p, L2_errs, H1_errs, Js, H1_J_ws, steps
end
