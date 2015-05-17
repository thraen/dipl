#using Images
#using ImageView

@everywhere using PyCall
#pygui(:gtk)
#pygui(:tk)

@everywhere using PyPlot
@everywhere pygui(false)

@everywhere isuff=".png"
#suff=".svg"
@everywhere vsuff=".dlm"

run(`mkdir -p $rootdir`)

function delete_output()
	run( `find $rootdir -name *$isuff -delete` )
end

@everywhere function save_image(im, pref, t)
	imshow(im[:,:,t], interpolation="none", origin="lower")
	#surf(im[:,:,t]',rstride=1, cstride=1) # rstride, cstride aufloesung
	savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	clf()
end

@everywhere function save_quiver(u,v, pref, t)
	quiver(u[:,:,t], v[:,:,t])
	savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	clf()
end

function save_quivers_(u,v, pref)
	tic()
	m, n, T = size(u)
	run(`mkdir -p $rootdir/$pref`)

	@sync @parallel for t=1:T 
		print("$t ")
		save_quiver(u, v, pref, t)
		#savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	end
	toc()
end

function save_images_(im, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)

	@sync @parallel for t=1:T 
		save_image(im, pref, t)
		#savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	end
	toc()
end

function save_images(im, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)
	for t=1:T 
		save_image(im, pref, t)
	end
	toc()
end

function save_all()
	save_images_(s, "s")
	save_images_(I, "I")
	#save_images_(p, "p")
	save_quivers_(u,v,"w")
	#save_quivers_(grd_u_J, grd_v_J,"grad_J")
	#save_w()
end

function save_value(M, pref)
	fn	= rootdir * pref * vsuff
	writedlm(fn, M)
end

function save_values(M, pref)
	_, _, T	= size(M)
	pref = rootdir * pref * "/"
	run(`mkdir -p $pref`)
	for t = 1:T
		#fn	= pref * lpad(t, 8,"0") * vsuff
		fn	= "$(pref)$(t)$vsuff"
		writedlm(fn, M[:,:,t])
	end
end

function report(I, s, L2, H1, J, alpha, norm_s)
end


#x = linspace(0,2*pi,1000); y = sin(3*x + 4*cos(2*x));
#plot(x, y, color="red", linewidth=2.0, linestyle="--")
#title("A sinusoidally modulated sinusoid")

