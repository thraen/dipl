#using Images
#using ImageView

#@everywhere using PyCall
#pygui(:gtk)
#pygui(:tk)

@everywhere using PyPlot
@everywhere pygui(false)

@everywhere isuff=".png"
#suff=".svg"
@everywhere vsuff=".dlm"

function delete_output()
	run(`rm $rootdir -r`)
end

@everywhere function save_image(im, pref, t)
	imshow(im[:,:,t], interpolation="none", origin="lower")#, cmap=gray)
	#imshow(im[:,:,t], interpolation="none", cmap="gray")
	savefig(rootdir * pref * "/img" * lpad(t, 8,"0") * isuff)
	clf()
	surf(im[:,:,t]',rstride=1, cstride=1) # rstride, cstride aufloesung
	savefig(rootdir * pref * "/srf" * lpad(t, 8,"0") * isuff)
	clf()
end

@everywhere function save_quiver(u,v, pref, t, mpad, npad)
	quiver( [zeros(m, npad) u[:,:,t]], [zeros(mpad, n); v[:,:,t]] )
	savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	clf()
end

function save_quivers_(u,v, pref)
	tic()
	run(`mkdir -p $rootdir`)
	mu, nu, T = size(u)
	mv, nv, T = size(v)
	mpad = abs(mu-mv)
	npad = abs(nu-nv)
	run(`mkdir -p $rootdir/$pref`)

	@sync @parallel for t=1:T 
		print("$t ")
		save_quiver(u, v, pref, t, mpad, npad)
	end
	toc()
end

function save_images_(im, pref)
	tic()
	run(`mkdir -p $rootdir`)
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)

	@sync @parallel for t=1:T 
		save_image(im, pref, t)
		#savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff)
	end
	toc()
end

function save_images(im, pref)
	run(`mkdir -p $rootdir`)
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

function extract_convergence_history()
# 	run( `sh -c "grep L2err  $(rootdir)/log | awk '{print $2}' "`)
	run( `sh -c "grep H1_J_  $(rootdir)/log | awk '{print \$2}' |uniq > $rootdir/H1hist"`)
	run( `sh -c "grep L2err  $(rootdir)/log | awk '{print \$2}' |uniq > $rootdir/L2hist"`)
	run( `sh -c "grep 'J\s'  $(rootdir)/log | awk '{print \$2}' |uniq > $rootdir/Jhist"`)

# 	L2hist	= readdlm("$rootdir/L2hist")
# 	H1hist	= readdlm("$rootdir/H1hist")
# 	Jhist	= readdlm("$rootdir/Jhist")
	return L2hist, H1hist, Jhist
end

# L2hist, H1hist, Jhist = extract_convergence_history()


#function report()
	#f		= open(rootdir *"report.txt")
	#write(f, 
#end
