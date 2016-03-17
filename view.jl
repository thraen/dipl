using PyPlot
@everywhere using PyPlot
@everywhere pygui(false)

@everywhere isuff=".png"
# @everywhere isuff=".eps"

# @everywhere isuff=".svg"
#suff=".svg"
@everywhere vsuff=".dlm"

# @everywhere dpi=1200
@everywhere dpi=60

function delete_output()
	run(`rm $rootdir -r`)
end

@everywhere function save_image_t(im, pref, t)
	imshow(im[:,:,t], interpolation="none", origin="lower", cmap="gray")
	#imshow(im[:,:,t], interpolation="none", cmap="gray")
	savefig(rootdir * pref * "/img" * lpad(t, 8,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

@everywhere function save_surf(im, pref, t)
	##
	## surf(X,Y,Z)!! man muss die arrays also transponieren!
	##
	surf(im[:,:,t]',rstride=1, cstride=1) # rstride, cstride aufloesung
	fig	= gcf()
	ax	= gca()
	ax[:set_zlim](0,1) # z-Achsen range auf 0,1 festlegen
	savefig(rootdir * pref * "/srf" * lpad(t, 8,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

@everywhere function save_quiver(u,v, pref, t, mpad, npad)
	quiver( [zeros(m, npad) u[:,:,t]], [zeros(mpad, n); v[:,:,t]] )
	savefig(rootdir * pref * "/" * lpad(t, 8,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

function save_quivers_(u,v, pref)
	tic()
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
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)
	@sync @parallel for t=1:T 
		save_image_t(im, pref, t)
	end
	toc()
end

function save_surfs_(im, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)
	@sync @parallel for t=1:T 
		save_surf(im, pref, t)
	end
	toc()
end

function save_images(im, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $rootdir/$pref`)
	for t=1:T 
		save_image_t(im, pref, t)
	end
	toc()
end

function save_all()
	println("save_images_ $rootdir") 
	save_images_(s, "s")

	save_images_(I, "I")
	save_surfs_(I,"I")

# 	save_images_(I_vorgabe, "I_given")
# 	save_surfs_(I_vorgabe, "I_given")

# 	save_images_(diff_vorgabe, "diff_vorgabe")
# 	save_surfs_(diff_vorgabe, "diff_vorgabe")

# 	save_images_(p, "p")
	save_quivers_(u,v,"w")
# 	save_quivers_(grd_u_J, grd_v_J,"grad_J")
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
	return L2hist, H1hist, Jhist
end
# L2hist, H1hist, Jhist = extract_convergence_history()

function make_output_dir(dir)
	run(`mkdir -p $dir/src`)
	run(`sh -c "cp *jl $dir/src"`)
	run(`sh -c "git log -1 > $dir/this_git_commit"`) 
end

function tabularline(l)
	lstr	=  string(l[1])
	for c in l[2:end]
		lstr *= " & $c"
	end
	return lstr *"\\\\\n"
end

function latextable_normal(caption, label, headl, lines...)
	tblhead		= "\\begin{table}[h]\n"
	captions	= "\\caption{$(caption)}\n"
	labels		= "\\label{table:$(label)}\n"
	hline		= "\\hline\n"
	nc			= length(headl)
	nl			= length(lines)

	tabularhead	= "\\begin{tabular}{l*{$(nc-1)}{c}}\n"
	
	headstr		= tabularline(headl)

	liness		= foldl(*, map(tabularline, lines))
	
	str	= tblhead * captions * labels * tabularhead *
			hline * headstr * hline *
			liness * hline *
			"\\end{tabular}\n\\end{table}"
	return str
end
# print(latextable_normal("capt", "lbl", ["h1", "h2", "h3"], [11, 12, 13], [21, 22, 23]))

function to_file(file, str)
	tblfile = open(file, "w")
	write(file, str); close(tblfile)
end


# Spezielle Reports
# 
# fuer die 2-Bilder-Demos
#
@everywhere function save_displacement(dir, isuff, dpi)
	clf()
	imshow(s[:,:,2]-s[:,:,1], cmap="gray_r", interpolation="none", origin="lower")
	savefig(dir * "displacement" *isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
end

# @everywhere cmpsurf = "gray_r"
@everywhere cmpsurf = "gray"
# @everywhere cmpsurf = "GnBu"
# @everywhere cmpsurf = "Blues"

@everywhere function save_auswahl_rot_disc(what, name, ts, dir, isuff, dpi, zlim=(0,1), azim=-30, elev=55)
	@sync @parallel for t in ts
		clf()
		imshow(what[:,:,t], cmap="gray_r", interpolation="none", origin="lower")
		savefig(dir * "img_$name" * lpad(Int(t), 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
		clf()

		surf(what[:,:,t]', cmap=cmpsurf, cstride=1, rstride=1, linewidth=0.2)
		fig	= gcf()
		ax	= gca()
		ax[:set_zlim](zlim) # z-Achsen range auf 0,1 festlegen
		ax[:view_init](azim=-39, elev=53)

		savefig(dir * "srf_$name" * lpad(Int(t), 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	end
end

@everywhere function save_verr(verr, name, dir, isuff, dpi)
	m,n,ts = size(verr)
	@sync @parallel for t in 1:ts
		imshow(abs(verr[:,:,t]), cmap="gray_r", interpolation="none", origin="lower")
		savefig(dir * "img_$name" * lpad(t, 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
		clf()
	end
end

function save_demo_rot_disc()
	save_displacement(rootdir, ".eps", 1200)
	save_displacement(rootdir, ".png", 100)
	vorgabe_frames	= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
	auswahl = collect(vorgabe_frames)
	save_auswahl_rot_disc(I, "I", auswahl, rootdir, ".png", 100)
	save_auswahl_rot_disc(I, "I", auswahl, rootdir, ".eps", 1200)
	save_verr(vorgabe_fehler, "verr", rootdir, ".eps", 1200)
	save_verr(vorgabe_fehler, "verr", rootdir, ".png", 100)
end

@everywhere function save_auswahl_taxi(what, name, ts, dir, isuff, dpi, zlim=(0,1), azim=-30, elev=55)
	@sync @parallel for t in ts
		clf()
		imshow(what[:,:,t], cmap="gray", interpolation="none", origin="lower")
		ax	= gca()
		ax[:set_axis_off]() 
		savefig(dir * "img_$name" * lpad(Int(t), 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)

		clf()
		surf(what[:,:,t]', cmap=cmpsurf, cstride=1, rstride=1, linewidth=0.2)
		fig	= gcf()
		ax[:set_zlim](zlim) # z-Achsen range auf 0,1 festlegen
		ax[:view_init](azim=-39, elev=53)
		savefig(dir * "srf_$name" * lpad(Int(t), 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	end
end

function save_demo_taxi()
	save_displacement(rootdir, ".eps", 1200)
	save_displacement(rootdir, ".png", 100)
	vorgabe_frames	= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
	auswahl = collect(vorgabe_frames)
	save_auswahl_taxi(I, "I", auswahl, rootdir, ".png", 100)
	save_auswahl_taxi(I, "I", auswahl, rootdir, ".eps", 1200)
	save_verr(vorgabe_fehler, "verr", rootdir, ".eps", 1200)
	save_verr(vorgabe_fehler, "verr", rootdir, ".png", 100)
end

function demo_table(capt, label)
# 	head	= ["\$\\alpha\$", "\$\\beta\$", "\$\\sum\\|I-s\\|_2^2\$", "Reg.-Fehler", "\$\\sum \\|V-I\\|_2^2\$"]
	head	= ["\$\\alpha\$"; "\$\\beta\$"; "\$\\e_{L^2}(I,s)\$"; "\$\\e_{reg_x}(w)\$"; "\$\\e_{reg_t}(w)\$"; "\$\\e_{L^2}(I,V)\$"]
	bet		= (time_regularization == false) ? 0 : beta
	res		= [alpha; bet; L2_errs[end]; H1_norm_w_noweight_space(u,v); H1_norm_w_noweight_time(u,v); L2norm(vorgabe_fehler)]
	to_file(rootdir*"table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )

	to_file(rootdir*"head_"*"errors"*".tex", tabularline(head))
	to_file(rootdir*"line_"*"errors"*".tex", tabularline(res))
	print( rootdir*"table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )
end

function nzw_table(capt, label)
	head	= ["\$# Zwischenbilder\$"; "\$\\e_{L^2}(I,s)\$"; "\$\\e_{reg_x}(w)\$"; "\$\\e_{L^2}(I,V)\$"]
	res		= [n_zwischensamples; L2_errs[end]; H1_norm_w_noweight_space(u,v); L2norm(vorgabe_fehler)]
	to_file(rootdir*"v1table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )
	to_file(rootdir*"v1head_"*"errors"*".tex", tabularline(head))
	to_file(rootdir*"v1line_"*"errors"*".tex", tabularline(res))
	#v2
	head	= ["\$# dummy frames\$"; "\$\\e_{L^2}(I,s)\$"; "\$\\e_{reg_x}(w)\$"; "\$\\e_{L^2}(I,V)\$"]
	res		= [zwischen_ausgelassen; L2_errs[end]; H1_norm_w_noweight_space(u,v); L2norm(vorgabe_fehler)]
	to_file(rootdir*"v2table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )
	to_file(rootdir*"v2head_"*"errors"*".tex", tabularline(head))
	to_file(rootdir*"v2line_"*"errors"*".tex", tabularline(res))
	#v3
	head	= ["\$# dummy frames\$";"\$# Bilder\$"; "\$\\e_{L^2}(I,s)\$"; "\$\\e_{reg_x}(w)\$"; "\$\\e_{L^2}(I,V)\$"]
	res		= [zwischen_ausgelassen; T; L2_errs[end]; H1_norm_w_noweight_space(u,v); L2norm(vorgabe_fehler)]
	to_file(rootdir*"v3table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )
	to_file(rootdir*"v3head_"*"errors"*".tex", tabularline(head))
	to_file(rootdir*"v3line_"*"errors"*".tex", tabularline(res))
end

function cmaptest(what, cmp)
	clf()
	surf(what[:,:,floor(3*T/4)]', cmap=cmp, cstride=1, rstride=1, linewidth=0.2)
	fig	= gcf()
	ax	= gca()
	ax[:set_zlim](0,1) # z-Achsen range auf 0,1 festlegen
	ax[:view_init](azim=-30, elev=55)
end

# pygui(true)
# cmaptest(I, "gray_r") #super scheisse
# cmaptest(I, "gray") #ok
# cmaptest(I, "GnBu")
# cmaptest(I, "GnBu_r")
# cmaptest(I, "Blues")
# cmaptest(I, "Blues_r")
# cmaptest(I, "coolwarm")
# cmaptest(I, "coolwarm_r")


function save_endergebnis(dir)
	try
		save("$(dir)res.jld", 
				"steps", steps,
				"I", I, 
				"u", u,
				"v", v, 
				"L2_errs", L2_errs,
				"H1_errs", H1_errs,
				"Js", Js,
				"H1_J_ws", H1_J_ws
				)
	catch e
		warn("ERGEBNIS KONNTE NICHT GESPEICHERT WERDEN!", e)
	end
end
