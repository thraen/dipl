using PyPlot
@everywhere using PyPlot
@everywhere pygui(false)

@everywhere vsuff=".dlm"

function delete_output()
	run(`rm $rootdir -r`)
end

# pygui(true)
# cmap = "gray_r"
# cmap = "gray"
# cmap = "GnBu"
# cmap = "GnBu_r"
# cmap = "Blues"
# cmap = "Blues_r"
# cmap = "coolwarm"
# cmap = "coolwarm_r"

@everywhere function save_image_t(im, t, dir, pref, isuff=".png", dpi=100, cmap="gray", vmin=0, vmax=1, axlbls=false)
	imshow(im[:,:,t], cmap=cmap, interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
	if ~axlbls
		fig	= gcf()
		ax	= gca()
		ax[:set_xticklabels]([])  #keine labels
		ax[:set_yticklabels]([])  #keine labels
	end

	savefig(dir * pref * lpad(t, 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

function save_images(im, ts, dir, pref, isuff, dpi, cmap="gray", vmin=0, vmax=1, axlbls=false)
	println("save_images", dir, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $dir`)
	@sync @parallel for t in ts
		save_image_t(im, t, dir, pref, isuff, dpi, cmap, vmin, vmax, axlbls)
	end
	toc()
end

@everywhere function save_surf_t(im, t, dir, pref, isuff=".png", dpi=100, cmap="gray", zlim=(0,1), azim=-39, elev=53, axlbls=false)
	##
	## surf(X,Y,Z)!! man muss die arrays also transponieren!
	##
	surf(im[:,:,t]', cmap=cmap, cstride=1, rstride=1, linewidth=0.2)
	ax	= gca()
	ax[:set_zlim](zlim) # z-Achsen range auf 0,1 festlegen
	ax[:view_init](azim=azim, elev=elev)
	if ~axlbls
		ax[:set_xticklabels]([])  #keine labels
		ax[:set_yticklabels]([])  #keine labels
	end

	savefig(dir * pref * lpad(t, 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

function save_surfs(im, ts, dir, pref, isuff, dpi, cmap="gray", zlim=(0,1), azim=-39, elev=53, axlbls=true)
	println("save_surfs", dir, pref)
	tic()
	m, n, T = size(im)
	run(`mkdir -p $dir`)
	@sync @parallel for t in ts
		save_surf_t(im, t, dir, pref, isuff, dpi, cmap, zlim, azim, elev, axlbls)
	end
	toc()
end


@everywhere function save_quiver_t(u, v, t, dir, pref, isuff=".png", dpi=100, alle=1, axlbls=true)
	mu, nu, T = size(u)
	mv, nv, T = size(v)
	mpad = abs(mu-mv)
	npad = abs(nu-nv)

# 	quiver( [zeros(m, npad) u[:,:,t]], [zeros(mpad, n); v[:,:,t]] )
# 	quiver( [zeros(m, npad) u[:,:,t]], [zeros(mpad, n); v[:,:,t]], width=0.0015 )
	quiver( [zeros(m, npad) u[:,:,t]][1:alle:end, 1:alle:end], [zeros(mpad, n); v[:,:,t]][1:alle:end, 1:alle:end], width=0.0015 )
	if ~axlbls
		ax	= gca()
		ax[:set_xticklabels]([])  #keine labels
		ax[:set_yticklabels]([])  #keine labels
	end

	savefig(dir * pref * lpad(t, 4,"0") * isuff, dpi=dpi, bbox_inches="tight", pad_inches=0)
	clf()
end

function save_quivers(u, v, ts, dir, pref, isuff, dpi, alle=1, axlbls=true)
	println("save_quivers", dir, pref)
	println(size(u), size(v))
	tic()
	run(`mkdir -p $dir`)

	@sync @parallel for t in ts
		save_quiver_t(u, v, t, dir, pref, isuff, dpi, alle, axlbls)
	end
	toc()
end

function save_all()
	save_images(s, 1:size(s,3), rootdir*"s/", "img", ".png", 60)
	save_images(I, 1:size(I,3), rootdir*"I/", "img", ".png", 60)
	save_surfs(I,1:size(I,3), rootdir*"I/", "srf", ".png", 60)
# 	save_images(I_vorgabe, 1:size(I,3), rootdir*"I_vorgabe/", "img", ".png", 60)
# 	save_surfs(I_vorgabe,1:size(I,3), rootdir*"I_vorgabe/", "img", ".png", 60)
# 	save_images(diff_vorgabe, 1:size(I,3), rootdir*"verr/", "img", ".png", 60)
# 	save_surfs(diff_vorgabe,1:size(I,3), rootdir*"verr/", "srf", ".png", 60)
	save_quivers(u, v, 1:size(u,3), rootdir*"w/", "", ".png", 60, 2)
# 	save_quivers(grd_u_J, grd_v_J,"grad_J", 1:size(u,3), rootdir*"grad/" ,"", ".png", 60)
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
		#fn	= pref * lpad(t, 4,"0") * vsuff
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

function save_demo_rot_disc(outputs)
	for op in outputs
		suff, dpi= op[1], op[2]
# 		save_displacement(rootdir, suff, dpi)
		vorgabe_frames	= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
		auswahl = collect(vorgabe_frames)
# 		save_images(I, auswahl, rootdir, "I", suff, dpi, "gray_r", 0, 1, false)
# 		save_surfs(I, auswahl, rootdir, "I_srf", suff, dpi, "gray", (0,1), -39, 53, true)
# 		save_images(abs(vorgabe_fehler), 1:size(vorgabe_fehler,3), rootdir, "verr", suff, dpi, "gray_r", 0, 1, false)
		#das Geschwindigkeitsfeld enthaelt einen Zeitpunkt weniger als die Bilder
		save_quivers(u, v, auswahl[1:end-1], rootdir, "w", suff, dpi, 2, false) # nur alle 2 pfeile plotten
	end
end

function save_demo_taxi(outputs)
	for op in outputs
		suff, dpi= op[1], op[2]
# 		save_displacement(rootdir, suff, dpi)
		vorgabe_frames	= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
		auswahl = collect(vorgabe_frames)
# 		save_images(I, auswahl, rootdir, "I", suff, dpi, "gray", 0, 1, false)
# 		save_surfs(I, auswahl, rootdir, "I_srf", suff, dpi, "gray", (0,1), -39, 53, true)
# 		save_images(abs(vorgabe_fehler), 1:size(vorgabe_fehler,3), rootdir, "verr", suff, dpi, "gray_r", 0, 1, false)
		save_quivers(u, v, auswahl[1:end-1], rootdir, "w", suff, dpi, 2, false) # nur alle 2 pfeile plotten
	end
end

function demo_table(capt, label)
# 	head	= ["\$\\alpha\$", "\$\\beta\$", "\$\\sum\\|I-s\\|_2^2\$", "Reg.-Fehler", "\$\\sum \\|V-I\\|_2^2\$"]
	head	= ["\$\\alpha\$"; "\$\\beta\$"; "\$\\e_{L^2}(I,s)\$"; "\$\\e_{reg_x}(w)\$"; "\$\\e_{reg_t}(w)\$"; "\$\\e_{L^2}(I,V)\$"]
	bet		= (time_regularization == false) ? 0 : beta
	res		= [alpha; bet; L2_errs[end]; H1_norm_w_noweight_space(u,v); H1_norm_w_noweight_time(u,v); L2norm(vorgabe_fehler)]
	to_file(rootdir*"table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )

	to_file(rootdir*"head_"*"errors"*".tex", tabularline(head))
	to_file(rootdir*"line_"*"errors"*".tex", tabularline(res))
# 	print( rootdir*"table_"*"errors"*".tex", latextable_normal(capt, label, head, res) )
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

