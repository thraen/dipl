using Grid

# coarsen image sequence by factor 2

dep	= 3

function constr_hierarchy(dx, dt, s, depth)
	ss	= cell(depth)
	ss[1]	= s
	for i=2:depth
		ss[i] = restrict(ss[i-1], [true, true, false])
		restrict(s, [true,true,false])
	end
	return ss
end

ss= constr_hierarchy(dx, dt, s, dep)

Is	= similar(ss)
us	= similar(ss)
vs	= similar(ss)
ps	= similar(ss)

_m, _n, _ = size(ss[dep])

function set_mn_on(worker, _m, _n)
	remotecall_wait(worker, (i,j)->(global m,n; @show m=i;@show n=j), _m, _n)
	#remotecall_wait(worker, (i,j)->(global m,n; @show m=i;@show n=j), _m, _n)
end

for worker = procs() set_mn_on(worker, _m, _n) end

echo("\nCoarsest Level", m, n, T, "\n")

u		= 0* ones( m, n, T-1 )
v		= 0* ones( m, n, T-1 )

maxsteps = 100

include("verfahren.jl")
@time Is[dep], us[dep], vs[dep], ps[dep], L2_err, H1_err, J, H1_J_w, _ = verfahren_grad(ss[dep], u, v, steps)


for i=length(ss)-1:-1:1
	echo("\nnext Level\n")
	global m,n
	@show i
	@show m, n
	@show _m, _n, _	= size(ss[i])


	for worker = procs() set_mn_on(worker, _m, _n) end

# 	@show m, n "warum? woher kommen diese variablen?"

	us[i] 	= prolong(us[i+1], [(_m,_n,T-1)...])/2
	vs[i] 	= prolong(vs[i+1], [(_m,_n,T-1)...])/2
	sleep(3)

	@show size(us[i])
	@show size(vs[i])

	global maxsteps
	maxsteps = i==1? 10000:100

	steps = 1
	include("verfahren.jl")
	@time Is[i], us[i], vs[i], ps[i], L2_err, H1_err, J, H1_J_w, _ = verfahren_grad(ss[i], us[i], vs[i], steps)
end

_="fertig"
