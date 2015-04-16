include("verfahren.jl")
tic()
l2 =0
L2 =0

for i=1:10000
	l2 = sample_err_l2(I,s,norm_s)[1]
end
toc()

tic()
for i=1:10000
	L2 = sample_err_L2(I,s,norm_s)[1]
end
toc()

tic()

echo(l2, L2)


