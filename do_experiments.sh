#!/bin/bash

# for i in '1.0' '0.5' '0.1' '0.01' '0.001' '0.0001' '0.00001' '0.000001';do
for i in '1.0' '0.5' '0.1' '0.001' '0.0001' '0.00001' '0.000001';do
	echo $i
# 	cat exp_rot_disc.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
# 	cat exp_taxi_klein.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
	cat exp_deform_disc.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
# 	~/julia/julia -p2 tmp.jl
	~/Julia/julia -p2 tmp.jl
	sleep 3
done

# for i in '0' '1' '2' '3' '4' '5' '6';do
# for i in '1.0' '0.5' '0.000001';do
# 	echo $i
# 	cat  exp_n_zwischensampls_quadrat.jl| sed "s/const zwischen_ausgelassen.*=.*/const zwischen_ausgelassen = $i/" > tmp.jl
# 	~/Julia/julia -p3 tmp.jl
# 	~/julia/julia -p2 tmp.jl
# 	sleep 3
# done
