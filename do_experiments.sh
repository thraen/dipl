#!/bin/bash

for i in '1.0' '0.5' '0.1' '0.01' '0.001' '0.0001' '0.00001' '0.000001';do
# for i in '1.0' '0.5' '0.000001';do
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
# 	cat exp_rot_disc.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
	cat exp_taxi_klein.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
# 	Julia tmp.jl ####NEIN!
# 	julia_new tmp.jl
	~/julia/julia -p5 tmp.jl
	sleep 3
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
	echo $i
done

