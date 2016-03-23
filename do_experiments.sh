#!/bin/bash

# for exp in 'exp_rot_disc.jl' 'exp_deform_disc.jl' 'exp_taxi_klein.jl';do
for exp in 'exp_taxi_klein.jl' ;do
	echo $exp
# 	for tr in 'false' 'true'; do
	for tr in 'false'; do
		echo time_reg $tr
		cat $exp | sed "s/time_regularization.*=.*/time_regularization = $tr/" > tmp.jl
		for i in '1.0' '0.5' '0.1' '0.01' '0.001' '0.0001' '0.00001' '0.000001';do
			echo $i
			cat tmp.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmptmp.jl

			~/Julia/julia -p2 tmptmp.jl >> mist.log
			sleep 6
		done
	done
done

# for i in '0' '1' '2' '3' '4' '5' '6';do
# for i in '1.0' '0.5' '0.000001';do
# 	echo $i
# 	cat  exp_n_zwischensampls_quadrat.jl| sed "s/const zwischen_ausgelassen.*=.*/const zwischen_ausgelassen = $i/" > tmp.jl
# 	~/Julia/julia -p3 tmp.jl
# 	~/julia/julia -p2 tmp.jl
# 	sleep 3
# done
