#!/bin/bash

for i in '1.0' '0.5' '0.1' '0.01' '0.001' '0.0001' '0.00001' '0.000001';do
	cat exp_rot_disc.jl | sed "s/const alpha.*/const alpha = $i/" | sed "s/const beta.*/const beta = $i/" > tmp.jl
	Julia tmp.jl
done
