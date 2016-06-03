git clone https://github.com/thraen/dipl

Requirements:

- julia, stable branch 0.4. installation:
	on the machine you want to run the program, julia has to be compiled locally:	
		git clone git://github.com/JuliaLang/julia.git
		git checkout release-0.4
		cd julia
		make
	
	make a link to ~/julia/julia:
		mkdir -p ~/bin; echo '~/julia/julia -p 6 $@' > ~/bin/julia
		(where -p 6 specifies the number of processes julia uses, only neccessary for parallel computations)
	open .bashrc and add line
		export PATH="$PATH:~/bin"

- julia packages:
	Pkg.add("PyPlot")
	Pkg.add("PyCall")
	Pkg.add("HDF5")
	Pkg.add("JLD")
	Pkg.add("Interpolations")
after installation it's often neccessary to restart julia.

- pyamg. installation:
	git clone https://github.com/pyamg/pyamg
	cd pyamg
	python setup.py install --user
