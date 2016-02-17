@everywhere using HDF5
@everywhere using JLD

@everywhere begin
	@load "zwischenergebnis_5000.jld"
	m,n,T = size(I)
end

@everywhere I = I
#@everywhere s = s
#@everywhere n_samples			= size(s)
@everywhere n_samples			= 5
@everywhere n_zwischensamples	= div(T, n_samples-1)
@everywhere rootdir = "../out/taxi_alle5/$(m)_x_$(n)_$(n_samples)_$(n_zwischensamples)_$(alpha)_$(beta)_dx$(dx)dt$(dt)/"


include("view.jl")
include("beispiele.jl")

#s		= inits(quadrat)
#s		= inits(rot_circle)
#s		= readtaxi()[:,:, 1:10:end]



#@everywhere const n_samples			= 5
#@everywhere const n_zwischensamples	= 40    # duerfen nicht zu wenige sein? abhaengig von dt?
# ...................... T, alle ZeitPUNKTE, also T-1 Zeitschritte von einem Punkt auf den naechsten
#@everywhere const T					= (n_samples-1)*n_zwischensamples+1
# Zuordnung Samplenummer zu Zeitpunkt 
#samples_to_frames		= [ (k+1, k*n_zwischensamples+1) for k in 0:n_samples-1 ]


#save_all()
save_images_(I, "I")
#save_images_(s, "s")

_="fertig"


