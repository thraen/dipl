function char_quadrat2d(y, x, cy, cx, ah)
	return (x >= cx-ah ).*(x <= cx+ah ).*(y >= cy-ah).*(y <= cy+ah);
end

function char_quadrat(m,n,T)
	ah	= floor( min(n,m)/ 6 )

	cx	= (n+1)/2
	cy	= (m+1)/2

	return (y,x,t) -> char_quadrat2d(y,x, cy+t-floor(T/2), cx, ah)
# 	return (y,x,t) -> char_quadrat2d(y,x, cy+t, cx, ah)
end

function init_vorgabe(f_generator,m,n,T)
	V	= zeros(m, n, T)
	f	= f_generator(m,n,T)

	for t = 0:T-1
		for j = 1:n
			for i = 1:m
				V[i,j,t+1] = f(i,j,t)
			end
		end
	end
	return V
end

function quadrat(y,x,t)
	cx	= floor(n/2)
	cy	= floor(m/2)

	d	= 2
	cl	= cx - 8
	cr	= cx + 8

	#return (x >= cl +t*d).*(x <= cr +t*d).*(y >= cl).*(y <= cr);
	return (x >= cl ).*(x <= cr ).*(y >= cl+t*d).*(y <= cr+t*d);
	#return (x >= cl -t*d).*(x <= cr -t*d) .* (y >= cl-t*d).*(y <= cr-t*d);
end

function __slotted_circle(y,x,cy=0,cx=0)
	slotd	= -9
	slotwh	=  5
	r		= 25
	y	-= cy
	x	-= cx
	return (sqrt(x^2+y^2) < r) * (1-( (abs(y)<slotwh)*(x>slotd) ))
end

function __rot_circle_ex(m,n,T) 
	# thr vorsicht mit y, x vertauscht!
	w=0.25 
	dings = function(y,x,t)
		rxy = rot(-1*w*t/pi/2 - pi/8)*[x-n/2;y-m/2]
		return __slotted_circle(rxy[2], rxy[1], 40, 0)
	end
	return dings
end

function slotted_circle(y,x,m,n,T,cy=0,cx=0)
# 	slotd	= -0.05*m #ich setze ein quadratisches Omega voraus
# 	slotwh	=  0.02*m
# 	r		= (1/8)*m 

	slotd	= -0.05*m #ich setze ein quadratisches Omega voraus
	slotwh	=  0.02*m
	r		= (1/10)*m 

	y	-= cy
	x	-= cx
	return (sqrt(x^2+y^2) < r) * (1-( (abs(y)<slotwh)*(x>slotd) ))
end

function rot(t)
	return [cos(t) sin(t);-sin(t) cos(t)]
end

function rot_circle_ex(m,n,T) 
	w=0.1
	fun = function(y,x,t)
		rxy = rot(-1*w*t/pi/2 - pi/4 +(T/2)*w/pi/2)*[x-n/2;y-m/2]
# 		rxy = rot(-1*w*t/pi/2 - pi/4 )*[x-n/2;y-m/2]
		return slotted_circle(rxy[2], rxy[1], m,n,T, m/3, 0)
	end
	return fun
end

function rot_circle(y,x,t) 
	# thr vorsicht mit y, x vertauscht!
	rxy = rot(-1*t/pi/2)*[x-n/2;y-m/2]
	return slotted_circle(rxy[2], rxy[1])
end

function _rot_circle_ex(y,x,t) 
	# thr vorsicht mit y, x vertauscht!
	#w=0.25 
	#w=0.5 
	w=1.0 
	#w=2.0
	#w=3.0 
	rxy = rot(-1*w*t/pi/2)*[x-n/2;y-m/2]
	return slotted_circle(rxy[2], rxy[1], 40, 0)
end

function inits(f)
	s	= zeros(m, n, n_samples)
	for t = 0:n_samples-1
		for j = 1:n
			for i = 1:m
				s[i,j,t+1] = f(i,j,t)
			end
		end
	end
	return s
end

# thr entfernen. so wurden die alten Taxidlms erstellt. 
function taxi_bmp2dlm_images_alt!() 
	taxi = zeros(256+4,256+4,41)
	for i=1:41
		#raus kommt ein fixedpoint grauwert bild
		a	= Images.imread("../taxi/taxi$(i).bmp")
		
		a_	= convert(Array, float64(a))
		a_	= Array{Float64}(a_)
		
		taxi[3:190+2,3:256+2,i]= a_[:,:]
		writedlm("../taxi/taxi$(i).dlm", taxi[:,:,i])
	end
end

# thr entfernen. so wurden die alten Taxidlms eingelesen
function readtaxi_alt()
	taxi = zeros(256+4,256+4,41)
	for i=1:41
		#@show a	= Images.imread("../taxi/taxi$(i).bmp")
		#a_	= convert(Array, a)
		#a_	= convert(Array, float64(a))
		#taxi[3:190+2,3:256+2,i]= a_[:,:]

		# von daheim vorbereitet, weil die Imagemagick-Bibliothek auf dem Uniserver zu alt ist fuer Julia
		a_	= readdlm("../taxialt/taxi$(i).dlm")
		#taxi[:,:,i] = a_[:,:]' #thr deshalb waren die Bilder gedreht. weiss nicht, warum ich das gemacht habe.
		taxi[:,:,i] = a_[:,:]
	end
	return taxi
end

# thr entfernen: die alte Methode zur Erstellung der DLMs ohne auffuellen auf 260x260
function taxi_bmp2dlm_images!()
	for i=1:41
		#raus kommt ein fixedpoint grauwert bild
		a	= Images.imread("../taxi/taxi$(i).bmp")
		a_	= convert(Array, float64(a))
		a_	= Array{Float64}(a_)
		writedlm("../taxi/taxi$(i).dlm", a_)
	end
end

function taxi_bmp2dlm_pyplot!()
	for i=1:41
		#raus kommt ein Uint8 rgb-bild normalisiert auf 0-255,
		#wobei r, g, b-Komponenten gleich sind
		a	= PyPlot.imread("../taxi/taxi$(i).bmp")

		# Da die Daten aus Bildern stammen, ist der Ursprung links OBEN
		# wir drehen die y-Achse um, damit die Image-Plots mit surf-Plots
		# kompatibel sind.
		a	= a[end:-1:1,:,1] # wir nehmen nur die erste Komponente und invertieren die y-Achse 
		a_	= Array{Float64}(a)/255
		writedlm("../taxi/taxi$(i).dlm", a_)
	end
end

# erwartet auf [0,1] normierte und von 1-T numerierte dlm-Dateien
# die Graustufenbilder enthalten. 
function image_seq_from_dlm(dir, pref, T, padu=0, padd=0, padl=0, padr=0)
	@show (dir, pref, T, padl, padr, padu, padd)
	m,n = size( readdlm(dir*pref*"1.dlm") )
	res = zeros(m+padu+padd, n+padl+padr,T)
	for i=1:T
		a_	= readdlm(dir*pref*"$(i).dlm")
		#res[:,:,i] = a_

		res[padd+1 : padd+m , padl+1 : padl+n , i]= a_[:,:]
	end
	return res
end


function load_taxi(m,n,T)
	taxi = image_seq_from_dlm("../taxi/", "taxi", T)
	@show tm, tn, tT = size(taxi)
	taxi = image_seq_from_dlm("../taxi/", "taxi", T, 2, (m-2)-tm, 2, 2)
	return taxi
end

function diff_vorgabe(I_vorgabe, I, auslassen, zwischen_ausgelassen)
	m,n,T_vorgabe	= size(I_vorgabe)
	vorgabe_frames	= (1:(zwischen_ausgelassen+1):(zwischen_ausgelassen+1)*T_vorgabe) 
	# 	vorgabe_to_frames		= [(k,vorgabe_frames[k]) for k in 1:T_vorgabe] #wird nicht wirklich gebraucht
	vorgabe_fehler	= zeros( size(I_vorgabe) )
	for t in 1:T_vorgabe
		j						= vorgabe_frames[t]
		vorgabe_fehler[:,:,t]	= I_vorgabe[:,:,t] - I[:,:,j]
	end
	return  vorgabe_fehler
end

@everywhere function auswahl_vorgabe(auslassen, n_samples)
	used_indices	= (1:(auslassen+1):(auslassen+1)*n_samples) 
	# 	samples_to_vorgabe	= [(k, used_indices[k]) for k in 1:n_samples]
	return used_indices
end

#taxi_bmp2dlm_pyplot!()

#taxipad = load_taxi(260, 260, 5)
#imshow(taxipad[:,:,1], origin="lower")
#@show size(taxipad)
#_ ="fertig"

#PyPlot.imshow(taxi_neu, origin="lower")
#PyPlot.imshow(taxi_alt, origin="lower")
#maximum(abs(taxi_neu-taxi_alt))

