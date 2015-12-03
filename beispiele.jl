#img				= imread("test/test-0.png")
#s					= convert(Array, img)

function quadrat(y,x,t)
	cx	= floor(n/2);
	cy	= floor(m/2);

	d	= 2;
	cl	= cx - 8
	cr	= cx + 8

	#return (x >= cl +t*d).*(x <= cr +t*d).*(y >= cl).*(y <= cr);
	return (x >= cl ).*(x <= cr ).*(y >= cl+t*d).*(y <= cr+t*d);
	#return (x >= cl -t*d).*(x <= cr -t*d) .* (y >= cl-t*d).*(y <= cr-t*d);
end

function slotted_circle(y,x)
	slotd	= -9
	slotwh	=  5
	r		= 25
	return (sqrt(x^2+y^2) < r) * (1-( (abs(y)<slotwh)*(x>slotd) ))
end

function rot(t)
	return [cos(t) sin(t);-sin(t) cos(t)]
end

function rot_circle(y,x,t) 
	# thr vorsicht mit y, x vertauscht!
	rxy = rot(-1*t/pi/2)*[x-n/2;y-m/2]
	return slotted_circle(rxy[2], rxy[1])
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

function readtaxi()
	taxi = zeros(256+4,256+4,41)
	for i=1:41
		#@show a	= Images.imread("../taxi/taxi$(i).bmp")
		#a_	= convert(Array, a)
		#a_	= convert(Array, float64(a))
		#taxi[3:190+2,3:256+2,i]= a_[:,:]

		# von daheim vorbereitet, weil die Imagemagick-Bibliothek auf dem Uniserver zu alt ist fuer Julia
		a_	= readdlm("../taxi/taxi$(i).dlm")
		taxi[:,:,i] = a_[:,:]'
	end
	return taxi
end

