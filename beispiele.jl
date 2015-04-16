#img				= imread("test/test-0.png")
#s					= convert(Array, img)

function quadrat(y,x,t)
	cx	= floor(n/2);
	cy	= floor(m/2);

	d	= 2;
	cl	= cx - 8
	cr	= cx + 8

	return (x >= cl +t*d).*(x <= cr +t*d).*(y >= cl).*(y <= cr);
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
	# vorsicht mit y, x vertauscht!
	rxy = rot(1*t/pi/2)*[x-n/2;y-m/2]
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
