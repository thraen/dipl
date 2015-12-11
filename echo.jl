#Base.Multimedia.writemime(stream,::MIME"text/plain",x::Float64)=@printf("%1.2f",x)

function _echo(args...)
	for x in args print( x, '\t') end 
	print('\n')
end

function _echolog(args...)
	_echo(args...)
	try
		log = open(rootdir * "log", "a")
		for x in args write(log,  string(x), '\t') end 
		write(log, ('\n'))
		close(log)
	catch e
		warn(e)
	end
end

echo=_echo
