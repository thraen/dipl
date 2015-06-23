#Base.Multimedia.writemime(stream,::MIME"text/plain",x::Float64)=@printf("%1.2f",x)
function echo(args...)
	for x in args print( x, '\t') end 
	print('\n')
end


