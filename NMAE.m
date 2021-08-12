function nmae=NMAE(Original,Recover,Omega)
nmae=(norm(Recover(:)-Original(:),1))/(norm(Original(:),1)-norm(Original(Omega),1));
end