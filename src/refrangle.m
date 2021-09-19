function [ delta2 ] = refrangle(delta1,alphat,n1,n2)
vz=1;
if mod(delta1-alphat,2*pi)>pi
   vz=-1; 
end

delta2=vz*acos(n1/n2*cos(vz*(delta1-alphat)))+alphat;
if isreal(delta2)
delta2=mod(delta2,2*pi);
end
end

