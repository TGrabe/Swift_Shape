function [delta_alpha] = angbtwn(alpha1,alpha2)
delta_alpha=min(alpha1-alpha2+2*pi*((alpha1-alpha2)<0),alpha2-alpha1+2*pi*((alpha2-alpha1)<0));
end