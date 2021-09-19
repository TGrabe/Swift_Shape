function [ihit] = intersection(Kx,Ky,Sx,Sy,delta)
imax=length(Kx);
ileft=1;
iright=imax;
ihit=0;

if delta>0
    delta=mod(delta,2*pi);
end

angmin=atan2(Ky(1)-Sy,Kx(1)-Sx);
angmax=atan2(Ky(imax)-Sy,Kx(imax)-Sx);
angdiff=angmax-angmin;

for ni=1:ceil(log2(imax))
    
    i=floor((ileft+iright)/2);
    
    ang1=atan2(Ky(i)-Sy,Kx(i)-Sx);
    ang2=atan2(Ky(i+1)-Sy,Kx(i+1)-Sx);
    
    angbtw12=min(ang1-ang2+2*pi*((ang1-ang2)<0),ang2-ang1+2*pi*((ang2-ang1)<0));
    angbtw1d=min(ang1-delta+2*pi*((ang1-delta)<0),delta-ang1+2*pi*((delta-ang1)<0));
    angbtw2d=min(ang2-delta+2*pi*((ang2-delta)<0),delta-ang2+2*pi*((delta-ang2)<0));
    
    if angbtw12>angbtw1d && angbtw12>angbtw2d
        ihit=i;
        break;
    else
        ihit=0;
        if (ang1-delta)*angdiff<0
            ileft=i;
        else
            iright=i;
        end
    end
end
end