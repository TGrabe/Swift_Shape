function [esurv] = energysurvived(delta,alphat,n1,n2,etype)
alpha=pi/2-acos(cos(delta)*cos(alphat)+sin(delta)*sin(alphat));
beta=asin(n1/n2*sin(alpha));

rs=(n1*cos(alpha)-n2*cos(beta))/(n1*cos(alpha)+n2*cos(beta));
rp=(n2*cos(alpha)-n1*cos(beta))/(n2*cos(alpha)+n1*cos(beta));
R=(rs^2+rp^2)/2;

if ~isreal(R)
    R=1;
end

if etype==1
    esurv=1-R;
else
    esurv=R;
end
end