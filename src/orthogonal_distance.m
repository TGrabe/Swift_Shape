function [ distance ] = orthogonal_distance(Ax,Ay,Bx,By,theta)
r=(Ay-By+tan(theta)*(Bx-Ax))/(1+tan(theta)^2);
xl=-r*tan(theta);
yl=r;
distance=sqrt(xl^2+yl^2);
if theta==pi/2
    distance=abs(Bx-Ax);
end
if theta==0
    distance=abs(By-Ay);
end
end