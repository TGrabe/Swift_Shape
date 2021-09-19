function [ distance ] = projected_distance( x1,y1,x2,y2,ang)

r=(y1-y2+tan(ang)*(x2-x1))/(1+tan(ang)^2);
xl=-r*tan(ang);
yl=r;
distance=sqrt(xl^2+yl^2);

if ang==pi/2
    distance=abs(x2-x1);
end
if ang==0
    distance=abs(y2-y1);
end
end