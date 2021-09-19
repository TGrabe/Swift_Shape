function [ sum ] = summe( a,b,c )
if ~exist('a','var')
a = 0;end
if ~exist('b','var')
b = 0;end
if ~exist('c','var')
c = 0;end

sum=a+b+c;

end