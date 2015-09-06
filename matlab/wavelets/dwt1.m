function [a,d] = dwt1(x,wv)

[Lo_D,Hi_D] = wfilters(wv,'d');
lf = length(Lo_D);
lx = length(x);

first = 2;
last = lx+lf-1;

% Half-point Symmetrization
I = [lf-1:-1:1 , 1:lx , lx:-1:lx-lf];
y = x(I);

% Approximation
z = conv1(y,Lo_D,'valid');
a = z(first:2:last);

% Details
z = conv1(y,Hi_D,'valid');
d = z(first:2:last);