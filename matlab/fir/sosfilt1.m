function [y] = sosfilt1(sos, x)
l = length(x);
x = x(:)';
y = zeros(1,l);
ls = size(sos);

for j=1:ls(1)
   a0 = sos(j,4);
   a1 = sos(j,5)/a0;
   a2 = sos(j,6)/a0;
   b0 = sos(j,1)/a0;
   b1 = sos(j,2)/a0;
   b2 = sos(j,3)/a0;
   
   v1 = 0;
   v2 = 0;
   
   for i=1:l
      v0 = x(i) - a1*v1 - a2*v2;
      y(i) = b0*v0 + b1*v1 + b2*v2;
      v2 = v1;
      v1 = v0;
   end
   
   x = y;
end