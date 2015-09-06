function y = upsconv(x,f,s)
% Dyadic upsampling.
p = 0;
rem2 = rem(p, 2);
addLEN = 2*rem2-1;
lux = 2*length(x) + addLEN;
ux = zeros(1,lux);
ux(1+rem2:2:end) = x;

y = conv1(ux,f,'full');

% Crop boundries
ly = length(y);
d = (ly-s)/2;
first = 1 + floor(d);
last = ly - ceil(d);
y = y(first:last);
