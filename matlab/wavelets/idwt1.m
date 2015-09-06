function x = idwt1(a,d,wv,lx)
[Lo_R,Hi_R] = wfilters(wv,'r');
x = upsconv(a,Lo_R,lx) + upsconv(d,Hi_R,lx);