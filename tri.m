function[out]=tri(x)
t=1-abs(x);
mask=abs(x)<=1;
out=t.*mask;
end