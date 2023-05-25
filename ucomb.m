function[out]=ucomb(x);
x=round(x*10^6)/10^6;
out=rem(x,1)==0;
end 