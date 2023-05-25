function[out]=jinc(x);
mask=(x~=0);
out=pi*ones(size(x));
out(mask)=besselj(1,2*pi*x(mask))./(x(mask));
end