function[uout]=focus(uin,L,lambda,zf)
[M,N]=size(uin); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber
x=-L/2:dx:L/2-dx; %coords
[X,Y]=meshgrid(x,x);
uout=uin.*exp(-j*k/(2*zf)*(X.^2+Y.^2)); %apply focus
end