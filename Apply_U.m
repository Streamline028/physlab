function[out]=Apply_U(U,del_U,lambda, z); %U=del_u matrix J=delJ=J(u+del_u)-J(u)
[M,N]=size(U); %get input field array size 
dx=0.5/M; %sample interval 
k=2*pi/lambda; %wavenumber 
fx=-1/(2*dx):1/0.5:1/(2*dx)-1/0.5; %freq coords
[FX,FY]=meshgrid(fx,fx); 
%H=exp(-j*pi*lambda*(FX.^2+FY.^2)); %trans func 

H=exp(-j*pi*2*lambda*z*(del_U.^2)); %trans func 
H=fftshift(H); %shift trans func 
U1=fft2(fftshift(U)); %shift, fft src field 
U2=H.*U1; %multiply 
out=ifftshift(ifft2(U2)).*exp(-j*k/(2*z)); %inv fft, center obs field

% [M,N]=size(u1); %get input field array size 
% dx=L/M; %sample interval 
% k=2*pi/lambda; %wavenumber 
% fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
% [FX,FY]=meshgrid(fx,fx); 
% H=exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func 
% H=fftshift(H); %shift trans func 
% U1=fft2(fftshift(u1)); %shift, fft src field 
% U2=H.*U1; %multiply 
% u2=ifftshift(ifft2(U2)); %inv fft, center obs field
end