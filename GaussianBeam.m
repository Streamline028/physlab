function[out]=GaussianBeam(X,L); 
[M,N]=size(X); %get input field array size 
dx=L/M; %sample interval 
fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
[FX,FY]=meshgrid(fx,fx); 
H=exp(-j*pi*(FX.^2+FY.^2)); %trans func 
H=fftshift(H); %shift trans func 
U1=fft2(fftshift(X)); %shift, fft src field 
U2=H.*U1; %multiply 
out=ifftshift(ifft2(U2)); %inv fft, center obs field
out=rescale(abs(out.^2).*2);
end