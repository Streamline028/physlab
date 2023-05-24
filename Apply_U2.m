function[out]=Apply_U2(U,del_U,lambda, z); %U=del_u matrix J=delJ=J(u+del_u)-J(u)
[M,N]=size(U); %get input field array size 
[M1,N1]=size(del_U); %get input field array size
dx=0.5/M; %sample interval 
k=2*pi/lambda; %wavenumber 
fx=-1/(2*dx):1/0.5:1/(2*dx)-1/0.5; %freq coords
out=zeros(M,N);

%H=exp(-j*pi*lambda*(FX.^2+FY.^2)); %trans func 
    for i = 1:M1
        for m = 1:N1
            [FX,FY]=meshgrid(del_U(N1,:),del_U(:,M1).'); 
            H=exp(-j*pi*2*lambda*z*(FX.^2+FY.^2)); %trans func 
            H=fftshift(H); %shift trans func 
            U1=fft2(fftshift(U)); %shift, fft src field 
            U2=H.*U1; %multiply 
            ADD=ifftshift(ifft2(U2)).*exp(-j*k/(2*z)); %inv fft, center obs field
            out=out+ADD;
        end
    end

end