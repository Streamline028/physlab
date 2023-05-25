function[u2]=propMOD(u1,del_U,L,lambda,z,zf);
    [u1]=focus(u1,L,lambda,zf);
    [M1,N1]=size(del_U); %get input field array size
    [M,N]=size(u1); %get input field array size
    
    u2=zeros(M,N);
    dx=L/M; %sample interval
	k=2*pi/lambda; %wavenumber 
    x=-L/2:dx:L/2-dx; %spatial coords 
    [X,Y]=meshgrid(x,x); 
    for i = 1:M1
        for m = 1:N1
            [alpha,theta]=find_at(i-M/2,m-N/2,z)
            [T]=tilt(u1,L,lambda,alpha,theta);
            
            h=1/(j*lambda*z)*exp(j*k/(2*z)*(X.^2+Y.^2).*(del_U.^2)); %impulse 
            H=fft2(fftshift(h))*dx^2; %create trans func 
            U1=fft2(fftshift(T)); %shift, fft src field 
            U2=H.*U1; %multiply 
            ADD=ifftshift(ifft2(U2)); %inv fft, center obs field 
            u2=(u2+rescale(abs(ADD.^2).*2));
%             
%             L1=0.5*10; 
%             M=100; 
%             dx1=L1/M; 
%             x1=-L1/2:dx1:L1/2-dx1; 
%             y1=x1;
%             figure(9) 
%             imagesc(x1,y1,u2); 
%             axis square; axis xy; 
%             colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
%             title('align'); 
%             
%             del_U(i,m)^2
        end
    end
end