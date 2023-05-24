clear all
clc

L1=0.5; 
M=250; 
dx1=L1/M; 
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=0.5*10^-6;
k=2*pi/lambda; 
w=0.051; 
z=2000; 
deg=pi/180; 
alpha=5.0e-5; %rad 
theta=5*deg; 
zf=2000;
[X1,Y1]=meshgrid(x1,y1);
u1=rect(X1/(2*w)).*rect(Y1/(2*w));
u2=zeros(size(u1));
I2=zeros(size(u1));
I1=abs(u1.^2); 

figure(1)
imagesc(x1,y1,I1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('z= 0 m');
for i = 0:1:360
    alpha=i*deg;
    for j = 0:0.1:360
        theta=j*deg;
        [T]=tilt(u1,L1,lambda,alpha,theta); 
        %[T]=focus(T,L1,lambda,zf);
        u2=propIR(T,L1,lambda,z); %focus
        I2=I2+abs(u2.^2);
    end
end

x2=x1; 
y2=y1; 

figure(2) 
imagesc(x2,y2,I2); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title(['z= ',num2str(z),' m']); 

figure(3) 
plot(x2,I2(M/2+1,:)); 
xlabel('x (m)'); ylabel('Irradiance'); 
title(['z= ',num2str(z),' m']); 

figure(4) 
plot(x2,abs(u2(M/2+1,:))); 
xlabel('x (m)'); ylabel('Magnitude'); 
title(['z= ',num2str(z),' m']); 

figure(5) 
plot(x2,unwrap(angle(u2(M/2+1,:)))); 
xlabel('x (m)'); ylabel('Phase (rad)'); 
title(['z= ',num2str(z),' m']);