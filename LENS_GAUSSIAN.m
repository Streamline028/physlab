clc
clear all

L1=50e-3; 
M=512; 
dx1=L1/M; 
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=0.5*10^-6;
k=2*pi/lambda; 
w=L1/4; 
z=100;
zf=100;
[X1,Y1]=meshgrid(x1,y1);
[theta,rho]=cart2pol(X1,Y1);
u1=zeros(M,M);
u1(rho<w)=1;
%u1=exp(-(rho.^2)/w^2);
P1=angle(u1);
I1=(abs(u1).^2);
R1=radCal_zero(I1);

figure(1)
imagesc(x1,y1,I1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title(['initial : ','radius = ',num2str(R1*dx1*1000),'mm']); 
colorbar

figure(2) 
imagesc(x1,y1,P1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('align'); 
colorbar

[u1]=focus(u1,L1,lambda,zf);
%[u2,L2]=propFF(u1,L1,lambda,z);
u2=propTF(u1,L1,lambda,z);
% dx2=L2/M;
P2=angle(u2);
x2=x1;%-L2/2:dx2:L2/2-dx2;
y2=y1;%x2;
I2=(abs(u2).^2);
R2=radCal_zero(I2);
figure(3) 
imagesc(x2,y2,I2);
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title(['result : ','radius = ',num2str(R2*dx1*1000),'mm']); 
colorbar

x3=x1; 
y3=y1; 
figure(4) 
imagesc(x3,y3,P2); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('align'); 
colorbar