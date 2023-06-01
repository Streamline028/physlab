clc
clear all

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

L1=3.84e-3; 
M=256; 
dx1=L1/M; 
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;

lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(16/2);
z=180e-3;
zf=z;
z2=8e-3;
zf2=z2;
[X1,Y1]=meshgrid(x1,y1);
[theta,rho]=cart2pol(X1,Y1);
u1=zeros(M,M);
u1(rho<w)=1;
u1=exp(-(rho.^2)/w^2);
P1=angle(u1);
I1=(abs(u1).^2);
R1=radCal(I1);

figure(1)
imagesc(x1/1e-3,y1/1e-3,I1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R1*0.0150),'mm']); 
colorbar

figure(2) 
imagesc(x1/1e-3,y1/1e-3,P1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%% 1
[u1]=focus(u1,L1,lambda,zf);
u2=propTF(u1,L1,lambda,z);
P2=angle(u2);
I2=(abs(u2).^2);
R2=radCal(I2);

figure(3) 
imagesc(x1/1e-3,y1/1e-3,I2);
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R2*0.0150),'mm']); 
colorbar

figure(4) 
imagesc(x1/1e-3,y1/1e-3,P2); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%% 2
u2p=propTF(u2,L1,lambda,z2);
[u2p]=focus(u2p,L1,lambda,zf2);
u3=propTF(u2p,L1,lambda,z2);
P3=angle(u3); 
I3=(abs(u3).^2);
R3=radCal(I3);

figure(5) 
imagesc(x1/1e-3,y1/1e-3,I3);
xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R3*0.0150),'mm']); 
colorbar

figure(6) 
imagesc(x1/1e-3,y1/1e-3,P3); 
xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar