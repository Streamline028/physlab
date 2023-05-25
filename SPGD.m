clc
clear all

L1=0.5*10; 
M=100; 
dx1=L1/M; 
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=0.5*10^-6;
k=2*pi/lambda; 
w=0.051*10; 
z=1000; 
r=1;
zf=500;
[X1,Y1]=meshgrid(x1,y1);
u1=rect(X1/(L1/2)).*rect(Y1/(L1/2)); 
I1=u1;%Initial_U(u1); %1-u1;
figure(1)
imagesc(x1,y1,I1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('z= 0 m');

%u2=propTF(u1,L1,lambda,z);
Before_U = zeros(M,M);
After_U = Initial_U(Before_U);
%u2=propTF(u1,L1,lambda,1); 
I2 = I1;

for i = 1:r
    I2=propMOD(I1,After_U,L1,lambda,z,zf); %Apply_U2(I1,After_U,lambda,z);
    %After_U apply to I2 >> get I2
    J(1,i+1) = chk_J(I2,2);
    u2=(After_U - Before_U).*((J(1,i+1) - J(1,i))/J(1,i+1)); %->multithread
    W(1,i)=(J(1,i+1) - J(1,i));
    Before_U = After_U;
    After_U = After_U + (u2);
    D(1,i) = sum(sum(u2));
end
x2=x1; 
y2=y1; 
figure(2) 
imagesc(x2,y2,I2); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title([num2str(r),' times']); 

x3=x1; 
y3=y1; 
figure(3) 
imagesc(x3,y3,After_U); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('align'); 

figure(4) 
hold off
plot(J,'b*-'); 
hold on
plot(W,'go-'); 
%plot(D,'r^-'); 
title('J'); 
hold off

figure(5) 
imagesc(x3,y3,u2); 
axis square; axis xy; 
colormap('gray'); xlabel('x (m)'); ylabel('y (m)'); 
title('weight'); 

% 
% figure(4) 
% plot(x2,abs(u2(M/2+1,:))); 
% xlabel('x (m)'); ylabel('Magnitude'); 
% title(['z= ',num2str(z),' m']); 
% 
% figure(5) 
% plot(x2,unwrap(angle(u2(M/2+1,:)))); 
% xlabel('x (m)'); ylabel('Phase (rad)'); 
% title(['z= ',num2str(z),' m']);