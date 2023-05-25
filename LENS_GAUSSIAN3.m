clc
clear all

L1=512*15e-6; %3.84e-3; 
M=512; %256
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

% L1=7.68e-3; 
% M=512; 
% dx1=L1/M; 
% x1=-L1/2:dx1:L1/2-dx1; 
% y1=x1;
% lambda=0.5*10^-6;
% k=2*pi/lambda; 
% w=dx1*(160/2); 

z=180e-3;
zf=z;
z2=8e-3;
zf2=z2;
[X1,Y1]=meshgrid(x1,y1);
[theta,rho]=cart2pol(X1,Y1);

RM= exp(-i*rand([16,16])*(2*pi));%

u1=zeros(M,M);
u1(rho<w)=1;
u1=exp(-(rho.^2)/w^2);

% u1=u1.*expand(RM);%/(2*pi)

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
% [u1]=focus(u1,L1,lambda,zf);
% u2=propTF(u1,L1,lambda,z);
% P2=angle(u2)/pi+1;
% I2=(abs(u2).^2);
% R2=radCal(I2);
% 
% figure(3) 
% imagesc(x1/1e-3,y1/1e-3,I2);
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title(['result : ','radius = ',num2str(R2*0.0150),'mm']); 
% colorbar
% 
% figure(4) 
% imagesc(x1/1e-3,y1/1e-3,P2); 
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title('align'); 
% colorbar

%% 2
% u3=propTF(u2,L1,lambda,z);
% P3=angle(u3)/pi+1; 
% I3=(abs(u3).^2);
% R3=radCal(I3);
% 
% figure(5) 
% imagesc(x1/1e-3,y1/1e-3,I3);
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title(['result : ','radius = ',num2str(R3*0.0150),'mm']); 
% colorbar
% 
% figure(6) 
% imagesc(x1/1e-3,y1/1e-3,P3); 
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title('align'); 
% colorbar

%% 3
% 
% u2p=propTF(u2,L1,lambda,z2);
% [u2p]=focus(u2p,L1,lambda,zf2);
% u4=propTF(u2p,L1,lambda,z2);
% P4=angle(u4)/pi+1; 
% I4=(abs(u4).^2);
% R4=radCal(I4);
% 
% figure(7) 
% imagesc(x1/1e-3,y1/1e-3,I4);
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title(['result : ','radius = ',num2str(R4*0.0150),'mm']); 
% colorbar
% 
% figure(8) 
% imagesc(x1/1e-3,y1/1e-3,P4); 
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title('align'); 
% colorbar

%% 2-modified

r=1000;
S=16*4;%16*4;

bu = exp(-i*ones(16,16).*(2*pi));
au = exp(-i*rand([16,16]).*(2*pi));

bu_ = padarray(expand(round(angle(bu))),[128 128],0,'both');

uu = u1;

% u5=;

u5=propTF(u1.*bu_,L1,lambda,2*z);
% [u5]=focus(u5,L1,lambda,zf);
u5=propTF(u5,L1,lambda,2*z);

P5=angle(u5)/pi+1; 
I5=(abs(u5).^2);
R5=radCal(I5);

JJ(1,1) = sum(sum(I5));
J(1,1) = chk_J(I5,S);
Jratio(1,1)=JJ(1,1)/J(1,1)-1;
% J(1,1) = J1(1,1)/JJ(1,1);
W(1,1)= 1;
D(1,1)= 1;
MU(1,1)=max(max(angle(bu)));
MaxI(1,1)=max(max((abs(I5).^2)));
MinI(1,1)=min(min((abs(I5).^2)));

figure(3) 
imagesc(x1/1e-3,y1/1e-3,I5);
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R5*0.0150),'mm']); 
colorbar

figure(4) 
imagesc(x1/1e-3,y1/1e-3,P5); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar
%u0=uu.*expand(au);[u0]=focus(u0,L1,lambda,zf);u0=propTF(u0,L1,lambda,2*z);
for ii = 1:1
    
%     u0=uu;
    au_ = padarray(expand(round(angle(au))),[128 128],0,'both');

    u0=propTF(uu.*au_,L1,lambda,2*z);
%     [u0]=focus(u0,L1,lambda,zf);
    u0=propTF(u0,L1,lambda,2*z);
    
    JJ(1,ii+1) = sum(sum((abs(u0).^2)));
    J(1,ii+1) = chk_J((abs(u0).^2),S);
    Jratio(1,ii+1)=JJ(1,ii+1)/J(1,ii+1)-1;
%     J(1,ii+1) = J1(1,ii+1)/JJ(1,ii+1);
    W(1,ii+1)=((J(1,ii+1) - J(1,ii))/J(1,ii+1));
    weight=(au - bu).*(W(1,ii+1)); %->multithread Jratio(1,i+i1)
    D(1,ii+1) = sum(sum(angle(weight)))/sum(sum(angle(au)))+1;
    MaxI(1,ii+1)=max(max((abs(u0).^2)));
    MinI(1,ii+1)=min(min((abs(u0).^2)));
    bu = au;
%     aum = (au + weight);%2
    aum = (angle(au) + angle(weight));%2
%     au = (au+(au.*weight))/2;%2
    MU(1,ii+1)=max(max(angle(au)))/(2*pi);
    au = exp(-i./(MU(1,ii+1))*(aum));%2 
%     au=au+(rand([16,16])-0.5);
    if W(1,ii+1) == 0 | MU(1,ii+1) == 0
        break
    end
%     if J(1,i+1) == inf | D(1,i+1) == inf
%         break
%     end
end

P3=angle(u0)/pi+1; 
I3=(abs(u0).^2);
R3=radCal(I3);

figure(5) 
imagesc(x1/1e-3,y1/1e-3,I3);
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R3*0.0150),'mm']); 
colorbar

figure(6)
imagesc(x1/1e-3,y1/1e-3,P3); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

figure(7)  
imagesc(angle(bu)); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

figure(8) 
hold off
plot(rescale(J),'b*-'); 
hold on
plot(rescale(W),'go-'); 
plot(rescale(D),'r^-'); 
plot(rescale(MU),'y^-'); 
title('J'); 
legend('J','del J','u diff','u max');%
hold off

figure(9)
subplot(2,1,1);
plot(MaxI,'b*-'); 
title('Max'); 
subplot(2,1,2);
plot(MinI,'ro-'); 
title('Min'); 

figure(10) 
% hold off
% plot((J),'b*-'); 
% hold on
% plot((W),'go-'); 
% % plot(rescale(D),'r^-'); 
% % plot(rescale(MU),'y^-'); 
% title('J'); 
% legend('J','del J','u diff','u max');%
% hold off

% returnVal = padarray(expand(bu),[128 128],0,'both');
returnVal = padarray(expand(round(angle(bu)*((2^16)/max(max(angle(bu)))))),[128 128],0,'both');
% returnVal = uint16(padarray(expand(round(bu)),[128 128],0,'both'));
returnVal = uint16(returnVal);
imagesc(returnVal); 
axis square; axis xy; 
colormap('gray');
colorbar
% imwrite(returnVal,"gray","ExportImage.tiff","BitDepth",16);
imwrite(returnVal,"ExportImage.tiff",'tiff');
