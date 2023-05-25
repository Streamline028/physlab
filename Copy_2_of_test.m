clc
clear all % focused gaussian beam calculator

load('initset.mat');

r=1000;
S=16*4;%16*4;
MM=16;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

z=180e-3;
zf=z;
z2=8e-3;
zf2=z2;
[X1,Y1]=meshgrid(x1,y1);
[theta,rho]=cart2pol(X1,Y1);

u1=zeros(M,M);
u1=exp(-(rho.^2)/w^2);

perturb = (binornd(ones(16,16),ones(16,16)./2)-0.5); %.*(4*pi)

P1=angle(u1);
I1=(abs(u1).^2);
R1=radCal(I1);

Jratio_init=1/(sum(sum(I1))/chk_J(I1,S)-1)*100;

hold off
figure(1)
subplot(4,4,1);
imagesc(x1/1e-3,y1/1e-3,I1); 
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,[num2str(Jratio_init)],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R1*0.0150),'mm']); 
colorbar
hold off

% figure(2) 
subplot(4,4,2);
imagesc(x1/1e-3,y1/1e-3,P1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%%

bu = zeros(MM,MM).*(2*pi);
% bbu = bbu;
au = rand([MM,MM]).*(2*pi);
au = preset6;

ausave(:,:,1) = bu;

initu = au;

bu_ = padarray(expand(bu),[128 128],0,'both');

% figure(3) 
% imagesc(bu_); 
% axis square; axis xy; 
% colormap('gray');
% colorbar

uu=u1;

u5=propTF(u1.*exp(-1i*bu_),L1,lambda,2*z);
% [u5]=focus(u5,L1,lambda,zf);
u5=propTF(u5,L1,lambda,2*z);

P5=angle(u5)/pi+1; 
I5=(abs(u5).^2);
R5=radCal(I5);

JJ(1,1) = sum(sum(I5));
J(1,1) = chk_J(I5,S);
maxJ = 1;
% J(1,1) = 0;
Jratio(1,1)=(J(1,1)/JJ(1,1))*100;
% J(1,1) = J1(1,1)/JJ(1,1);
W(1,1)= 1;
D(1,1)= 1;
MU(1,1)=max(max(angle(bu)));
MaxI(1,1)=max(max((abs(I5).^2)));
MinI(1,1)=min(min((abs(I5).^2)));
diffUsave(:,:,1) = abs(au-bu);

hold off
% figure(3)
subplot(4,4,3);
imagesc(x1/1e-3,y1/1e-3,I5);
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,[num2str(Jratio(1,1))],'Color','red');
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R5*0.0150),'mm']); 
colorbar
hold off

% figure(4) 
subplot(4,4,4);
imagesc(x1/1e-3,y1/1e-3,P5); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar
%u0=uu.*expand(au);[u0]=focus(u0,L1,lambda,zf);u0=propTF(u0,L1,lambda,2*z);
for ii = 1:r
    
%     u0=uu;
    ausave(:,:,ii+1) = au;
    au_ = padarray(expand(au),[128 128],0,'both');

    u0=propTF(uu.*exp(-1i*au_),L1,lambda,2*z);
%     [u0]=focus(u0,L1,lambda,zf);
    u0=propTF(u0,L1,lambda,2*z);
    Intensity = abs(u0).^2;
    Isave(:,:,ii) = Intensity;
    
    JJ(1,ii+1) = sum(sum((Intensity)));
    J(1,ii+1) = chk_J((Intensity),S);
    
    Jratio(1,ii+1)=(J(1,ii+1)/JJ(1,ii+1))*100;
%     J(1,ii+1) = J1(1,ii+1)/JJ(1,ii+1);
    W(1,ii+1)=((J(1,ii+1) - J(1,ii)));% mean(J)
    weight=(W(1,ii+1)); %->multithread Jratio(1,ii+1)-Jratio(1,ii)
    BB=(au-bu);
%     diffUsave(:,:,ii+1) = (au-bu)./weight;
%     diffU = (diffUsave(:,:,ii+1));
%     diffU(isnan(diffU)) = 0;
%     if sum(sum(BB))<=0.01
%         break
%     end
    WM = (weight.*BB);
%     WM = (weight.*(BB+diffU)).*rand([MM,MM]); %.*perturb
%     diffU = (diffUsave(:,:,ii)-diffUsave(:,:,ii+1))/sign(W(1,ii)-W(1,ii+1));%
%     diffU = (diffUsave(:,:,ii+1))/sign(W(1,ii+1));
%     WM = (weight.*diffUsave(:,:,ii+1)); %.*rand([MM,MM])
%     if (ii == 1)
%         WM = rand([MM,MM]).*(2*pi);
%     end
    WM = rem((WM), (2*pi));   %.*(binornd(ones(16,16),ones(16,16)./2)-0.5).*(4*pi)
%     WM = WM./(max(max(WM))/(2*pi));
    D(1,ii+1) = sum(sum(WM))/sum(sum(au))+1;
    MaxI(1,ii+1)=max(max((Intensity)));
    MinI(1,ii+1)=min(min((Intensity)));
    bu = au;
%     bu = rem(bu, (2*pi));
%     aum = (au + weight);%2
    au = (au - WM);%2
%     au = rem(au, (2*pi));
%     au = (au+(au.*weight))/2;%2
    MU(1,ii+1)=max(max(au))/(2*pi);
    au = au./ MU(1,ii+1);
%     if (min(min(au))<0)
%         au = au - min(min(au));
%     end
%     au = exp(-i./(MU(1,ii+1))*(au));%2 
%     au = au - min(min(au));
%     if rem(ii,30) == 0
%         au=au+((rand([MM,MM])-0.5)*(2*pi));%| MU(1,ii+1) == 0 30
%     end
    if W(1,ii+1) == 0 %& half == 0
        break
    end
%     if J(1,i+1) == inf | D(1,i+1) == inf
%         break
%     end
    if ii==1
        hold off
        P3=angle(u0)/pi+1; 
        I3=(abs(u0).^2);
%         figure(11) 
        subplot(4,4,5);
        imagesc(x1/1e-3,y1/1e-3,I3);
        hold on
        rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
        text(S/2*dx1/1e-3, -S/2*dx1/1e-3,[num2str(Jratio(1,ii+1))],'Color','red');
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('result'); 
        colorbar
        hold off
        
        subplot(4,4,6);
        imagesc(x1/1e-3,y1/1e-3,P3); 
        % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('align'); 
        colorbar
    end
end

P3=angle(u0)/pi+1; 
I3=(abs(u0).^2);
R3=radCal(I3);

hold off
% figure(5) 
subplot(4,4,7);
imagesc(x1/1e-3,y1/1e-3,I3);
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,[num2str(Jratio(1,ii+1))],'Color','red');
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('result'); 
colorbar
hold off

% figure(6)
subplot(4,4,8);
imagesc(x1/1e-3,y1/1e-3,P3); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

% figure(7)  
subplot(4,4,9);
imagesc(initu); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(7)  
subplot(4,4,10);
imagesc(bu); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(8)                                                                                         % rescaled figure of J, delta J, change amount of au, cycle of au
subplot(4,4,11);
hold off
plot(rescale(J),'b*-');                                                                             % cacluated J value during iterration
hold on
plot(rescale(W),'go-');                                                                             % improvement of output during iterration
plot(rescale(D),'r^-');                                                                             % expectation of improvement of au
plot(rescale(MU),'y^-');                                                                            % au's cycle count
title('J'); 
legend('J','del J','u diff','u max');%
hold off

% figure(9)
subplot(4,4,12);
% subplot(2,1,1);
plot(MaxI,'b*-'); 
title('Max'); 
subplot(4,4,13);
% subplot(2,1,2);
plot(MinI,'ro-'); 
title('Min'); 

% figure(10) 
subplot(4,4,14);
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
returnVal = padarray(expand(round(bu*((2^16)/max(max(bu))))),[128 128],0,'both');
% returnVal = uint16(padarray(expand(round(bu)),[128 128],0,'both'));
returnVal = uint16(returnVal);
imagesc(returnVal); 
axis square; axis xy; 
colormap('gray');
colorbar
% imwrite(returnVal,"gray","ExportImage.tiff","BitDepth",16);
imwrite(returnVal,"ExportImage.tiff",'tiff');

% figure(12) 
subplot(4,4,15);
hold off
plot(J,'b*-'); 
hold on
title('J');  xlabel('iteration count'); ylabel('J value'); 
hold off

% figure(13) 
subplot(4,4,16);
hold off
plot(Jratio,'b*-'); 
hold on
title('PBR');  xlabel('iteration count'); ylabel('percentage'); 
hold off

figure(2) 
subplot(2,2,1);
plot(J,'b*-'); 
title('J'); 
subplot(2,2,2);
plot(W,'go-'); 
title('(J(1,ii+1) - J(1,ii))/J(1,ii+1)'); 
subplot(2,2,3);
plot(D,'r^-');
title('expectation of improvement of au'); 
subplot(2,2,4);
plot(MU,'y^-'); 
title('cycle count'); 


figure(100); subplot(2,1,2); plot(J,'b*-'); title('J');  xlabel('iteration count'); ylabel('J value');