clc
clear all % focused gaussian beam calculator
close all;

r=100;
S=16;%16*4;
MM=16;
GG =0.5;
filename = 'test_reduced2.gif';

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);
test_x=[0:0.001:2*pi]; test_y=[0:0.001:2*pi];
[X,Y] = meshgrid(test_x,test_y); test_z=abs(exp(-1i*X)+exp(-1i*Y)).^2;
% figure(1); 
% imagesc(test_x,test_y,test_z);
% axis image; axis xy;
% colormap('gray')
% xlabel('x');ylabel('y');
% hold on; for (ij = 1:size(test_x,2)-1) test_dy(ij)=(test_y(ij+1)-test_y(ij)); test_dy(ij)= test_dy(ij).*test_dy(ij); test_dx(ij) = test_x(ij); end; plot(test_dx, test_dy*10000+2, 'r'); hold off;

%%

bu = rand(1,2)*(2*pi);
au = rand(1,2)*(2*pi);

u5 = exp(-1i*bu(1,1))+exp(-1i*bu(1,2));
ausave(:,1) = bu;
 
I5=(abs(u5).^2);

J(1,1) = I5;
% J(1,1) = J1(1,1)/JJ(1,1);
W(1,1)= 1;
D(1,1)= 1;
MU(1,1)=max(max(angle(bu)));
MaxI(1,1)=max(max((abs(I5).^2)));
MinI(1,1)=min(min((abs(I5).^2)));
diffUsave(:,:,1) = abs(au-bu);

for ii = 1:r
    
    u0 = exp(-1i*au(1,1))+exp(-1i*au(1,2));
    Intensity = abs(u0).^2;
    ausave(:,ii+1) = au;
    
    J(1,ii+1) = Intensity;
    
%     J(1,ii+1) = J1(1,ii+1)/JJ(1,ii+1);
    W(1,ii+1)=((J(1,ii+1) - J(1,ii)));% /mean(J)
    weight=(W(1,ii+1)); %->multithread 1/Jratio(1,ii+1)
    BB=(au-bu);
    BB = BB/abs(sum(sum(BB)));
    dusave(:,ii) = BB;
%     diffUsave = (au-bu)./weight;
%     diffU = (diffUsave);
%     diffU(isnan(diffU)) = 0;
%     diffUsave(isnan(diffUsave(:,:,ii+1)),ii+1) = 0;
    gamma = (pi)./(J(1,ii+1)+1);
%     if sum(sum(BB))<=0.01
%         break
%     end
    WM = (weight.*(BB)); %.*perturb.*rand([MM,MM])gamma .* 
%     diffU = (diffUsave(:,:,ii)-diffUsave(:,:,ii+1))/sign(W(1,ii)-W(1,ii+1));%
%     diffU = (diffUsave(:,:,ii+1))/sign(W(1,ii+1));
%     WM = weight.*(diffU);
%     WM = (weight.*diffUsave(:,:,ii+1)); %.*rand([MM,MM])
%     if (ii == 1)
%         WM = rand([MM,MM]).*(2*pi);
%     end
    WM = rem((WM), (2*pi));   %.*(binornd(ones(16,16),ones(16,16)./2)-0.5).*(4*pi)
%     WM = WM./(max(max(WM))/(2*pi));
    D(1,ii+1) = sum(sum(WM))/sum(sum(au))+1;
    value(1,ii+1)=Intensity;
    bu = au;
%     bu = rem(bu, (2*pi));
%     aum = (au + weight);%2
    au = (au - WM);%2
    au = mod(au, (2*pi));
%     au = (au+(au.*weight))/2;%2
%     MU(1,ii+1)=max(max(au))/(2*pi);
%     au = au./ MU(1,ii+1);
    if (min(min(au))<0)
        au = au - min(min(au));
    end
    figure(100);
    set(gcf,'position',[456,411,560,420]);
    clf;
    imagesc(test_x,test_y,test_z);
    line(0.1*cos(0:0.01:2*pi)+au(1),0.1*sin(0:0.01:2*pi)+au(2),'Color','r','Linewidth',2);
    axis image; axis xy;
    colormap('gray')
    title(sprintf("ii = %d, J = %.2f\nx = %.2f, y = %.2f, d\\phi = %.3f",ii,J(ii),au(1)/2/pi,au(2)/2/pi,(au(2)-au(1))/2/pi))
    xlabel('x');ylabel('y');
    xlim([0 2*pi]); ylim([0 2*pi]);
    drawnow
    frame = getframe(100); 
    img = frame2im(frame);
    [imind cm] = rgb2ind(img,256); 
    if ii == 1
        imwrite(imind,cm,filename,'gif','Loopcount',1,'DelayTime',1/5); 
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/5); 
    end
% %     figure(101);
% %     set(gcf,'position',[1024,410,560,420]);
% %     clf
% %     subplot(2,2,1);
% %     plot(J,'b*-'); 
% %     title('J'); 
% %     subplot(2,2,2);
% %     plot(W,'go-'); 
% %     title('dJ'); 
% %     subplot(2,2,3);
% %     hold on;
% %     plot((ausave(1,:)./(2*pi)),'r');
% %     plot((ausave(2,:)./(2*pi))+2,'b');
% %     hold off;
% %     title('u'); 
% %     subplot(2,2,4);
% %     hold on;
% %     plot((dusave(1,:)./(2*pi)),'r');
% %     plot((dusave(2,:)./(2*pi)),'b');
% %     hold off;
% %     title('du');
% %     drawnow
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
    
end

I3=(abs(u0).^2);

% % figure(8)                                                                                         % rescaled figure of J, delta J, change amount of au, cycle of au
% subplot(4,4,11);
% hold off
% plot(rescale(J),'b*-');                                                                             % cacluated J value during iterration
% hold on
% plot(rescale(W),'go-');                                                                             % improvement of output during iterration
% plot(rescale(D),'r^-');                                                                             % expectation of improvement of au
% plot(rescale(MU),'y^-');                                                                            % au's cycle count
% title('J'); 
% legend('J','del J','u diff','u max');%
% hold off

% figure(10) 
% subplot(4,4,14);
% hold off
% plot((J),'b*-'); 
% hold on
% plot((W),'go-'); 
% % plot(rescale(D),'r^-'); 
% % plot(rescale(MU),'y^-'); 
% title('J'); 
% legend('J','del J','u diff','u max');%
% hold off

figure(2)
set(gcf,'position',[1024,410,560,420]);
clf
subplot(2,2,1);
plot(J,'b*-'); 
title('beam 중심의 intensity'); 
subplot(2,2,2);
plot(W,'go-'); 
title('intensity의 변화량'); 
subplot(2,2,3);
hold on;
plot((ausave(1,:)./(2*pi)),'r');
plot((ausave(2,:)./(2*pi)),'b');
hold off;
title('각 beam의 위상'); 
subplot(2,2,4);
hold on;
plot((dusave(1,:)./(2*pi)),'r');
plot((dusave(2,:)./(2*pi)),'b');
hold off;
title('각 beam의 위상의 변화량'); 