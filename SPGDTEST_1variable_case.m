clc
clear all % focused gaussian beam calculator

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

r=1000;
S=16;%16*4;
MM=16;
GG =0.5;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

figure(1); test_x=[0:0.001:4]; test_y=abs(exp(-1i*test_x*(2*pi))+exp(-1i*GG)).^2; plot(test_x, test_y);
hold on; for (ij = 1:size(test_x,2)-1) test_dy(ij)=(test_y(ij+1)-test_y(ij)); test_dy(ij)= test_dy(ij).*test_dy(ij); test_dx(ij) = test_x(ij); end; plot(test_dx, test_dy*10000+2, 'r'); hold off;

%%

bu = 0.2*(2*pi);
au = rand()*(2*pi);

u5 = exp(-1i*GG)+exp(-1i*bu);
ausave(1,1) = bu;
 
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
    
    u0 = exp(-1i*au)+exp(-1i*GG);
    Intensity = abs(u0).^2;
    ausave(1,ii+1) = au;
    
    J(1,ii+1) = Intensity;
    
%     J(1,ii+1) = J1(1,ii+1)/JJ(1,ii+1);
    W(1,ii+1)=((J(1,ii+1) - J(1,ii)));% /mean(J)
    weight=(W(1,ii+1)); %->multithread 1/Jratio(1,ii+1)
    BB=(au-bu);
    BB = BB/abs(BB) *2;
    dusave(1,ii) = BB;
%     diffUsave = (au-bu)./weight;
%     diffU = (diffUsave);
%     diffU(isnan(diffU)) = 0;
%     diffUsave(isnan(diffUsave(:,:,ii+1)),ii+1) = 0;
    ga(1,ii) = 2.6/J(1,ii+1);
    ga(isnan(ga)) = 1;
%     if sum(sum(BB))<=0.01
%         break
%     end
    WM = ga(1,ii)*(weight.*(BB)); %.*perturb.*rand([MM,MM])
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
    au = rem(au, (2*pi));
%     au = (au+(au.*weight))/2;%2
%     MU(1,ii+1)=max(max(au))/(2*pi);
%     au = au./ MU(1,ii+1);
    if (min(min(au))<0)
        au = au - min(min(au));
    end
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
subplot(2,2,1);
plot(J,'b*-'); 
title('beam 중심의 intensity'); 
subplot(2,2,2);
plot(W,'go-'); 
title('intensity의 변화량'); 
subplot(2,2,3);
plot(ausave./(2*pi),'r^-');
title('바꾸어 주는 빛의 위상'); 
subplot(2,2,4);
plot(dusave,'y^-'); 
title('위상의 변화량'); 

figure(10) 
plot(ga,'*-'); 
title('gamma'); 