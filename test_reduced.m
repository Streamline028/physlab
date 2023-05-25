clc
clear all % focused gaussian beam calculator

r=2000;
S=16;%16*4;
MM=16;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

%%

bu = 0.1*exp(-1i*0*(2*pi));
au = exp(-1i*rand()*(2*pi));

u5 = bu+bu;
 
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
    
    u0 = au+bu;
    Intensity = abs(u0).^2;
    
    J(1,ii+1) = Intensity;
    
%     J(1,ii+1) = J1(1,ii+1)/JJ(1,ii+1);
    W(1,ii+1)=((J(1,ii+1) - J(1,ii)));% /mean(J)
    weight=(W(1,ii+1)); %->multithread 1/Jratio(1,ii+1)
    BB=abs(au-bu);
    diffUsave = (au-bu)./weight;
    diffU = (diffUsave);
    diffU(isnan(diffU)) = 0;
%     diffUsave(isnan(diffUsave(:,:,ii+1)),ii+1) = 0;
    
%     if sum(sum(BB))<=0.01
%         break
%     end
    WM = (weight.*(BB+diffU)); %.*perturb.*rand([MM,MM])
%     diffU = (diffUsave(:,:,ii)-diffUsave(:,:,ii+1))/sign(W(1,ii)-W(1,ii+1));%
%     diffU = (diffUsave(:,:,ii+1))/sign(W(1,ii+1));
%     WM = weight.*(diffU);
%     WM = (weight.*diffUsave(:,:,ii+1)); %.*rand([MM,MM])
%     if (ii == 1)
%         WM = rand([MM,MM]).*(2*pi);
%     end
%     WM = rem((WM), (2*pi));   %.*(binornd(ones(16,16),ones(16,16)./2)-0.5).*(4*pi)
%     WM = WM./(max(max(WM))/(2*pi));
    D(1,ii+1) = sum(sum(WM))/sum(sum(au))+1;
    value(1,ii+1)=Intensity;
    bu = au;
%     bu = rem(bu, (2*pi));
%     aum = (au + weight);%2
    au = (au - WM);%2
%     au = rem(au, (2*pi));
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
title('J'); 
subplot(2,2,2);
plot(W,'go-'); 
% title('(J(1,ii+1) - J(1,ii))/J(1,ii+1)'); 
% subplot(2,2,3);
% plot(D,'r^-');
% title('expectation of improvement of au'); 
% subplot(2,2,4);
% plot(MU,'y^-'); 
% title('cycle count'); 