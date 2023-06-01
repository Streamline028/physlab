clc
clear all % focused gaussian beam calculator
close all;

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

Iterration_Count=100;

M=512;
L1=M*15e-6; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

%%

Before_U = rand(1,3)*(2*pi);
After_U = rand(1,3)*(2*pi);

Before_Beam_Flow = exp(-1i*Before_U(1,1))+exp(-1i*Before_U(1,2))+exp(-1i*Before_U(1,3));
Usave(:,1) = Before_U;
 
Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

Target_Intensity_Sum(1,1) = Before_Beam_Intensity;
dJ(1,1)= 1;
dUsave(:,:,1) = abs(After_U-Before_U);

for ii = 1:Iterration_Count
    
    After_Beam_Flow = exp(-1i*After_U(1,1))+exp(-1i*After_U(1,2));
    After_Beam_Intensity = abs(After_Beam_Flow).^2;
    Usave(:,ii+1) = After_U;
    
    Target_Intensity_Sum(1,ii+1) = After_Beam_Intensity;
    
    dJ(1,ii+1)=((Target_Intensity_Sum(1,ii+1) - Target_Intensity_Sum(1,ii)));
    Gamma = (pi)./(Target_Intensity_Sum(1,ii+1));
    Gamma(isnan(Gamma)) = 1;
    weight=Gamma .* (dJ(1,ii+1));
    dU=(After_U-Before_U);
    dUsave(:,ii) = dU;
    J_prime = (weight.*(dU)); 
    J_prime = rem((J_prime), (2*pi)); 
    Intensity_save(1,ii+1)=After_Beam_Intensity;
    Before_U = After_U;
    After_U = (After_U - J_prime);
    After_U = mod(After_U, (2*pi));
    if (min(min(After_U))<0)
        After_U = After_U - min(min(After_U));
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

    if dJ(1,ii+1) == 0
        break
    end
    
end

figure(2)
set(gcf,'position',[1024,410,560,420]);
clf
subplot(2,2,1);
plot(Target_Intensity_Sum,'b*-'); 
title('beam 중심의 intensity'); 
subplot(2,2,2);
plot(dJ,'go-'); 
title('intensity의 변화량'); 
subplot(2,2,3);
hold on;
plot((Usave(1,:)./(2*pi)),'r');
plot((Usave(2,:)./(2*pi)),'b');
hold off;
title('각 beam의 위상'); 
subplot(2,2,4);
hold on;
plot((dUsave(1,:)./(2*pi)),'r');
plot((dUsave(2,:)./(2*pi)),'b');
hold off;
title('각 beam의 위상의 변화량'); 