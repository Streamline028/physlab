clc
clear all

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

Iterration_Count=100;
init_phase =pi;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

figure(1); test_x=[0:0.001:4]; test_y=abs...
(exp(-1i*test_x*(2*pi))+exp(-1i*init_phase)).^2; plot(test_x, test_y);
hold on; for (ij = 1:size(test_x,2)-1) test_dy(ij)=(test_y(ij+1)-test_y(ij)); test_dy(ij)= test_dy(ij).*test_dy(ij); test_dx(ij) = test_x(ij); end; plot(test_dx, test_dy*10000+2, 'r'); hold off;

%%

Before_U = 0.2*(2*pi);
After_U = rand()*(2*pi);

Before_Beam_Flow = exp(-1i*init_phase)+exp(-1i*Before_U);
Usave(1,1) = Before_U;
 
Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

Target_Intensity_Sum(1,1) = Before_Beam_Intensity;
dJ(1,1)= 1;
dUsave(:,:,1) = abs(After_U-Before_U);

for ii = 1:Iterration_Count
    
    After_Beam_Flow = exp(-1i*After_U)+exp(-1i*init_phase);
    After_Beam_Intensity = abs(After_Beam_Flow).^2;
    Usave(1,ii+1) = After_U;
    
    Target_Intensity_Sum(1,ii+1) = After_Beam_Intensity;
    
    dJ(1,ii+1)=((Target_Intensity_Sum(1,ii+1) - Target_Intensity_Sum(1,ii)));
    Gamma(1,ii) = 2.6/Target_Intensity_Sum(1,ii+1);
    Gamma(isnan(Gamma)) = 1;
    weight=Gamma(1,ii)*(dJ(1,ii+1));
    dU=(After_U-Before_U);
    dUsave(1,ii) = dU;
    J_prime = (weight.*(dU)); 
    J_prime = rem((J_prime), (2*pi));  
    Intensity_save(1,ii+1)=After_Beam_Intensity;
    Before_U = After_U;
    After_U = (After_U - J_prime);
    After_U = rem(After_U, (2*pi));
    if (min(min(After_U))<0)
        After_U = After_U - min(min(After_U));
    end
    if dJ(1,ii+1) == 0 
        break
    end
    
end

figure(2) 
subplot(2,2,1);
plot(Target_Intensity_Sum,'b*-'); 
title('beam 중심의 intensity'); 
subplot(2,2,2);
plot(dJ,'go-'); 
title('intensity의 변화량'); 
subplot(2,2,3);
plot(Usave./(2*pi),'r^-');
title('바꾸어 주는 빛의 위상'); 
subplot(2,2,4);
plot(dUsave,'y^-'); 
title('위상의 변화량'); 

figure(10) 
plot(Gamma,'*-'); 
title('gamma'); 