clc
clear all

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

Intensity_manual_store = zeros(1,5);
count = 0;

Iterration_Count=300;
init_phase = pi;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

figure(1); test_x=[0:0.001:4]; test_y=abs(exp(-1i*test_x*(2*pi))+exp(-1i*init_phase)).^2; plot(test_x, test_y);
hold on; for (ij = 1:size(test_x,2)-1) test_dy(ij)=(test_y(ij+1)-test_y(ij)); test_dy(ij)= test_dy(ij).*test_dy(ij); test_dx(ij) = test_x(ij); end; plot(test_dx, test_dy*10000+2, 'r'); hold off;

%%

for jj = 1:10000

Before_U = 0.2*(2*pi);
After_U = rand()*(2*pi);

Before_Beam_Flow = exp(-1i*init_phase)+exp(-1i*Before_U);
Usave(1,1) = Before_U;
 
Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

Target_Intensity_Sum(jj,1) = Before_Beam_Intensity;
dJ(1,1)= 1;
dUsave(:,:,1) = abs(After_U-Before_U);

for ii = 1:Iterration_Count
    
    After_Beam_Flow = exp(-1i*After_U)+exp(-1i*init_phase);
    After_Beam_Intensity = abs(After_Beam_Flow).^2;
    Usave(1,ii+1) = After_U;
    
    Target_Intensity_Sum(jj,ii+1) = After_Beam_Intensity;
    
    dJ(1,ii+1)=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
    dU(1,ii)=(After_U-Before_U);
    if (var(dU) == 0) 
        Variance_dU = 0.000001;
    else
        Variance_dU = abs(var(dU))./(ii);
    end
    Gamma(1,ii) = (1-(ii/1000)+(50/(ii^1.2)))/max(Target_Intensity_Sum(jj,:));
    weight=Gamma(1,ii)*(dJ(1,ii+1))./(Variance_dU);
    dUsave(1,ii) = dU(1,ii);
    J_prime = (weight.*(dU(1,ii))); 
    if (J_prime >= (2*pi))
        J_prime = J_prime - (2*pi);
    end
    Intensity_save(1,ii+1)=After_Beam_Intensity;
    Before_U = After_U;
    After_U = (After_U + J_prime);
    if (After_U >= (2*pi))
        After_U = After_U - (2*pi);
    end
    if abs(dJ(1,ii+1)) <= 1e-6
        break
    end
    
end

Intensity_Storing(1,jj) = Target_Intensity_Sum(jj,ii+1);
IS = Target_Intensity_Sum(jj,ii+1);

Intensity_manual_store(1,round(IS)+1) = Intensity_manual_store(1,round(IS)+1)+1;

if (Target_Intensity_Sum(jj,ii+1)>=3.8)
    count = count+1;
end

end

Store = sort(Intensity_Storing);

figure(2) 
subplot(2,2,1);
plot(Target_Intensity_Sum(jj,:),'b*-'); 
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

figure(20) 
cdfplot(Intensity_Storing);
title([num2str(jj),' iteration']); 

figure(30) 
bar([0:4],Intensity_manual_store); 
ylim([0 jj]);
title([num2str(jj),' iteration']); 
text(-0.3, Intensity_manual_store(1,1)+500,[num2str((Intensity_manual_store(1,1)/jj)*100),'%'],'Color','m');
text(0.7, Intensity_manual_store(1,2)+500,[num2str((Intensity_manual_store(1,2)/jj)*100),'%'],'Color','m');
text(1.7, Intensity_manual_store(1,3)+500,[num2str((Intensity_manual_store(1,3)/jj)*100),'%'],'Color','m');
text(2.7, Intensity_manual_store(1,4)+500,[num2str((Intensity_manual_store(1,4)/jj)*100),'%'],'Color','m');
text(3.7, Intensity_manual_store(1,5)+500,[num2str((Intensity_manual_store(1,5)/jj)*100),'%'],'Color','m');

100*count/jj
