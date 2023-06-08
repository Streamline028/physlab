clc
clear all % focused gaussian beam calculator
close all;

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

mode = 8;
Intensity_manual_store = zeros(1,((mode^2)+1));
count = 0;

Iterration_Count = 1000;
total_Iterration = 100;

filename = 'test_reduced4.gif';

M=512;
L1=M*15e-6; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 

%%

for jj = 1:total_Iterration

    Before_U = rand(1,mode)*(2*pi);
    Before_Beam_Flow = 0;
    
    for pp=1:mode
        Before_Beam_Flow = Before_Beam_Flow + exp(-1i*Before_U(1,pp));
    end
    
    Usave(:,1) = Before_U;

    Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

    Target_Intensity_Sum(jj,1) = Before_Beam_Intensity;
    dJ(1,1)= 1;
    dUsave(:,1) = Before_U;

    for ii = 1:Iterration_Count

        After_U = rand(1,mode)*(2*pi);
        
        After_Beam_Flow = 0;
        for pp=1:mode
            After_Beam_Flow = After_Beam_Flow + exp(-1i*After_U(1,pp));
        end
        
        After_Beam_Intensity = abs(After_Beam_Flow).^2;
        Usave(:,ii+1) = After_U;

        Target_Intensity_Sum(jj,ii+1) = After_Beam_Intensity;

        dJ(1,ii+1)=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
        dU = (After_U-Before_U);
        if (var(dU) == 0) 
            Variance_dU = 0.000001;
        else
            Variance_dU = abs(var(dU));
        end
        Gamma(1,ii) = (1-((ii)/1000)+(100/(ii^1.2)))/max(Target_Intensity_Sum(jj,:));%9
        Gamma(isnan(Gamma)) = 1;
        weight=Gamma(1,ii) .* (dJ(1,ii+1))./(Variance_dU);
        dUsave(:,ii) = (After_U-Before_U);
        J_prime = (weight.*(dU)); 
        Intensity_save(1,ii+1)=After_Beam_Intensity;
        Before_U = (Before_U + J_prime);
        Before_U = mod(Before_U, (2*pi));
        for ll = 1:mode
            if (Before_U(1,ll) < 0)
                Before_U(1,ll) = Before_U(1,ll) - (2*pi)*(fix(Before_U(1,ll)/(2*pi))-1);
            end
        end
        
        Before_Beam_Flow = 0;
        for pp=1:mode
            Before_Beam_Flow = Before_Beam_Flow + exp(-1i*Before_U(1,pp));
        end
        Before_Beam_Intensity = abs(Before_Beam_Flow).^2;
        Usave(:,ii+1) = Before_U;
        Target_Intensity_Sum(jj,ii+1) = Before_Beam_Intensity;
        
        if dJ(1,ii+1) == 0
            break
        end

    end

    Intensity_Storing(1,jj) = Target_Intensity_Sum(jj,ii+1);
    IS = Target_Intensity_Sum(jj,ii+1);

    Intensity_manual_store(1,round(IS)+1) = Intensity_manual_store(1,round(IS)+1)+1;

    if (Target_Intensity_Sum(jj,ii+1)>=(mode^2)-0.2)
        count = count+1;
    end

end

figure(2)
set(gcf,'position',[1024,410,560,420]);
clf
subplot(2,2,1);
plot(Target_Intensity_Sum(jj,:),'b*-'); 
xlim([0 ii+1]);
title('beam 중심의 intensity'); 
xlabel('반복 횟수'); ylabel('intensity');
subplot(2,2,2);
plot(dJ,'go-'); 
xlim([0 ii+1]);
title('intensity의 변화량'); 
xlabel('반복 횟수'); ylabel('intensity');
subplot(2,2,3);
hold on;
for pp=1:mode
    plot((Usave(pp,:)));
end
xlim([0 ii+1]);
hold off;
title('각 beam의 위상'); 
xlabel('반복 횟수'); ylabel('Phase');
subplot(2,2,4);
hold on;
for pp=1:mode
    plot((dUsave(pp,:)));
end
xlim([0 ii+1]);
hold off;
title('각 beam의 위상의 변화량'); 
xlabel('반복 횟수'); ylabel('Phase');

figure(20) 
cdfplot(Intensity_Storing);
title([num2str(jj),' iteration']); 
xlabel('intensity'); ylabel('누적 발견 횟수');

figure(30) 
bar([0:mode^2],Intensity_manual_store); 
ylim([0 jj]);
title([num2str(jj),' iteration']); 
xlabel('intensity'); ylabel('발견 횟수');
for qq = 1:((mode^2)+1)
    text(qq-1.3, Intensity_manual_store(1,qq)+(jj/20),[num2str((Intensity_manual_store(1,qq)/jj)*100),'%'],'Color','m');
end

figure(52)
plot(Gamma,'b*-'); 
xlim([0 ii+1]);
title('Gamma'); 
xlabel('반복 횟수'); ylabel('Gamma');

100*count/jj