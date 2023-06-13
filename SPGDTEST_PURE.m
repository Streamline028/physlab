clc
clear all % focused gaussian beam calculator
close all;

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

iterval = [750 1000 1250 1500 1750 2000 2250 2500];
modval = [60 50 40 30 20 10 20 30 40 50 60 70 80];

for LL = 1:6
    Iterration_Count = 2000;%iterval(1,LL);
    total_Iterration = 1000;
    mode = modval(1,LL);

    Intensity_manual_store = zeros(1,((mode^2)+1));
    Iterration_save_dummy = zeros(Iterration_Count, ((mode^2)+1));
    Iterration_save = zeros(1,Iterration_Count);
    Iteration_save = 0;
    count = 0;

    %%

    for jj = 1:total_Iterration

        clear dJ dU

        Before_U = rand(1,mode)*(2*pi);
        Before_Beam_Flow = 0;

        for pp=1:mode
            Before_Beam_Flow = Before_Beam_Flow + exp(-1i*Before_U(1,pp));
        end

        Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

        Target_Intensity_Sum(jj,1) = Before_Beam_Intensity;
        dJ= 1;

        for ii = 1:Iterration_Count
            tic
            After_U = rand(1,mode)*(2*pi);

            After_Beam_Flow = 0;
            for pp=1:mode
                After_Beam_Flow = After_Beam_Flow + exp(-1i*After_U(1,pp));
            end

            After_Beam_Intensity = abs(After_Beam_Flow).^2;

            Target_Intensity_Sum(jj,ii+1) = After_Beam_Intensity;

            dJ=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
            dU = (After_U-Before_U);
            if (var(dU) == 0) 
                Variance_dU = 0.000001;
            else
                Variance_dU = abs(var(dU));
            end
            Gamma = (1-((ii)/Iterration_Count)+(100/(ii^1.2)))/max(Target_Intensity_Sum(jj,:));%9
            weight=Gamma .* (dJ)./(Variance_dU);
            J_prime = (weight.*(dU)); 
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
            Target_Intensity_Sum(jj,ii+1) = Before_Beam_Intensity;
            dJ=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
            T(jj,ii)=toc;
            if abs(dJ) <= 1e-4
                break
            end

    %         if (Before_Beam_Intensity >= ((mode^2)*0.99))
    %             break
    %         end

        end

        Intensity_Storing(1,jj) = Target_Intensity_Sum(jj,ii+1);
        IS = Target_Intensity_Sum(jj,ii+1);

        Intensity_manual_store(1,round(IS)+1) = Intensity_manual_store(1,round(IS)+1)+1;

    %     Iterration_save_dummy(ii,round(IS)+1) = Iterration_save_dummy(ii,round(IS)+1) + 1;
    %     Iterration_save(1,ii) = Iterration_save(1,ii) + 1;

        if (Target_Intensity_Sum(jj,ii+1)>=(mode^2)*0.9)
            Iterration_save_dummy(ii,round(IS)+1) = Iterration_save_dummy(ii,round(IS)+1) + 1;
            Iterration_save(1,ii) = Iterration_save(1,ii) + 1;
            Iteration_save(1, count+1) = ii;
            count = count+1;
        end

    end
    figure(1200+LL)
    plot(Target_Intensity_Sum(jj,:),'b*-'); 
    xlim([0 ii+1]);
    title('beam 중심의 intensity'); 
    xlabel('반복 횟수'); ylabel('intensity');
    figure(1300+LL)
    plot(dJ,'go-'); 
    xlim([0 ii+1]);
    title('intensity의 변화량'); 
    xlabel('반복 횟수'); ylabel('intensity');

    figure(3000+LL) 
    bar([0:mode^2],Intensity_manual_store); 
    ylim([0 jj]);
    title([num2str(jj),' iteration']); 
    xlabel('intensity'); ylabel('발견 횟수');
    
    save_max_percent(1,LL) = count;
    variter(1,LL) = std(Iteration_save);
    meaniter(1,LL) = mean(Iteration_save);
    meanIterationTime(LL) = mean(mean(T));
end

figure(2)
plot(Target_Intensity_Sum(jj,:),'b*-'); 
xlim([0 ii+1]);
title('beam 중심의 intensity'); 
xlabel('반복 횟수'); ylabel('intensity');

figure(20) 
cdfplot(Intensity_Storing);
title([num2str(jj),' iteration']); 
xlabel('intensity'); ylabel('누적 발견 횟수');

figure(30) 
bar([0:mode^2],Intensity_manual_store); 
ylim([0 jj]);
title([num2str(jj),' iteration']); 
xlabel('intensity'); ylabel('발견 횟수');

figure(31) 
bar([1:Iterration_Count],Iterration_save); 
ylim([0, max(Iterration_save)+10]);
title([num2str(jj),' iteration']); 
xlabel('iteration'); ylabel('count');

figure(39) 
plot(modval(1:LL),save_max_percent,'o-'); 
ylim([0 1000]);
title([num2str(jj),' iteration']); 
xlabel('mode'); ylabel('발견 횟수');

figure(300) 
bar(modval(1:LL),meanIterationTime); 
title([num2str(jj),' mode']); 
xlabel('iteration'); ylabel('time');

figure(301) 
hold on
bar(modval(1:LL),meaniter); 
errorbar(modval(1:LL),meaniter,variter,"o"); 
hold off
ylim([0, Iterration_Count+100]);
title([num2str(jj),' mode']); 
xlabel('iteration'); ylabel('iteration count');

100*count/jj
% 
% figure(39) 
% plot(modval,save_max_percent,'o-'); 
% ylim([0 1000]);
% xlim([10 90]);
% title([num2str(jj),' iteration']); 
% xlabel('mode'); ylabel('발견 횟수');