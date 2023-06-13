clc
clear
close all;

addpath(genpath('.\testing_result'))

load(strcat('.\testing_result\0613\individual\',num2str(10),'.mat'),'mode');
load(strcat('.\testing_result\0613\individual\',num2str(10),'.mat'),'jj');
modval = [10 20 30 40 50 60];

for ii=10:10:60
    percentage(1,ii/10) = load(strcat('.\testing_result\0613\individual\',num2str(ii),'.mat'),'ans');
    meanVal(1,ii/10) = load(strcat('.\testing_result\0613\individual\',num2str(ii),'.mat'),'meaniter');
    variVal(1,ii/10) = load(strcat('.\testing_result\0613\individual\',num2str(ii),'.mat'),'variter');
    TimeVal(1,ii/10) = load(strcat('.\testing_result\0613\individual\',num2str(ii),'.mat'),'meanIterationTime');
    pureTimeVal(1,ii/10) = load(strcat('.\testing_result\0613\individual_pure\',num2str(ii),'.mat'),'meanIterationTime');
end
percentage = table2array(struct2table(percentage));
meanVal = table2array(struct2table(meanVal));
variVal = table2array(struct2table(variVal));
TimeVal = table2array(struct2table(TimeVal));
pureTimeVal = table2array(struct2table(pureTimeVal));

pTV=pureTimeVal(6,1)-pureTimeVal(1,1);
TV=TimeVal(6,1)-TimeVal(1,1);

channel_per_usingtime = pureTimeVal.*meanVal./[10;20;30;40;50;60]

figure(39) 
plot(modval,percentage,'o-'); 
ylim([0 100]);
title([num2str(jj),' iteration']); 
xlabel('mode'); ylabel('%');

figure(300) 
bar(modval,TimeVal); 
title([num2str(jj),' iteration']); 
xlabel('mode'); ylabel('time');

figure(400) 
bar(modval,pureTimeVal); 
title([num2str(jj),' iteration']); 
xlabel('mode'); ylabel('time');

figure(301) 
hold on
bar(modval,meanVal); 
errorbar(modval,meanVal,variVal,"o"); 
hold off
ylim([0 2100]);
title([num2str(jj),' iteration']); 
xlabel('mode'); ylabel('iteration count');

withGamma = load(strcat('.\testing_result\0613\individual\withGamma.mat'),'Target_Intensity_Sum');
withGamma = load(strcat('.\testing_result\0613\individual\withGamma.mat'),'weight_save');
withoutGamma = load(strcat('.\testing_result\0613\individual\withoutGamma.mat'),'Target_Intensity_Sum');
withoutGamma = load(strcat('.\testing_result\0613\individual\withoutGamma.mat'),'weight_save');
