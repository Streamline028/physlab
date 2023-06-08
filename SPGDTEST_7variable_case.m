clc
clear all % focused gaussian beam calculator
close all;

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

mode = 7;
Intensity_manual_store = zeros(1,((mode^2)+1));
count = 0;

Iterration_Count = 100;
total_Iterration = 2000;

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
    After_U = rand(1,mode)*(2*pi);
    Before_Beam_Flow = 0;
    
    for pp=1:mode
        Before_Beam_Flow = Before_Beam_Flow + exp(-1i*Before_U(1,pp));
    end
    
    Usave(:,1) = Before_U;

    Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

    Target_Intensity_Sum(jj,1) = Before_Beam_Intensity;
    dJ(1,1)= 1;
    dUsave(:,1) = abs(After_U-Before_U);
    ii = 1;
    while (dJ(1,ii) ~= 0) %for ii = 1:Iterration_Count
        
        After_Beam_Flow = 0;
        for pp=1:mode
            After_Beam_Flow = After_Beam_Flow + exp(-1i*After_U(1,pp));
        end
        
        After_Beam_Intensity = abs(After_Beam_Flow).^2;
        Usave(:,ii+1) = After_U;

        dU = rand(1,mode)*(2*pi);

        Target_Intensity_Sum(jj,ii+1) = After_Beam_Intensity;

        dJ(1,ii+1)=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
        dU = dU.*(After_U-Before_U);
        if (var(dU) == 0) 
            Variance_dU = 0.000001;
        else
            Variance_dU = abs(var(dU));
        end
        Gamma(1,ii) = (1+((ii)/10)-(1/(ii^10)))/max(Target_Intensity_Sum(jj,:));%9
        Gamma(isnan(Gamma)) = 1;
        weight=Gamma(1,ii) .* (dJ(1,ii+1))./(Variance_dU);
        dUsave(:,ii) = dU;
        J_prime = (weight.*(dU)); 
        Intensity_save(1,ii+1)=After_Beam_Intensity;
        Before_U = After_U;
        After_U = (After_U + J_prime);
        After_U = mod(After_U, (2*pi));
        for ll = 1:mode
            if (After_U(1,ll) >= (2*pi))
                After_U(1,ll) = (After_U(1,ll) - (2*pi));
            end
            if (After_U(1,ll) < 0)
                After_U(1,ll) = After_U(kk,ll) - (2*pi)*(fix(After_U(1,ll)/(2*pi))-1);
            end
        end
    % %     figure(100);
    % %     set(gcf,'position',[456,411,560,420]);
    % %     clf;
    % %     imagesc(test_x,test_y,test_z);
    % %     line(0.1*cos(0:0.01:2*pi)+After_U(1),0.1*sin(0:0.01:2*pi)+After_U(2),'Color','r','Linewidth',2);
    % %     axis image; axis xy;
    % %     colormap('gray')
    % %     title(sprintf("ii = %d, J = %.2f\nx = %.2f, y = %.2f, d\\phi = %.3f",ii,Target_Intensity_Sum(jj,ii),After_U(1)/2/pi,After_U(2)/2/pi,(After_U(2)-After_U(1))/2/pi))
    % %     xlabel('x');ylabel('y');
    % %     xlim([0 2*pi]); ylim([0 2*pi]);
    % %     drawnow
    % %     frame = getframe(100); 
    % %     img = frame2im(frame);
    % %     [imind cm] = rgb2ind(img,256); 
    % %     if ii == 1
    % %         imwrite(imind,cm,filename,'gif','Loopcount',1,'DelayTime',1/5); 
    % %     else
    % %         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/5); 
    % %     end
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
%         figure(200+jj)
%         set(gcf,'position',[1024,410,560,420]);
%         clf
%         subplot(2,2,1);
%         plot(Target_Intensity_Sum(jj,:),'b*-'); 
%         title('beam 중심의 intensity'); 
%         subplot(2,2,2);
%         plot(dJ,'go-'); 
%         title('intensity의 변화량'); 
%         subplot(2,2,3);
%         hold on;
%         for pp=1:mode
%             plot((Usave(mode,:)./(2*pi)),'r');
%         end
%         hold off;
%         title('각 beam의 위상'); 
%         subplot(2,2,4);
%         hold on;
%         for pp=1:mode
%             plot((dUsave(mode,:)./(2*pi)),'r');
%         end
%         hold off;
%         title('각 beam의 위상의 변화량'); 
% 
%         if dJ(1,ii+1) == 0
%             pause(2);
%             break
%         end
        ii = ii+1;

    end

    Intensity_Storing(1,jj) = Target_Intensity_Sum(jj,ii);
    IS = Target_Intensity_Sum(jj,ii);

    Intensity_manual_store(1,round(IS)+1) = Intensity_manual_store(1,round(IS)+1)+1;

    if (Target_Intensity_Sum(jj,ii)>=(mode^2)-0.2)
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

100*count/jj