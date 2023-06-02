clc
clear all % focused gaussian beam calculator
close all;

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

Intensity_manual_store = zeros(1,5);
count = 0;

Iterration_Count=100;

filename = 'test_reduced2.gif';

M=512;
L1=M*15e-6; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 

test_x=[0:0.001:2*pi]; test_y=[0:0.001:2*pi];
[X,Y] = meshgrid(test_x,test_y); test_z=abs(exp(-1i*X)+exp(-1i*Y)).^2;

%%

for jj = 1:10000

    Before_U = rand(1,2)*(2*pi);
    After_U = rand(1,2)*(2*pi);

    Before_Beam_Flow = exp(-1i*Before_U(1,1))+exp(-1i*Before_U(1,2));
    Usave(:,1) = Before_U;

    Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);

    Target_Intensity_Sum(jj,1) = Before_Beam_Intensity;
    dJ(1,1)= 1;
    dUsave(:,1) = abs(After_U-Before_U);

    for ii = 1:Iterration_Count

        After_Beam_Flow = exp(-1i*After_U(1,1))+exp(-1i*After_U(1,2));
        After_Beam_Intensity = abs(After_Beam_Flow).^2;
        Usave(:,ii+1) = After_U;

        dU = rand(1,2)*(2*pi);

        Target_Intensity_Sum(jj,ii+1) = After_Beam_Intensity;

        dJ(1,ii+1)=((Target_Intensity_Sum(jj,ii+1) - Target_Intensity_Sum(jj,ii)));
        dU = dU.*(After_U-Before_U);
        if (var(dU) == 0) 
            Variance_dU = 0.000001;
        else
            Variance_dU = abs(var(dU));
        end
        Gamma(1,ii) = (1-(ii/1000)+(100/(ii^1.2)))/max(Target_Intensity_Sum(jj,:));%9
        Gamma(isnan(Gamma)) = 1;
        weight=Gamma(1,ii) .* (dJ(1,ii+1))./(Variance_dU);
        dUsave(:,ii) = dU;
        J_prime = (weight.*(dU)); 
        Intensity_save(1,ii+1)=After_Beam_Intensity;
        Before_U = After_U;
        After_U = (After_U + J_prime);
        After_U = mod(After_U, (2*pi));
        for ll = 1:2
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

        if dJ(1,ii+1) == 0
            break
        end

    end

    Intensity_Storing(1,jj) = Target_Intensity_Sum(jj,ii+1);
    IS = Target_Intensity_Sum(jj,ii+1);

    Intensity_manual_store(1,round(IS)+1) = Intensity_manual_store(1,round(IS)+1)+1;

    if (Target_Intensity_Sum(jj,ii+1)>=8.8)
        count = count+1;
    end

end

figure(2)
set(gcf,'position',[1024,410,560,420]);
clf
subplot(2,2,1);
plot(Target_Intensity_Sum(jj,:),'b*-'); 
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
