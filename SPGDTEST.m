clc
clear all

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

load('initset.mat');

Iterration_Count=1000;
Checking_Size=16*4;%16*4;
Super_Pixels=8;

M=512; %256
L1=M*15e-6; %3.84e-3; 
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; 
y1=x1;
lambda=1.064e-6;
k=2*pi/lambda; 
w=dx1*(160/2);

z=180e-3;
zf=z;
z2=8e-3;
zf2=z2;
[X1,Y1]=meshgrid(x1,y1);
[theta,rho]=cart2pol(X1,Y1);

Initial_Beam_Flow=zeros(M,M);
Initial_Beam_Flow=exp(-(rho.^2)/w^2);

perturb = rand(Super_Pixels, Super_Pixels);%(binornd(ones(Super_Pixels,Super_Pixels),ones(Super_Pixels,Super_Pixels)./2)-0.5); %.*(4*pi)

Initial_Beam_Phase=angle(Initial_Beam_Flow);
Initial_Beam_Intensity=(abs(Initial_Beam_Flow).^2);
Initial_Beam_Radius=radCal(Initial_Beam_Intensity);

Jratio_init=1/(sum(sum(Initial_Beam_Intensity))/chk_J(Initial_Beam_Intensity,Checking_Size)-1)*100;

%% plot area

hold off
figure(1)
subplot(4,4,1);
imagesc(x1/1e-3,y1/1e-3,Initial_Beam_Intensity); 
hold on
rectangle('Position',[-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,[num2str(Jratio_init)],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(Initial_Beam_Radius*0.0150),'mm']); 
colorbar
hold off

% figure(2) 
subplot(4,4,2);
imagesc(x1/1e-3,y1/1e-3,Initial_Beam_Phase); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%% SPGD initializing

Before_U = zeros(Super_Pixels,Super_Pixels).*(2*pi);
After_U = rand([Super_Pixels,Super_Pixels]).*(2*pi);
% au = preset1;

Usave(:,:,1) = After_U;
Initial_U = After_U;

Before_U_expand = padarray(expand(Before_U),[128 128],0,'both');

Changing_Beam_Flow=Initial_Beam_Flow;

Before_Beam_Flow=propTF(Initial_Beam_Flow.*exp(-1i*Before_U_expand),L1,lambda,2*z);
Before_Beam_Flow=propTF(Before_Beam_Flow,L1,lambda,2*z);

Before_Beam_Phase=angle(Before_Beam_Flow)/pi+1; 
Before_Beam_Intensity=(abs(Before_Beam_Flow).^2);
Before_Beam_Radius=radCal(Before_Beam_Intensity);

Image_Intensity_Sum(1,1) = sum(sum(Before_Beam_Intensity));
Target_Intensity_Sum(1,1) = chk_J(Before_Beam_Intensity,Checking_Size);
Target_Intensity_Ratio(1,1)=1/(Image_Intensity_Sum(1,1)/Target_Intensity_Sum(1,1)-1)*100;
dJ(1,1)= 1;
MaxValue_U(1,1)=max(max(angle(Before_U)));
Max_Intensity(1,1)=max(max((abs(Before_Beam_Intensity).^2)));
Min_Intensity(1,1)=min(min((abs(Before_Beam_Intensity).^2)));
dUsave(:,:,1) = (After_U-Before_U);

%% plot area

hold off
% figure(3)
subplot(4,4,3);
imagesc(x1/1e-3,y1/1e-3,Before_Beam_Intensity);
hold on
rectangle('Position',[-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,[num2str(Target_Intensity_Ratio(1,1))],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(Before_Beam_Radius*0.0150),'mm']); 
colorbar
hold off

% figure(4) 
subplot(4,4,4);
imagesc(x1/1e-3,y1/1e-3,Before_Beam_Phase); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%% SPGD Loop

for ii = 1:Iterration_Count
    
    Usave(:,:,ii+1) = After_U;
    After_U_expand = padarray(expand(After_U),[128 128],0,'both');
    
    perturb = rand(Super_Pixels, Super_Pixels).*(2*pi);

    After_Beam_Flow=propTF(Changing_Beam_Flow.*exp(-1i*After_U_expand),L1,lambda,2*z);
    After_Beam_Flow=propTF(After_Beam_Flow,L1,lambda,2*z);
    After_Beam_Intensity = abs(After_Beam_Flow).^2;
    Intensity_save(:,:,ii) = After_Beam_Intensity;
    
    Image_Intensity_Sum(1,ii+1) = sum(sum((After_Beam_Intensity)));
    Target_Intensity_Sum(1,ii+1) = chk_J((After_Beam_Intensity),Checking_Size);
    Target_Intensity_Ratio(1,ii+1)=1/(Image_Intensity_Sum(1,ii+1)/Target_Intensity_Sum(1,ii+1)-1)*100;
    dJ(1,ii+1)=((Target_Intensity_Sum(1,ii+1) - Target_Intensity_Sum(1,ii)));
    dU=(After_U-Before_U);
    if (var(var(dU)) == 0)
        Variance_dU = 0.000001;
    else
        Variance_dU = abs(var(var(dU)))./(ii);
    end
    Gamma(1,ii) = max(Target_Intensity_Sum)/Target_Intensity_Sum(1,ii+1);  %(1-(ii/1000)+(400/(ii^1.25)))
    Weight=Gamma(1,ii).*(dJ(1,ii+1))./(Variance_dU);
    J_prime = (Weight.*(perturb));
    Max_Intensity(1,ii+1)=max(max((After_Beam_Intensity)));
    Min_Intensity(1,ii+1)=min(min((After_Beam_Intensity)));
    Before_U = After_U;
    After_U = (After_U + J_prime);
    MaxValue_U(1,ii+1)=max(max(After_U))/(2*pi);
    for kk = 1:Super_Pixels
        for ll = 1:Super_Pixels
            if (After_U(kk,ll) >= (2*pi))
                After_U(kk,ll) = rem(After_U(kk,ll), (2*pi));
            end
            if (After_U(kk,ll) < 0)
                After_U(kk,ll) = After_U(kk,ll) - (2*pi)*(fix(After_U(kk,ll)/(2*pi)));
            end
        end
    end
    if abs(dJ(1,ii+1)) <= 1e-24
        break
    end
    
    %% plot area
    
    if ii==1
        hold off
        After_Beam_Phase=angle(After_Beam_Flow)/pi+1; 
        After_Beam_Intensity=(abs(After_Beam_Flow).^2);
%         figure(11) 
        subplot(4,4,5);
        imagesc(x1/1e-3,y1/1e-3,After_Beam_Intensity);
        hold on
        rectangle('Position',[-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
        text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,[num2str(Target_Intensity_Ratio(1,ii+1))],'Color','red');
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('result'); 
        colorbar
        hold off
        
        subplot(4,4,6);
        imagesc(x1/1e-3,y1/1e-3,After_Beam_Phase); 
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('align'); 
        colorbar
    end
    figure(200)
    
    hold off
    subplot(1,2,1);
    imagesc(x1/1e-3,y1/1e-3,After_Beam_Intensity);
    hold on
    rectangle('Position',[-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
    text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,[num2str(Target_Intensity_Ratio(1,ii+1))],'Color','red');
    axis square; axis xy; 
    colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
    title('result'); 
    colorbar
    hold off
    
    subplot(1,2,2);
    imagesc(Before_U); 
    axis square; axis xy; 
    colormap('gray');
    title('U'); 
    colorbar
    drawnow
end

After_Beam_Phase=angle(After_Beam_Flow)/pi+1; 
After_Beam_Radius=radCal(After_Beam_Intensity);

%% plot area
figure(1)

hold off
% figure(5) 
subplot(4,4,7);
imagesc(x1/1e-3,y1/1e-3,After_Beam_Intensity);
hold on
rectangle('Position',[-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,[num2str(Target_Intensity_Ratio(1,ii+1))],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('result'); 
colorbar
hold off

% figure(6)
subplot(4,4,8);
imagesc(x1/1e-3,y1/1e-3,After_Beam_Phase); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

% figure(7)  
subplot(4,4,9);
imagesc(Initial_U); 
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(7)  
subplot(4,4,10);
imagesc(Before_U); 
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(8)                                                                                         % rescaled figure of J, delta J, change amount of au, cycle of au
subplot(4,4,11);
hold off
plot(rescale(Target_Intensity_Sum),'b*-');                                                                             % cacluated J value during iterration
hold on
plot(rescale(dJ),'go-');                                                                             % improvement of output during iterration
title('J'); 
legend('J','del J');%
hold off

% figure(9)
subplot(4,4,12);
plot(Max_Intensity,'b*-'); 
title('Max'); 
subplot(4,4,13);
plot(Min_Intensity,'ro-'); 
title('Min'); 

% figure(10) 
subplot(4,4,14);
returnVal = padarray(expand(round(Before_U*((2^Super_Pixels)/max(max(Before_U))))),[128 128],0,'both');
returnVal = uint16(returnVal);
imagesc(returnVal); 
axis square; axis xy; 
colormap('gray');
colorbar
imwrite(returnVal,"ExportImage.tiff",'tiff');

% figure(12) 
subplot(4,4,15);
hold off
plot(Target_Intensity_Sum,'b*-'); 
hold on
title('J');  xlabel('iteration count'); ylabel('J value'); 
hold off

% figure(13) 
subplot(4,4,16);
hold off
plot(Target_Intensity_Ratio,'b*-'); 
hold on
title('intensity ratio');  xlabel('iteration count'); ylabel('percentage'); 
hold off

figure(2) 
subplot(2,2,1);
plot(Target_Intensity_Sum,'b*-'); 
title('J'); 
subplot(2,2,2);
plot(dJ,'go-'); 
title('dJ'); 
subplot(2,2,4);
plot(Gamma,'y^-'); 
title('gamma'); 