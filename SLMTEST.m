clear all % focused gaussian beam calculator
close all
clc

% addpath(genpath('C:\Users\a0105\Documents\MATLAB\TM'))

%% SLM calibration
BNS_OpenSLMs()                                                                                      % SLM 함수 열기
BNS_LoadLUTFile(1,'C:\Users\a0105\Documents\MATLAB\TM\LUT_Files\slm5240_at1064_P16.lut')            % SLM calibration LUT 함수(1064 nm용) 삽입
Image.optimization_data = libpointer('uint8Ptr', zeros(524288,1));                                  % SLM image 삽입 포인터 지정 (Calibration용)
BNS_ReadTIFF('C:\PCIe16MatlabSDK\Image_Files\White.tiff', Image.optimization_data);                 % Calibration 용 white 이미지 삽입 (2 pi phase (2^16) 전체 삽입)
BNS_WriteCal(Image.optimization_data, Image);                                                       % Calibration
%% SLM on
BNS_SetPower(true)                                                                                  % SLM 전원 on
Image.Image0 = libpointer('uint8Ptr', zeros(524288,1));                                             % SLM 이미지 데이터 삽입용 포인터 지정
name_initial = ['C:\Users\a0105\Documents\MATLAB\TM\Tilt_basis_200_waves_DR_16_bits.tiff'];       
% SLM 이미지 파일 사전 지정 -> 초기 Blank 이미지 삽입
BNS_ReadTIFF(name_initial,Image.Image0);                                                            % 이미지 포인터에 삽입
BNS_WriteImage(Image.Image0, Image);                                                                % Calibration
%% Fiber input carmera setup
Camera = TM_camera_renewal;                                                                         % image acquisition toolbox에서 Export한 m-file
%% Basis dir
Basis_folder = uigetdir('C:\Users\a0105\Documents\MATLAB\TM\');
% Basis_folder = 'C:\Users\a0105\Documents\MATLAB\TM\Square_basis_352_32x32';
Basis_folder_name = [Basis_folder,'\'];
%% initial setting
load('initset.mat');                                                                                % get some matrix for testing

r=1000;                                                                                             % iterration count
S=16;                                                                                               % searching size
MM=16;                                                                                              % SLM superpixel count

M=512;                                                                                              % screen size (pixel)
L1=M*15e-6;                                                                                         % screen size (m)
dx1=L1/M;                                                                                           % pixel interval (m)
x1=-L1/2:dx1:L1/2-dx1;                                                                              % generate x field
y1=x1;                                                                                              % generate y field
lambda=1.064e-6;                                                                                    % laser wave length (m)
k=2*pi/lambda;                                                                                      % wave number
w=dx1*(160/2);                                                                                      % laser radius

z=180e-3;                                                                                           % first lens distance
zf=z;                                                                                               % first lens's focal length
z2=8e-3;                                                                                            % second lens distance
zf2=z2;                                                                                             % second lens's focal length
%% initial beam

[X1,Y1]=meshgrid(x1,y1);                                                                            % generate xy field
[theta,rho]=cart2pol(X1,Y1);                                                                        % cartesian coordinate to spherical coordinate

u1=zeros(M,M);                                                                                      % generate zero field
u1=exp(-(rho.^2)/w^2);                                                                              % generate gaussian beam

P1=angle(u1);                                                                                       % calculate phase
I1=(abs(u1).^2);                                                                                    % extract image
R1=radCal(I1);                                                                                      % calculate radius

Jratio_init=1/(sum(sum(I1))/chk_J(I1,S))*100;                                                     % ratio of total to target

hold off                                                                                            % initial gaussian image
figure(1)
subplot(4,4,1);
imagesc(x1/1e-3,y1/1e-3,I1); 
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,['J=',num2str(Jratio_init)],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R1*0.0150),'mm']); 
colorbar
hold off

% figure(2)                                                                                         % phase of initial gaussian image
subplot(4,4,2);
imagesc(x1/1e-3,y1/1e-3,P1); 
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('align'); 
colorbar

%% before SPGD

tic;
Shot_initial = getsnapshot(Source.Parent);                                  % 초기 speckle image
Shot_merge_initial = resize(picavg(Shot_initial,screen_super_pixel_width), screen_super_pixel_width);% super pixel 단위 카메라 output 사전 지정
start(Camera)

bu = zeros(MM,MM).*(2*pi);                                                                          % generate initial before u (since then bu)
au = rand([MM,MM]).*(2*pi);                                                                         % generate initial after u (since then au)
% au = preset6;

bu_ = padarray(expand(bu),[128 128],0,'both');                                                      % generate SLM mask

uu=u1;                                                                                              % save before SLM phase

Image.Image0 = SM_tiffencoder(expand(bu));
%
for i=1:1:Number_of_super_pixel*4
    Basis_name=[Basis_folder_name,'Basis_',num2str(i),'.tiff'];                                     % SLM image 위치

    BNS_ReadTIFF(Basis_name,Image.Image0);                                        % SLM에 image 불러옴
    BNS_WriteImage(Image.Image0, Image);                                    % 만들어진 Hadamard basis matrix를 SLM에 삽입
    pause(0.03);                                                            % SLM 및 카메라 속도 조절
    BNS_WriteImage(Image.Image0, Image);                                    % 만들어진 Hadamard basis matrix를 SLM에 삽입
    pause(0.03);                                                            % SLM 및 카메라 속도 조절
    trigger(Camera)                                                         % Camera 찍기
    pause(str2num(Camera_parameter{1})/1000);
end
Shot_total_init = getdata(Camera,(Number_of_SLM_super_pixel_length^2)*4);
%

% u5=propTF(u1.*exp(-1i*bu_),L1,lambda,2*z);                                                          % simulate SLM
% [u5]=focus(u5,L1,lambda,zf);
% u5=propTF(u5,L1,lambda,2*z);

% P5=angle(u5)/pi+1;                                                                                  % calculate phase
I5=Shot_total_init;                                                                                    % extract image
R5=radCal(I5);                                                                                      % calculate radius

JJ(1,1) = sum(sum(I5));                                                                             % calculate total value
J(1,1) = chk_J(I5,S);                                                                               % calculate target value
Jratio(1,1)=1/(JJ(1,1)/J(1,1))*100;                                                                 % ratio of total to target
W(1,1)= 1;                                                                                          % value of change in J - initially 0
D(1,1)= 1;                                                                                          % expectation of improvement - initially 0
MU(1,1)=max(max(angle(bu)))/(2*pi);                                                                 % calculate bu's cycle
MaxI(1,1)=max(max((abs(I5).^2)));                                                                   % find initial image's maximum value
MinI(1,1)=min(min((abs(I5).^2)));                                                                   % find initial image's minimum value

hold off                                                                                            % before SPGD image
% figure(3)
subplot(4,4,3);
imagesc(x1/1e-3,y1/1e-3,I5);
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,['J=',num2str(Jratio(1,1))],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(R5*0.0150),'mm']); 
colorbar
hold off

% figure(4)                                                                                         % phase of before SPGD image
% subplot(4,4,4);
% imagesc(x1/1e-3,y1/1e-3,P5); 
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title('align'); 
% colorbar
%% SPGD Loop

for ii = 1:r                                                                                        % run until changement is 0
    
    au_ = padarray(expand(au),[128 128],0,'both');                                                  % generate SLM mask
    Image.Image0 = SM_tiffencoder(expand(au));
    %
    for i=1:1:Number_of_super_pixel*4
        Basis_name=[Basis_folder_name,'Basis_',num2str(i),'.tiff'];                                     % SLM image 위치

        BNS_ReadTIFF(Basis_name,Image.Image0);                                        % SLM에 image 불러옴
        BNS_WriteImage(Image.Image0, Image);                                    % 만들어진 Hadamard basis matrix를 SLM에 삽입
        pause(0.03);                                                            % SLM 및 카메라 속도 조절
        BNS_WriteImage(Image.Image0, Image);                                    % 만들어진 Hadamard basis matrix를 SLM에 삽입
        pause(0.03);                                                            % SLM 및 카메라 속도 조절
        trigger(Camera)                                                         % Camera 찍기
        pause(str2num(Camera_parameter{1})/1000);
    end
    Shot_total = getdata(Camera,(Number_of_SLM_super_pixel_length^2)*4);
    %
    
%     u0=propTF(uu.*exp(-1i*au_),L1,lambda,2*z);                                                      % simulate SLM
%     [u0]=focus(u0,L1,lambda,zf);
%     u0=propTF(u0,L1,lambda,2*z);
    
    JJ(1,ii+1) = sum(sum(Shot_total));                                                            % calculate total value
    J(1,ii+1) = chk_J(Shot_total,S);                                                              % calculate target value
    Jratio(1,ii+1)=1/(JJ(1,ii+1)/J(1,ii+1))*100;                                                    % ratio of total to target
    W(1,ii+1)=((J(1,ii+1) - J(1,ii))/J(1,ii+1));                                                    % value of change in J
    half=mean(au-bu);                                                                               % get average of change amount of u
    weight=(W(1,ii+1)); %->multithread Jratio(1,ii+1)                                               % calcualte J'
    BB=(abs(au - bu)./(max(max(abs(au - bu)))/(2*pi)));                                             % (choose one of two) normallize of change of u
    BB=abs(rem(au,(2*pi)) - rem(bu,(2*pi)));                                                        % (choose one of two) change of actual phase of u
    WM = rem((BB.* weight), (2*pi));                                                                % change factor of au
    D(1,ii+1) = sum(sum(WM))/sum(sum(au))+1;                                                        % expectation of au's improvement
    MaxI(1,ii+1)=max(max(Shot_total));                                                            % find image's maximum value
    MinI(1,ii+1)=min(min(Shot_total));                                                            % find image's minimum value
    bu = au;                                                                                        % save au in bu
    au = (au + WM);                                                                                 % update au
    MU(1,ii+1)=max(max(au))/(2*pi);                                                                 % calculate au's cycle
    au = rem(au, (2*pi));                                                                           % simplify au
    if ii==1                                                                                        % after first iterration figure
        hold off
        I3=Shot_total;
        R3=radCal(I3);
%         figure(11) 
        subplot(4,4,5);
        imagesc(x1/1e-3,y1/1e-3,I3);
        hold on
        rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
        text(S/2*dx1/1e-3, -S/2*dx1/1e-3,['J=',num2str(Jratio(1,ii+1))],'Color','red');
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('result'); 
        colorbar
        hold off
    end
end

% P3=angle(u0)/pi+1;                                                                                  % calculate phase
I3=Shot_total;                                                                                    % extract image
R3=radCal(I3);                                                                                      % calculate radius

hold off                                                                                            % after SPGD image
% figure(5) 
subplot(4,4,7);
imagesc(x1/1e-3,y1/1e-3,I3);
hold on
rectangle('Position',[-S/2*dx1/1e-3 -S/2*dx1/1e-3 S*dx1/1e-3 S*dx1/1e-3], 'EdgeColor','r','LineWidth',1);
text(S/2*dx1/1e-3, -S/2*dx1/1e-3,['J=',num2str(Jratio(1,ii+1))],'Color','red');
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('result'); 
colorbar
hold off

% figure(6)                                                                                         % phase of after SPGD image
% subplot(4,4,8);
% imagesc(x1/1e-3,y1/1e-3,P3); 
% % xlim([-0.2 0.2]); ylim([-0.2 0.2]);
% axis square; axis xy; 
% colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
% title('align'); 
% colorbar

% figure(7)                                                                                         % pure image of slm control matrix
subplot(4,4,9);
imagesc(bu); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(8)                                                                                         % rescaled figure of J, delta J, change amount of au, cycle of au
subplot(4,4,10);
hold off
plot(rescale(J),'b*-');                                                                             % cacluated J value during iterration
hold on
plot(rescale(W),'go-');                                                                             % improvement of output during iterration
plot(rescale(D),'r^-');                                                                             % expectation of improvement of au
plot(rescale(MU),'y^-');                                                                            % au's cycle count
title('J'); 
legend('J','del J','u diff','u max');%
hold off

% figure(9)                                                                                         % output max/min value figure
subplot(4,4,11);
% subplot(2,1,1);
plot(MaxI,'b*-'); 
title('Max'); 
subplot(4,4,12);
% subplot(2,1,2);
plot(MinI,'ro-'); 
title('Min'); 

% figure(10)                                                                                        % tiff image figure
subplot(4,4,13);
returnVal = padarray(expand(round(bu*((2^16)/max(max(bu))))),[128 128],0,'both');                   % generate SLM modulate matrix
returnVal = uint16(returnVal);                                                                      % image matrix normalization
imagesc(returnVal); 
axis square; axis xy; 
colormap('gray');
colorbar
imwrite(returnVal,"ExportImage.tiff",'tiff');                                                       % tiff image export

% figure(12)                                                                                        % J Value figure
subplot(4,4,14);
hold off
plot(J,'b*-'); 
hold on
title('J');  xlabel('iteration count'); ylabel('J value'); 
hold off

% figure(13)                                                                                        % PBR figure
subplot(4,4,15);
hold off
plot(Jratio,'b*-'); 
hold on
title('PBR');  xlabel('iteration count'); ylabel('percentage'); 
hold off