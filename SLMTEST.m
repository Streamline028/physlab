clear all 
close all
clc

addpath(genpath('.\Book_Reference_Code'))
addpath(genpath('.\module'))

% addpath(genpath('C:\Users\a0105\Documents\MATLAB\TM'))

%% SLM calibration
BNS_OpenSLMs()                                                             % SLM 함수 열기
BNS_LoadLUTFile(1, ... 
    'C:\Users\a0105\Documents\MATLAB\TM\LUT_Files\slm5240_at1064_P16.lut') % SLM calibration LUT 함수(1064 nm용) 삽입
Image.optimization_data = libpointer('uint8Ptr', zeros(524288,1));         % SLM image 삽입 포인터 지정 (Calibration용)
BNS_ReadTIFF('C:\PCIe16MatlabSDK\Image_Files\White.tiff', ...
    Image.optimization_data);                                              % Calibration 용 white 이미지 삽입 (2 pi phase (2^16) 전체 삽입)
BNS_WriteCal(Image.optimization_data, Image);                              % Calibration
%% SLM on
BNS_SetPower(true)                                                         % SLM 전원 on
Image.Image0 = libpointer('uint8Ptr', zeros(524288,1));                    % SLM 이미지 데이터 삽입용 포인터 지정
name_initial = [...
'C:\Users\a0105\Documents\MATLAB\TM\Tilt_basis_200_waves_DR_16_bits.tiff'];       
% SLM 이미지 파일 사전 지정 -> 초기 Blank 이미지 삽입
BNS_ReadTIFF(name_initial,Image.Image0);                                   % 이미지 포인터에 삽입
BNS_WriteImage(Image.Image0, Image);                                       % Calibration
%% Fiber input carmera setup
Camera = TM_camera_renewal;                                                % image acquisition toolbox에서 Export한 m-file
Camera_size = Camera.ROIPosition;
prompt = {'x_offset: ', 'y_offset: ',['cam_width (max:', ... 
    num2str(Camera_size(3)),'): '],['cam_height (max:', ... 
    num2str(Camera_size(4)),'): '],'# of screen super pixel: ','Exposure'};
dlg_title = 'Screenshot parameter';
num_lines = 1;
def = {'768','768','512','512','16','-2.5'};                               % Default setting
% # of screen super pixel 16일 때 120
% # of screen super pixel 32일 때 496 497 528 529
% # of screen super pixel 64일 때 1951 1952 1953 1954 2015 2016 2017 2018
% 2079 2080 2081 2082 2143 2144 2145 2146
screenshot_parameter = inputdlg(prompt,dlg_title,num_lines,def);

x_offset=str2num(screenshot_parameter{1});                                 % 카메라 x 시작점
y_offset=str2num(screenshot_parameter{2});                                 % 카메라 y 시작점
cam_width=str2num(screenshot_parameter{3});                                % 카메라로 보는 사이즈 길이 Minimum 32
cam_height=str2num(screenshot_parameter{4});                               % 카메라로 보는 사이즈 높이 Minimum 2
number_of_screen_super_pixels=str2num(screenshot_parameter{5});            % 카메라 화면 나눈 갯수(정사각형, 한 변 갯수)
Source.Exposure = str2num(screenshot_parameter{6});
% Source.Shutter = str2num(screenshot_parameter{6});
screen_super_pixel_width=fix(cam_width/number_of_screen_super_pixels);     % 카메라 화면 나눈 조각 크기(길이, 정사각형)
Source = Camera.Source;
Source.ShutterMode = 'auto';
Source.GainMode = 'auto';
Source.Exposure = -2.5;
preview(Source.Parent)
%% initial setting
load('initset.mat');                                                       % get some matrix for testing

Iterration_Count=1000;                                                     % iterration count
Checking_Size=screen_super_pixel_width*2;                                  % searching size
Super_Pixels=16;                                                           % SLM superpixel count

M=512;                                                                     % screen size (pixel)
L1=M*15e-6;                                                                % screen size (m)
dx1=L1/M;                                                                  % pixel interval (m)
x1=-L1/2:dx1:L1/2-dx1;                                                     % generate x field
y1=x1;                                                                     % generate y field
lambda=1.064e-6;                                                           % laser wave length (m)
k=2*pi/lambda;                                                             % wave number
w=dx1*(160/2);                                                             % laser radius

z=180e-3;                                                                  % first lens distance
zf=z;                                                                      % first lens's focal length
z2=8e-3;                                                                   % second lens distance
zf2=z2;                                                                    % second lens's focal length

%% Add Tilt
Tilt_reference_mode = ones(512,512);
% Set ramp numbers
% Caution: Considering the maximum number of waves in SLM size. 
Number_of_wave = 200;
% Set your SLM dynamic range as bit-number.
for i = 1:1:512
    Tilt_reference_mode(:,i) = ((2^16-1)/(512/Number_of_wave))*(i-1);
    Tilt_reference_mode = mod(Tilt_reference_mode,2^16-1);
end

%% before SPGD

tic;
Shot_initial = getsnapshot(Source.Parent);                                 % 초기 speckle image
Shot_merge_initial = resize(picavg(Shot_initial,...
    screen_super_pixel_width), screen_super_pixel_width);                  % super pixel 단위 카메라 output 사전 지정

Before_U = zeros(Super_Pixels,Super_Pixels).*(2*pi);                       % generate initial before u (since then bu)
After_U = rand([Super_Pixels,Super_Pixels]).*(2*pi);                       % generate initial after u (since then au)
% au = preset6;

Before_U_expand_ori = padarray(expand(Before_U),[128 128],0,'both');       % generate SLM mask
Before_U_expand = Tilt_reference_mode+(Before_U_expand_ori*2^16-1);
Before_U_expand = mod(Before_U_expand,2^16-1);
Before_U_expand = uint16(Before_U_expand);

Image.Image0 = SM_tiffencoder(Before_U_expand);

BNS_WriteImage(Image.Image0, Image);                                       % 만들어진 Hadamard basis matrix를 SLM에 삽입
pause(0.03);                                                               % SLM 및 카메라 속도 조절
BNS_WriteImage(Image.Image0, Image);                                       % 만들어진 Hadamard basis matrix를 SLM에 삽입
pause(0.03);                                                               % SLM 및 카메라 속도 조절
Shot_total_init = getsnapshot(Source.Parent);     
% Camera 찍기
pause(3*Source.Shutter/1000);

Before_Beam_Intensity=Shot_total_init;                                     % extract image
Before_Beam_Radius=radCal(Before_Beam_Intensity);                          % calculate radius

Image_Intensity_Sum(1,1) = sum(sum(Before_Beam_Intensity));                % calculate total value
Target_Intensity_Sum(1,1) = chk_J(Before_Beam_Intensity,Checking_Size);    % calculate target value
Target_Intensity_Ratio(1,1)=1/(Image_Intensity_Sum(1,1)/...
    Target_Intensity_Sum(1,1))*100;                                        % ratio of total to target
dJ(1,1)= 1;                                                                % value of change in J - initially 0
MaxValue_U(1,1)=max(max(angle(Before_U)))/(2*pi);                          % calculate bu's cycle
Max_Intensity(1,1)=max(max((abs(Before_Beam_Intensity).^2)));              % find initial image's maximum value
Min_Intensity(1,1)=min(min((abs(Before_Beam_Intensity).^2)));              % find initial image's minimum value

hold off                                                                   % before SPGD image
% figure(3)
subplot(4,4,3);
imagesc(x1/1e-3,y1/1e-3,Before_Beam_Intensity);
hold on
rectangle('Position',...
    [-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 ...
    Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3],...
    'EdgeColor','r','LineWidth',1);
text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,['J=',...
    num2str(Target_Intensity_Ratio(1,1))],'Color','red');
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title(['result : ','radius = ',num2str(Before_Beam_Radius*0.0150),'mm']); 
colorbar
hold off

%% after u expand
After_U_expand_ori = padarray(expand(After_U),[128 128],0,'both');         % generate SLM mask
%     bu_ori = padarray(expand(bu),[128 128],0,'both');                    % generate SLM mask
After_U_expand = Tilt_reference_mode+(After_U_expand_ori*2^16-1);
After_U_expand = mod(After_U_expand,2^16-1);
After_U_expand = uint16(After_U_expand);
Image.Image0 = SM_tiffencoder(After_U_expand);
    BNS_WriteImage(Image.Image0, Image);                                   % 만들어진 Hadamard basis matrix를 SLM에 삽입
pause(0.03);                                                               % SLM 및 카메라 속도 조절
BNS_WriteImage(Image.Image0, Image);                                       % 만들어진 Hadamard basis matrix를 SLM에 삽입
pause(0.03);  

%% SPGD Loop

for ii = 1:Iterration_Count                                                % run until changement is 0
    
    After_U_expand_ori = padarray(expand(After_U),[128 128],0,'both');     % generate SLM mask
%     bu_ori = padarray(expand(bu),[128 128],0,'both');                    % generate SLM mask
    After_U_expand = Tilt_reference_mode+(After_U_expand_ori*2^16-1);
    After_U_expand = mod(After_U_expand,2^16-1);
    After_U_expand = uint16(After_U_expand);
    Image.Image0 = SM_tiffencoder(After_U_expand);
    
    BNS_WriteImage(Image.Image0, Image);                                   % 만들어진 Hadamard basis matrix를 SLM에 삽입
    pause(0.03);                                                           % SLM 및 카메라 속도 조절
    BNS_WriteImage(Image.Image0, Image);                                   % 만들어진 Hadamard basis matrix를 SLM에 삽입
    pause(0.03);                                                           % SLM 및 카메라 속도 조절
    Shot_total = getsnapshot(Source.Parent);                               % Camera 찍기
    pause(1);
    
    Image_Intensity_Sum(1,ii+1) = sum(sum(Shot_total));                    % calculate total value
    Target_Intensity_Sum(1,ii+1) = chk_J(Shot_total,Checking_Size);        % calculate target value
    Target_Intensity_Ratio(1,ii+1)=...
        Target_Intensity_Sum(1,ii+1)/Image_Intensity_Sum(1,ii+1)*100;      % ratio of total to target
    dJ(1,ii+1)=(Target_Intensity_Sum(1,ii+1) - Target_Intensity_Sum(1,ii));% value of change in J
    Gamma(1,ii) = ...
        (max(Target_Intensity_Sum)*pi)/(Target_Intensity_Sum(1,ii+1));     % calculate Gamma value
    weight=Gamma(1,ii).*(dJ(1,ii+1));                                      % calcualte Gamma X dJ
    dU=(After_U-Before_U);                                                 % change of actual phase of u
    J_prime = rem((BB.* weight), (2*pi));                                  % change factor of au
    Max_Intensity(1,ii+1)=max(max(Shot_total));                            % find image's maximum value
    Min_Intensity(1,ii+1)=min(min(Shot_total));                            % find image's minimum value
    Before_U = After_U;                                                    % save au in bu
    After_U = (After_U + J_prime);                                         % update au
    MaxValue_U(1,ii+1)=max(max(After_U))/(2*pi);                           % calculate au's cycle
    After_U = rem(After_U, (2*pi));                                        % simplify au
    figure(51)
    subplot(4,1,1);
    imagesc(Before_U)
    axis image
    title ('Phase before')
    xticklabels ([])
    yticklabels ([])
    colorbar
    subplot(4,1,2);
    imagesc(After_U)
    axis image
    title (['Phase after ',num2str(ii)])
    xticklabels ([])
    yticklabels ([])
    colorbar
    subplot(4,1,3);
    imagesc(BB)
    axis image
    title ('Phase difference')
    xticklabels ([])
    yticklabels ([])
    colorbar
    subplot(4,1,4);
    imagesc(J_prime)
    axis image
    title ('Phase changing')
    xticklabels ([])
    yticklabels ([])
    colorbar
    if ii==1                                                               % after first iterration figure
        hold off
        I3=Shot_total;
        R3=radCal(I3);
%         figure(11) 
        subplot(4,4,5);
        imagesc(x1/1e-3,y1/1e-3,I3);
        hold on
        rectangle('Position',...
            [-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 ...
            Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], ...
            'EdgeColor','r','LineWidth',1);
        text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,...
            ['J=',num2str(Target_Intensity_Ratio(1,ii+1))],'Color','red');
        axis square; axis xy; 
        colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
        title('result'); 
        colorbar
        hold off
    end
end
toc;
% P3=angle(u0)/pi+1;                                                       % calculate phase
I3=Shot_total;                                                             % extract image
R3=radCal(I3);                                                             % calculate radius

hold off                                                                   % after SPGD image
% figure(5) 
subplot(4,4,7);
imagesc(x1/1e-3,y1/1e-3,I3);
hold on
rectangle('Position',...
    [-Checking_Size/2*dx1/1e-3 -Checking_Size/2*dx1/1e-3 ...
    Checking_Size*dx1/1e-3 Checking_Size*dx1/1e-3], ...
    'EdgeColor','r','LineWidth',1);
text(Checking_Size/2*dx1/1e-3, -Checking_Size/2*dx1/1e-3,...
    ['J=',num2str(Target_Intensity_Ratio(1,ii+1))],'Color','red');
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray'); xlabel('x (mm)'); ylabel('y (mm)'); 
title('result'); 
colorbar
hold off

% figure(7)                                                                % pure image of slm control matrix
subplot(4,4,9);
imagesc(Before_U); 
% xlim([-0.2 0.2]); ylim([-0.2 0.2]);
axis square; axis xy; 
colormap('gray');
title('U'); 
colorbar

% figure(8)                                                                % rescaled figure of J, delta J, change amount of au, cycle of au
subplot(4,4,10);
hold off
plot(rescale(Target_Intensity_Sum),'b*-');                                 % cacluated J value during iterration
hold on
plot(rescale(dJ),'go-');                                                   % improvement of output during iterration
plot(rescale(MaxValue_U),'y^-');                                           % au's cycle count
title('J'); 
legend('J','del J','u diff','u max');%
hold off

% figure(9)                                                                % output max/min value figure
subplot(4,4,11);
% subplot(2,1,1);
plot(Max_Intensity,'b*-'); 
title('Max'); 
subplot(4,4,12);
% subplot(2,1,2);
plot(Min_Intensity,'ro-'); 
title('Min'); 

% figure(10)                                                               % tiff image figure
subplot(4,4,13);
returnVal = padarray(expand(...
    round(Before_U*((2^16)/max(max(Before_U))))),[128 128],0,'both');      % generate SLM modulate matrix
returnVal = uint16(returnVal);                                             % image matrix normalization
imagesc(returnVal); 
axis square; axis xy; 
colormap('gray');
colorbar
imwrite(returnVal,"ExportImage.tiff",'tiff');                              % tiff image export

% figure(12)                                                               % J Value figure
subplot(4,4,14);
hold off
plot(Target_Intensity_Sum,'b*-'); 
hold on
title('J');  xlabel('iteration count'); ylabel('J value'); 
hold off

% figure(13)                                                               % PBR figure
subplot(4,4,15);
hold off
plot(Target_Intensity_Ratio,'b*-'); 
hold on
title('PBR');  xlabel('iteration count'); ylabel('percentage'); 
hold off

%%
BNS_SetPower(false)
BNS_CloseSLM()
endtime = toc;
minute = fix(endtime/60);   sec = endtime-60*minute;
fprintf("Total time : %dm %.2fs\n",minute,sec)