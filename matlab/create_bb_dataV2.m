%%% Air US using 1Tx and multi-channel Rx
%%%  Revised on 8/11/2013 by Gency Jeng
 path(path,'C:\DSP計畫\DSP計畫\3D\polar\0318\0318\fieldII');
 
field_init
clc
% close all
% clear all

%---------------------------------------------------------------------
%%
% 設定5個Tx座標

Rap=0.118/2
temp=0.0783+15/1000;
right_x=temp;
left_x=-temp;
top_y=temp;
bot_y=-temp;

vs_x = [0  left_x right_x   0 0];
vs_z = 0;
vs_y=[0 0 0 top_y bot_y];

%---------------------------------------------------------------------
%%
%設定模擬參數
f0 = 40*10^3;              %center of frequency    
c = 340;                   %  Speed of sound [m/s]
focal_depth = 10^10*10^10; %FocalValue(1)*10^(-3);  
lambda = c/f0;
pitch = lambda/2;
width = 1.85/1000;            %pitch - kerf;                   
kerf = pitch-width;          %  Kerf [m]
height = 2.75/1000;            %  height of element
focus=[0 0 focal_depth];     %  Fixed focal point [m]
Rfocus = 10^20/1000;%Rconvex; %%%%%35mm
fs = 3.125e6/19;   % sampling rate for raw data 
% fs = 3.125e6/16;   % sampling rate for raw data 
fs_simu = fs*4; % simulation frequency
set_sampling(fs_simu)
set_field('c',c);

%---------------------------------------------------------------------
%%
%設定Tx探頭格式丟給模擬用

R=12.15/1000;%TX半徑
ele_size=1/1000;
xmit_aperture = xdc_piston (R,ele_size);
data=xdc_get(xmit_aperture,'rect');
[M,N]=size(data);
rect1=zeros(5*N,19);
rect1(:,1)=1;
rect1(:,14)=1;
rect1(:,15)=1/1000;
rect1(:,16)=1/1000;
for ii=1:1:N
    rect1(ii,1)=1;

    rect1(ii,2)=data(11,ii);
    rect1(ii,3)=data(12,ii);
    rect1(ii,4)=data(13,ii);
    rect1(ii,5)=data(14,ii);
    rect1(ii,6)=data(15,ii);
    rect1(ii,7)=data(16,ii);
    rect1(ii,8)=data(17,ii);
    rect1(ii,9)=data(18,ii);
    rect1(ii,10)=data(19,ii);
    rect1(ii,11)=data(20,ii);
    rect1(ii,12)=data(21,ii);
    rect1(ii,13)=data(22,ii);
    
    rect1(ii,17)=data(24,ii);
    rect1(ii,18)=data(25,ii);
    rect1(ii,19)=data(26,ii);
end
    iii=1;


for ii=N+1:1:2*N
    rect1(ii,1)=2;
    rect1(ii,2)=data(11,iii)+left_x;
    rect1(ii,3)=data(12,iii);
    rect1(ii,4)=data(13,iii);
    rect1(ii,5)=data(14,iii)+left_x;
    rect1(ii,6)=data(15,iii);
    rect1(ii,7)=data(16,iii);
    rect1(ii,8)=data(17,iii)+left_x;
    rect1(ii,9)=data(18,iii);
    rect1(ii,10)=data(19,iii);
    rect1(ii,11)=data(20,iii)+left_x;
    rect1(ii,12)=data(21,iii);
    rect1(ii,13)=data(22,iii);
    
    rect1(ii,17)=data(24,iii)+left_x;
    rect1(ii,18)=data(25,iii);
    rect1(ii,19)=data(26,iii);
    iii=iii+1;
end


iii=1;


for ii=2*N+1:1:3*N
    rect1(ii,1)=3;
    rect1(ii,2)=data(11,iii)+right_x;
    rect1(ii,3)=data(12,iii);
    rect1(ii,4)=data(13,iii);
    rect1(ii,5)=data(14,iii)+right_x;
    rect1(ii,6)=data(15,iii);
    rect1(ii,7)=data(16,iii);
    rect1(ii,8)=data(17,iii)+right_x;
    rect1(ii,9)=data(18,iii);
    rect1(ii,10)=data(19,iii);
    rect1(ii,11)=data(20,iii)+right_x;
    rect1(ii,12)=data(21,iii);
    rect1(ii,13)=data(22,iii);
    
    rect1(ii,17)=data(24,iii)+right_x;
    rect1(ii,18)=data(25,iii);
    rect1(ii,19)=data(26,iii);
    iii=iii+1;
end
iii=1;


for ii=3*N+1:1:4*N
    rect1(ii,1)=4;
    rect1(ii,2)=data(11,iii);
    rect1(ii,3)=data(12,iii)+top_y;
    rect1(ii,4)=data(13,iii);
    rect1(ii,5)=data(14,iii);
    rect1(ii,6)=data(15,iii)+top_y;
    rect1(ii,7)=data(16,iii);
    rect1(ii,8)=data(17,iii);
    rect1(ii,9)=data(18,iii)+top_y;
    rect1(ii,10)=data(19,iii);
    rect1(ii,11)=data(20,iii);
    rect1(ii,12)=data(21,iii)+top_y;
    rect1(ii,13)=data(22,iii);
    
    rect1(ii,17)=data(24,iii);
    rect1(ii,18)=data(25,iii)+top_y;
    rect1(ii,19)=data(26,iii);
    iii=iii+1;
end

iii=1;
for ii=4*N+1:1:5*N
    rect1(ii,1)=5;
    rect1(ii,2)=data(11,iii);
    rect1(ii,3)=data(12,iii)+bot_y;
    rect1(ii,4)=data(13,iii);
    rect1(ii,5)=data(14,iii);
    rect1(ii,6)=data(15,iii)+bot_y;
    rect1(ii,7)=data(16,iii);
    rect1(ii,8)=data(17,iii);
    rect1(ii,9)=data(18,iii)+bot_y;
    rect1(ii,10)=data(19,iii);
    rect1(ii,11)=data(20,iii);
    rect1(ii,12)=data(21,iii)+bot_y;
    rect1(ii,13)=data(22,iii);
    
    rect1(ii,17)=data(24,iii);
    rect1(ii,18)=data(25,iii)+bot_y;
    rect1(ii,19)=data(26,iii);
    iii=iii+1;
end


center=[0 vs_y(1) 0];
xmit_aperture = xdc_rectangles (rect1, ones(5,1)*center, focus);
%---------------------------------------------------------------------
%%
%設定RX座標
N_enabled=64
M=N_enabled;

jj=5

for m=1:64
    RR(m)=Rap*sqrt(m/64);
    if RR(m)<=15/1000
        RR(m)=18/1000;
    end
    phii(m)=2*pi*(m-1)*(1+sqrt(5))/2;
    chx(m)=RR(m)*cos(phii(m));
    chy(m)=RR(m)*sin(phii(m));
    chz=zeros(1,M);

    
end
%---------------------------------------------------------------------
%%
%設定RX格式丟給模擬用

rx=zeros(M,19);
for ii=1:M
rx(ii,:) =            [ii,chx(ii)-width/2,chy(ii)-width/2,chz(ii), ...
                          chx(ii)-width/2,chy(ii)+width/2,chz(ii), ...
                          chx(ii)+width/2,chy(ii)+width/2,chz(ii), ...
                          chx(ii)+width/2,chy(ii)-width/2,chz(ii), ...
                          1,width,height,chx(ii),chy(ii),chz(ii)];
end
receive_aperture= xdc_rectangles(rx,[0,0,0],focus);
%---------------------------------------------------------------------
%%
%設定發射波型
excitation=1*cos(2*pi*f0*(0:1/fs_simu:1/f0));
pwm_t = 0:1/fs_simu:50e-6;
pwm1 = double(mod(pwm_t, 2/80e3) < 1/80e3); % 第1組PWM：從 0 µs 開始
pwm1 = [ pwm1 ];
excitation_fft=abs(fftshift(fft(pwm1)));
% figure;
% plot( pwm1, 'b', 'LineWidth', 1.5); hold on;
% figure
% plot([-floor(length(excitation))/2 : floor(length(excitation))/2-1]*(fs_simu/length(excitation)),abs(fftshift(fft(excitation))))
excitation=pwm1;
%---------------------------------------------------------------------
%%
%設定TX探頭響應
impulse_response_TX=cos(2*pi*f0*(0:1/fs_simu:10/f0));
N=length(impulse_response_TX);
NN=N-1;
f=[-NN/2:NN/2]/NN*fs_simu;
% figure
% FFT_TX=abs(fftshift(fft(impulse_response_TX)));
% FFT_TX=FFT_TX/max(FFT_TX);
% FFT_TX_DB=20*log10(FFT_TX);
% ff=interp(f,100);
% plot(f,FFT_TX_DB)
% hold on
% plot(ff,-6*ones(1,(length(ff))))
%---------------------------------------------------------------------
%%
%設定RX探頭響應
impulse_response_RX=cos(2*pi*28000*(0:1/fs_simu:2.5/f0));
% N=length(impulse_response_RX);
% if mod(N,2)==0
%     NN=N;
%     f=[-NN/2:NN/2-1]/NN*fs_simu;
% 
% else
%     NN=N-1;
%     f=[-NN/2:NN/2]/NN*fs_simu;
% 
% end
% NN=N-1;
% figure
% FFT_RX=abs(fftshift(fft(impulse_response_RX)));
% FFT_RX=FFT_RX/max(FFT_RX);
% FFT_RX_DB=20*log10(FFT_RX);
% ff=interp(f,100);
% 
% plot(f,FFT_RX_DB)
% hold on
% 
% plot(ff,-6*ones(1,(length(ff))))
%%
% figure
% plot(impulse_response_TX)
xdc_impulse (xmit_aperture, impulse_response_TX);
xdc_impulse (receive_aperture, impulse_response_RX);
xdc_excitation (xmit_aperture, excitation);
%---------------------------------------------------------------------
%%
%設定scatterers座標/強度

scatterers=[
                0    0       1.0   ;% FOR 1ball 50cm
                % -0.05        0.02       0.37   ;% FOR 1ball 50cm

           ];
 

amp = ones(size(scatterers,1),1);

%---------------------------------------------------------------------
%%
%模擬幾何目標

% x_size=0.02; %寬度
% y_size=0.02;%高度
% z_size=0.02;%厚度
% z_start=0.5;%離探頭距離
% 
% pht_x_num=x_size/(2*lambda);
% pht_y_num=y_size/(2*lambda);
% pht_z_num=z_size/(2*lambda);
% N=round(10*(pht_x_num)*(pht_y_num)*(pht_z_num));
% aver=zeros(Nr,no_lines,no_lines_phi);
% aver_idx=1;
% 
% x = (rand (N,1)-0.5)*x_size;
% y = (rand (N,1)-0.5)*y_size;
% z = (rand (N,1)-0.5)*z_size;
% amp = randn(N,1)*1;
% r=x_size/2;      %  Radius of cyst [mm]
% xc= 0;    %  Place of cyst [mm]
% yc= 0;
% % zc=z_start+r;  %50
% %%%目標為圓柱體 
%  inside = ( ((x).^2 + (z).^2  ) < r^2);
%  amp = amp .* (inside);
% amp(amp==0)=NaN;
% 
% %%%
% xx=[ x;];
% zz=[z+z_start; ];
% yy=[y];
% num=size(xx,1);
% amp=(randn(num,1));
% amp=[amp;];
% scatterers=[xx yy zz ];
% figure
%      % scatter3(xx,zz,yy,5,amp,'filled')
%      % hold on
%      scatter3(chx,chz,chy,5,'filled')
%      hold on
%      scatter3(vs_x,vs_z,vs_y,10,'blue','filled')
%      amp(isnan(amp))=min(amp);
%      colormap('hot')
%      alpha(1)
%      axis equal
%      view([45 22.5])
%      aa=gca;aa.FontSize=24
%      xlabel('X (m)')
%      ylabel('Z (m)')
%      zlabel('Y (m)')

%---------------------------------------------------------------------
%%
%設定非聚焦和無加權
xdc_center_focus (receive_aperture, [0 0 0]);
xdc_times_focus(receive_aperture, 0,zeros(1,N_enabled));
apo_rx = ones(1, N_enabled);
xdc_apodization(receive_aperture, 0, apo_rx);       
apo_tx = ones(1, 5);
xdc_apodization(xmit_aperture, 0, apo_tx);        
%---------------------------------------------------------------------
%%
%生成raw data
[rrf_data, toff] = calc_scat_all(xmit_aperture, receive_aperture, scatterers, amp,1);

%---------------------------------------------------------------------
%%
% downsampling
rrf_data = rrf_data(1:floor(fs_simu/fs):end,:,:);
AA=zeros(2000,5*N_enabled);
toff_idx=round(toff*fs);
AA(toff_idx:toff_idx+size(rrf_data,1)-1,:)=rrf_data;
rrf_data=AA;
rrf_data=rrf_data/max(rrf_data(:));
toff=0;
%---------------------------------------------------------------------


Rx_enable = [1:1:N_enabled];
len_Rx = length(Rx_enable);
[Nz, Nx] = size (rrf_data);
Nx=N_enabled;
%%% Demodulation filter
lpf=[1 3 5 6 5 3 1];

tt = [0:Nz-1]/fs;
tt = tt(ones(Nx,1),:);
data_SL=reshape(rrf_data,[size(rrf_data,1) M 5]);      


% Demodulation
for ii=1:1:5
data_ = squeeze(data_SL(1:max(size(data_SL,1)),:,ii)); 
bb_data_temp = conv2(1,lpf,(data_.').*exp(-j*2*pi*f0.*tt),'same').';
d(:,:,ii) = [bb_data_temp ];

end



