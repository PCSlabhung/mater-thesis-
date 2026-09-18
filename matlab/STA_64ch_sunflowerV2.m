%%% Air US using 5Tx and multi-channel Rx
%%%  Revised on 8/11/2013 by Gency Jeng
 

clc
close all
clear all
% load SA_5data_V2
toff=0;
%%
% 設定5個Tx座標
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp=0.0783+15/1000;
right_x=temp;
left_x=-temp;
top_y=temp;
bot_y=-temp;
vs_x = [0  left_x right_x 0 0];
vs_z = 0;
vs_y=[0 0 0 top_y bot_y];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 設定模擬參數
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the transducer apertures for send and receive
f0 = 40*10^3;                  %  Transducer center frequency [Hz]
c = 340;                   %  Speed of sound [m/s]
focal_depth = 10^10*10^10; %FocalValue(1)*10^(-3);  
lambda = c/f0;
%%% Array geometry defined by Field II.
pitch = lambda/2;
fs = 3.125e6/19;   % sampling rate for raw data 
% [Nz M no_tx]=size(d);
Nz=2000;
M=64;
no_tx=5;

%%
% 設定RX座標
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=5;
Rap=0.059;
for m=1:64
    RR(m)=Rap*sqrt(m/M);
    if RR(m)<=15/1000
        RR(m)=18/1000;
    end
    phi(m)=2*pi*(m-1)*(1+sqrt(V))/2;
    chx(m)=RR(m)*cos(phi(m));
    chy(m)=RR(m)*sin(phi(m));   
end
chz=zeros(1,M);



%%
% scanline數量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_lines = 69;
no_lines_phi = 69;
d_th = 0.025574095887433;
d_th_phi = 0.025574095887433;

%%
%深度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Depth information for B-mode
% dz = 1/(f0)*c/2;
% range = [0.3:dz:1.5];  
range=[300:6:2000]/fs*c/2;
Rmin = range(1);
Rmax = range(end);
Nr = max(size(range));
Rx_enable = [1:1:M];

% write_range_LUT(range, 'range_LUT.h');
%% Show array and target in clearer 2D views

% Example target
target = [0, 0, 0.5];   % [x y z] in meters

figure;

%----------------------------
% Top view: X-Y plane
%----------------------------
subplot(1,2,1);
hold on; grid on; axis equal;

scatter(chx, chy, 36, 'b', 'filled', 'DisplayName', 'Rx elements');
scatter(vs_x, vs_y, 90, 'r', 'filled', 'DisplayName', 'Tx positions');
scatter(target(1), target(2), 120, 'k', '*', 'DisplayName', 'Target');

xlabel('X (m)');
ylabel('Y (m)');
title('Top View (X-Y Plane)');
legend('Location','best');

%----------------------------
% Side view: X-Z plane
%----------------------------
subplot(1,2,2);
hold on; grid on; axis equal;

scatter(chx, chz, 36, 'b', 'filled', 'DisplayName', 'Rx elements');
scatter(vs_x, vs_z*ones(size(vs_x)), 90, 'r', 'filled', 'DisplayName', 'Tx positions');
scatter(target(1), target(3), 120, 'k', '*', 'DisplayName', 'Target');

xlabel('X (m)');
ylabel('Z (m)');
title('Side View (X-Z Plane)');
legend('Location','best');
%%
% STA DAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bf = zeros(Nr,no_lines,no_lines_phi); 
apo_vector = ones(1,M);
all_subbf= zeros(Nr,no_lines,no_lines_phi);
all_CF= ones(Nr,no_lines,no_lines_phi);
CF= ones(Nr,no_lines,no_lines_phi);
subbf= zeros(Nr,no_lines,no_lines_phi);
subbf_paf= zeros(Nr,no_lines,no_lines_phi);
create_bb_dataV2
index_mode = 'hardware';

% ==============================
% Delay approximation error stats
% ==============================
cnt_err = 0;

sum_abs_dist_mm = 0;
sum_sq_dist_mm  = 0;
max_abs_dist_mm = 0;

sum_abs_delay_ns = 0;
sum_sq_delay_ns  = 0;
max_abs_delay_ns = 0;

sum_abs_sample = 0;
sum_sq_sample  = 0;
max_abs_sample = 0;

sum_abs_phase_deg = 0;
sum_sq_phase_deg  = 0;
max_abs_phase_deg = 0;

max_info.sample_value = -inf;
max_info.tx = NaN;
max_info.theta_idx = NaN;
max_info.phi_idx = NaN;
max_info.depth_idx = NaN;
max_info.rx_idx = NaN;
max_info.depth_m = NaN;
max_info.theta_deg = NaN;
max_info.phi_deg = NaN;
max_info.rx_ch = NaN;

for ii=1:1:no_tx
% for ii=1:1:1

    bb_data = squeeze(d(1:max(size(d,1)),:,ii)); % for each ii element

     % Demodulation+LPF
    bb_data = [bb_data ; zeros(1,M)];
    [N_rx, N] = size(bb_data);
    

        for ms = 1:no_lines_phi
                % sin_phii = -1*(no_lines_phi-1)/2*d_th_phi  + (ms-1)*d_th_phi;
                sin_phi =  -1*(no_lines_phi-1)/2*d_th_phi  + (ms-1)*d_th_phi;
        
                cos_phi = cos(asin(sin_phi));
            for jj = 1: no_lines  %% for each scan line
                 
                 sin_th = -1*(no_lines-1)/2*d_th  + (jj-1)*d_th;
                 cos_th = cos(asin(sin_th));
               
                 x_org=range.'*sin_th*cos_phi;
                 y_org=range.'*sin_phi;
                 z_org=range.'*cos_th.*cos_phi;
        
                 %Tx
                 Tx_dist = sqrt(   (x_org*ones(1,N)-vs_x(ii)).^2 + ...
                                    (z_org*ones(1,N)-vs_z).^2 + ...
                                    ( y_org*ones(1,N)-vs_y(ii)).^2    )  ;
                
                 %Rx
                 onewaydist = sqrt((x_org*ones(1,N)-ones(Nr,1)*chx(Rx_enable)).^2 + ...
                                    (z_org*ones(1,N)-ones(Nr,1)*chz(Rx_enable)).^2 + ...
                                    (y_org*ones(1,N)-ones(Nr,1)*chy(Rx_enable)).^2); 
                 T_rd = sqrt((range'.*ones(1,N)).^2 + (chx(Rx_enable).*ones(1,N)).^2 +(chy(Rx_enable).*ones(1,N)).^2);
                 T_sd_rd = -(chx(Rx_enable).*cos_phi*sin_th.*ones(1,N) + chy(Rx_enable).*sin_phi.*ones(1,N));
                 
                 totdelay_sw = (onewaydist + Tx_dist)/c - toff;
                 totdelay_hw = (Tx_dist + T_rd + T_sd_rd) / c;

                 
% ============================================================
% Accumulate error statistics without storing all error samples
% ============================================================

rx_exact  = onewaydist;
rx_approx = T_rd + T_sd_rd;

dist_error_mm = (rx_approx - rx_exact) * 1e3;   % mm
delay_error_ns = dist_error_mm * 1e-3 / c * 1e9; % ns
sample_error = delay_error_ns * 1e-9 * fs;       % samples
phase_error_deg = delay_error_ns * 1e-9 * f0 * 360; % degree

abs_dist_mm = abs(dist_error_mm);
abs_delay_ns = abs(delay_error_ns);
abs_sample = abs(sample_error);
abs_phase_deg = abs(phase_error_deg);

num_now = numel(dist_error_mm);
cnt_err = cnt_err + num_now;

sum_abs_dist_mm = sum_abs_dist_mm + sum(abs_dist_mm(:));
sum_sq_dist_mm  = sum_sq_dist_mm  + sum(dist_error_mm(:).^2);
max_abs_dist_mm = max(max_abs_dist_mm, max(abs_dist_mm(:)));

sum_abs_delay_ns = sum_abs_delay_ns + sum(abs_delay_ns(:));
sum_sq_delay_ns  = sum_sq_delay_ns  + sum(delay_error_ns(:).^2);
max_abs_delay_ns = max(max_abs_delay_ns, max(abs_delay_ns(:)));

sum_abs_sample = sum_abs_sample + sum(abs_sample(:));
sum_sq_sample  = sum_sq_sample  + sum(sample_error(:).^2);
max_abs_sample = max(max_abs_sample, max(abs_sample(:)));

sum_abs_phase_deg = sum_abs_phase_deg + sum(abs_phase_deg(:));
sum_sq_phase_deg  = sum_sq_phase_deg  + sum(phase_error_deg(:).^2);
max_abs_phase_deg = max(max_abs_phase_deg, max(abs_phase_deg(:)));

% Find local maximum sample error in this scanline
[local_max_sample, local_linear_idx] = max(abs_sample(:));

% If this local max is larger than global max, update global max info
if local_max_sample > max_info.sample_value

    [depth_idx, rx_local_idx] = ind2sub(size(abs_sample), local_linear_idx);

    max_info.sample_value = local_max_sample;
    max_info.tx = ii;
    max_info.theta_idx = jj;
    max_info.phi_idx = ms;
    max_info.depth_idx = depth_idx;
    max_info.rx_idx = rx_local_idx;

    max_info.depth_m = range(depth_idx);

    max_info.theta_deg = asin(sin_th) * 180/pi;
    max_info.phi_deg   = asin(sin_phi) * 180/pi;

    rx_list = find(Rx_enable);
    max_info.rx_ch = rx_list(rx_local_idx);
end

max_abs_sample = max(max_abs_sample, local_max_sample);
                 switch(index_mode)
                    case'software'
                        totdelay = (onewaydist + Tx_dist)/c - toff;
                    case'hardware'
                        totdelay = (Tx_dist + T_rd + T_sd_rd) / c; 
                end
                 %totdelay = (onewaydist + Tx_dist)/c-toff;
                 delayindx = (totdelay*fs)+1;
                 flagdelay = (delayindx < 1) | (delayindx > N_rx);
                 delayindx = delayindx.*(1-flagdelay) + N_rx*flagdelay;         
                 delayindx_low = floor(delayindx) + N_rx*ones(Nr,1)*(0:N-1);
                 delay_out = find(delayindx_low >= (N_rx*N));
                 delayindx_low(delay_out) = N_rx*N-1; 
                 interp_bbdata = bb_data(delayindx_low);
                 chandata = (ones(Nr,1)*apo_vector).*interp_bbdata.*exp(j*2*pi*f0*totdelay);
                 subbf(:,jj,ms) = sum((ones(Nr,1)*apo_vector).*interp_bbdata.*exp(j*2*pi*f0*totdelay),2) ;
        
                 CF(:,jj,ms)=( (abs(subbf(:,jj,ms))).^2 )  ./  (   M*  sum(abs (chandata.^2) , 2 ) )+0.00000001  ;
        
                 alph =0;  
                 subbf_paf(:, jj, ms) = (1 - alph) * subbf(:, jj, ms) + alph * (subbf(:, jj, ms) .^2);
                 

            end
        end


    all_subbf=all_subbf+subbf_paf.*CF.^2;

end

mean_dist_mm = sum_abs_dist_mm / cnt_err;
rmse_dist_mm = sqrt(sum_sq_dist_mm / cnt_err);

mean_delay_ns = sum_abs_delay_ns / cnt_err;
rmse_delay_ns = sqrt(sum_sq_delay_ns / cnt_err);

mean_sample = sum_abs_sample / cnt_err;
rmse_sample = sqrt(sum_sq_sample / cnt_err);

mean_phase_deg = sum_abs_phase_deg / cnt_err;
rmse_phase_deg = sqrt(sum_sq_phase_deg / cnt_err);

fprintf('\n========== RX Delay Approximation Error ==========\n');

fprintf('Distance error : mean = %.6f mm, RMSE = %.6f mm, max = %.6f mm\n', ...
    mean_dist_mm, rmse_dist_mm, max_abs_dist_mm);

fprintf('Delay error    : mean = %.6f ns, RMSE = %.6f ns, max = %.6f ns\n', ...
    mean_delay_ns, rmse_delay_ns, max_abs_delay_ns);

fprintf('Sample error   : mean = %.6f samples, RMSE = %.6f samples, max = %.6f samples\n', ...
    mean_sample, rmse_sample, max_abs_sample);

fprintf('Phase error    : mean = %.6f deg, RMSE = %.6f deg, max = %.6f deg\n', ...
    mean_phase_deg, rmse_phase_deg, max_abs_phase_deg);

fprintf('==================================================\n');

fprintf('\n========== Maximum Sample Error Location ==========\n');
fprintf('Max sample error = %.6f samples\n', max_info.sample_value);
fprintf('TX index         = %d\n', max_info.tx);
fprintf('Theta index      = %d\n', max_info.theta_idx);
fprintf('Phi index        = %d\n', max_info.phi_idx);
fprintf('Theta angle      = %.3f deg\n', max_info.theta_deg);
fprintf('Phi angle        = %.3f deg\n', max_info.phi_deg);
fprintf('Depth index      = %d\n', max_info.depth_idx);
fprintf('Depth            = %.6f m\n', max_info.depth_m);
fprintf('RX local index   = %d\n', max_info.rx_idx);
fprintf('RX channel index = %d\n', max_info.rx_ch);
fprintf('===================================================\n');

%% Plot the image
DR = 40;
bfn_data = abs(all_subbf/max(all_subbf(:))); %normalize
bfn_data(isnan(bfn_data))=1e-20;

qqq = 20*log10(bfn_data); %取log
sin_th = -1*(no_lines-1)/2*d_th + ([1:1:no_lines]-1)*d_th; 
cos_th = cos(asin(sin_th));
sin_phi =  -1*(no_lines_phi-1)/2*d_th + ([1:1:no_lines_phi]-1)*d_th; 
cos_phi = cos(asin(sin_phi));


       
Z=range.'*cos_th;
X=range.'*sin_th;

QQ=((qqq(:,:,(no_lines_phi+1)/2))) ;
%%
%XZ平面
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
pcolor(X, Z, QQ)
axis('ij');
shading interp;
caxis([-DR 0])
colormap hot
xlabel('azimuth [m]')
ylabel('range [m]')
colorbar


%%
%XY平面
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



range_index=find(qqq==0);
cmode_range=mod(range_index,size(range,2));
cmode_range=cmode_range(1); 

XX=range(cmode_range)*sin_th'*cos_phi;
YY=range(cmode_range)*ones(no_lines,1)*sin_phi;
temp=12;
r_min=cmode_range-temp;
r_max=cmode_range+temp;
if(r_min<1)
    r_min=1;
end
if(r_max>size(range,2))
    r_max=size(range,2);
end
 c_qqq=squeeze(max(qqq(r_min:3:r_max,:,:),[],1));
[x, y] = meshgrid(1:no_lines, 1:no_lines_phi);
[xq, yq] = meshgrid(linspace(1, no_lines, 1000), linspace(1, no_lines_phi, 1000));
A_interp = interp2(x, y, c_qqq, xq, yq, 'linear'); 
XX_interp = interp2(x, y, XX, xq, yq, 'linear'); 
YY_interp = interp2(x, y, YY, xq, yq, 'linear'); 

temp=linspace(-60,60,no_lines);
[XX_polar, YY_polar] = meshgrid(temp,temp);
XX_polar = interp2(x, y, XX_polar, xq, yq, 'linear'); 
YY_polar = interp2(x, y, YY_polar, xq, yq, 'linear'); 
A_interp_temp=A_interp;
% A_interp_temp(650:850,650:850)=-100;%30 30
A_interp_temp(350:650,350:650)=-100;%0 0
A_interp_temp(:,1:100)=-100;
max_idx=find(A_interp==max(A_interp_temp(:)));
max_point=zeros(1,1000);
max_point_raw=floor(max_idx/1000)+1;
if(max_point_raw==1001)
    max_point_raw=1000;
end
max_point_col=max_idx-floor(max_idx/1000)*1000;
if max_point_col==0
    max_point_col=1000
end

figure
pcolor(XX_polar,YY_polar,A_interp')
hold on
caxis([-DR 0])
shading interp;
axis equal
colormap("hot")
ttemp=zeros(1,1000);
ttemp(max_point_col)=YY_polar(max_point_raw,1);
ttemp(ttemp==0)=NaN;
% plot(fliplr(XX_polar(1,:)),fliplr(ttemp),'ro')
plot((XX_polar(1,:)),(ttemp),'yo')

max_sidelobe_peak=A_interp(max_idx)
check_SSL=max(A_interp_temp(:))
filename = sprintf('10 project polar line=%d sidelobe=%.2f', no_lines,max_sidelobe_peak);

% title(filename)
xlabel('theta[°]')
ylabel('phi[°]')
aa=gca;aa.FontSize=24
IMAG=zeros(Nr,100,100);




%%
%beampattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bp=c_qqq(:,(no_lines+1)/2);
% bp=c_qqq(:,38);%實際
% bp=c_qqq(:,37);%模擬

inter_no_lines=ones(1,no_lines);
inter_no_lines=interp1([1: length(inter_no_lines)],inter_no_lines,linspace(1, length(inter_no_lines), length(inter_no_lines)*100),'spline');
bp3=interp1([1: length(bp)],bp,linspace(1, length(bp), length(bp)*100),'spline');
range3=interp1([1: length(sin_th)],range(cmode_range)*sin_th,linspace(1, length(sin_th), ...
length(sin_th)*100),'linear');

figure
str = sprintf('%d', no_lines)
plot(range3,bp3,'DisplayName',str )
hold on 
plot(range3,-6*inter_no_lines)

L=find(bp3>-6);
res_=range3(L(end))-range3(L(1))
aa=gca;aa.FontSize=24;
ylabel('dB')
xlabel('X(cm)')
title('beampattern@50cm')

%%
%3D imag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    x_range = 3.5;
    y_range = 3.5;
    z_range = 2;

    x_pixel=256;
    y_pixel=256;
    z_pixel=100;

    pos_vec_x_new = (0:1 / (x_pixel-1):1) .* x_range - x_range / 2;
    pos_vec_y_new = (0:1 / (y_pixel-1):1) .* y_range - y_range / 2;
    pos_vec_z_new = linspace(0.3,1.5,z_pixel);
    rr=range;
    [pos_mat_z_new,pos_mat_y_new, pos_mat_x_new] = ndgrid(pos_vec_z_new,pos_vec_x_new, pos_vec_y_new);
    [R,THETA, PHI] = ndgrid(rr,asin(sin_th), asin(sin_phi));
    % convert new points to polar coordinates
    [th_cart, r_cart,z_cart] = cart2pol(pos_mat_z_new,pos_mat_x_new, pos_mat_y_new);
    phi_cart=atan2(z_cart,r_cart);
    r_cart=r_cart./cos(phi_cart);

    figure   
    bb_mode = interp3( THETA,R,PHI, qqq,phi_cart, r_cart, th_cart, 'linear');
    % bb_mode(isnan(bb_mode))=0;
    bb_mode(bb_mode<=-DR)=NaN;
    aa =slice(pos_vec_y_new,  pos_vec_z_new,pos_vec_x_new, bb_mode,[],pos_vec_z_new, []);
    set(aa,'EdgeColor','none');
    caxis([-DR 0])
    colormap("hot")
     axis equal
     xlabel('X (m)')
     ylabel('Z (m)')
     zlabel('Y (m)')
    view([20 10])
     title('C scan');
    % axis equal
     pos2 = get(gca, 'Position'); % 獲取當前子圖位置
set(gca, 'Position', [pos2(1)-0.05, pos2(2), pos2(3), pos2(4)]); % 調整大小和位置
      % xlim([-1 1])
      % zlim([-1 1])
      % xlim([-2 2])
      % zlim([-0.5 0.5])
     ylim([0 1.9])
     % colormap('jet')
function write_range_LUT(range, filename)
    % If the user forgets ".h", automatically add it
    if ~endsWith(filename, '.h')
        filename = [filename '.h'];
    end

    N = length(range);

    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot create file: %s', filename);
    end

    % Header guard name
    [~, name] = fileparts(filename);
    guard = upper([name '_H']);
    guard = regexprep(guard, '\W', '_'); % replace non-word chars

    fprintf(fid, '#ifndef %s\n', guard);
    fprintf(fid, '#define %s\n', guard);

    fprintf(fid, 'const TX_t range[%d] = {', N);

    % Write each value with full precision
    for i = 1:N
        fprintf(fid, '%.16f', range(i));
        if i ~= N
            fprintf(fid, ',');
        end
    end

    fprintf(fid, '};\n');
    fprintf(fid, '#endif\n');

    fclose(fid);
end