%% STA_64ch_sunflowerV2 完整自動化分析程式 (最終修正版)
% 功能：自動偵測目標、繪製 B-mode, C-mode, Beampattern, 3D Slice
% 修正：MIP 投影、防崩潰機制、3D 維度校正

clc;
close all;
clear all;

%% 1. 設置參數
% 檔案名稱 (請確認路徑正確)
filename = 'no_CF_0_0_150/SA_image_sum.txt'; 

% 模擬常數
c = 340;                   % 聲速 [m/s]
fs = 3.125e6/19;           % 採樣率 

% 掃描線參數
no_lines = 69;             % 方位角掃描線數量
no_lines_phi = 69;         % 仰角掃描線數量
d_th = 0.025574095887433;  % 角度間隔 (方位角)
d_th_phi = 0.025574095887433; % 角度間隔 (仰角)

% 深度 (Range) 資訊
range = [300:6:2000]/fs*c/2; % 距離向量
Nr = max(size(range));     % 深度點數量 (NR)
Rmax = range(end);

% 繪圖參數
DR = 40;                   % 動態範圍 [dB]

%% 2. 載入數據並前處理
fprintf('正在載入數據: %s ...\n', filename);

% 讀取 TXT 檔案
try
    loaded_data = dlmread(filename);
catch ME
    error('無法讀取檔案 "%s"。請確認檔案是否存在且路徑正確。', filename);
end

% 檢查維度
expected_rows = no_lines * no_lines_phi;
if size(loaded_data, 1) ~= expected_rows || size(loaded_data, 2) ~= Nr
    error('維度錯誤。預期 %d 行 x %d 列，實際 %d x %d。', ...
          expected_rows, Nr, size(loaded_data, 1), size(loaded_data, 2));
end

% 重塑為 3D 陣列: (Nr, no_lines, no_lines_phi)
% 原始數據排列假設為 (Line_Total x Depth)，需轉置後 Reshape
bf = reshape(loaded_data.', [Nr, no_lines, no_lines_phi]);

% 正規化與 Log Compression
bfn_data = abs(bf / max(abs(bf(:)))); 
bfn_data(isnan(bfn_data)) = 1e-20;   
qqq = 20 * log10(bfn_data);          

%% 3. 全域變數準備與自動目標鎖定
% 建立角度向量
sin_th = -1 * (no_lines - 1) / 2 * d_th + ([1:1:no_lines] - 1) * d_th;
sin_phi = -1 * (no_lines_phi - 1) / 2 * d_th_phi + ([1:1:no_lines_phi] - 1) * d_th_phi;

% --- 自動尋找全域最強點 ---
[max_val, max_idx_3d] = max(qqq(:));
[idx_r, idx_az, idx_el] = ind2sub(size(qqq), max_idx_3d);

% 換算物理座標
target_range = range(idx_r);
target_sin_th = sin_th(idx_az);
target_sin_phi = sin_phi(idx_el);

% 計算笛卡爾座標 (用於 3D 繪圖)
target_x = target_range * target_sin_th;
target_y = target_range * target_sin_phi;
% 深度 Z 使用投影校正，若為虛數則直接用 Range
target_z_depth = target_range * sqrt(1 - target_sin_th^2 - target_sin_phi^2);
if ~isreal(target_z_depth), target_z_depth = target_range; end

fprintf('\n=== 自動偵測目標位置 ===\n');
fprintf('  Range Index: %d (%.4f m)\n', idx_r, target_range);
fprintf('  Azimuth Index: %d (%.2f deg)\n', idx_az, asin(target_sin_th)*180/pi);
fprintf('  Elevation Index: %d (%.2f deg)\n', idx_el, asin(target_sin_phi)*180/pi);
fprintf('  物理座標 (X, Y, Z): (%.2f, %.2f, %.2f)\n', target_x, target_y, target_z_depth);
fprintf('========================\n\n');

%% 4. XZ 平面 B-mode 圖像 (自動切片)
disp('繪製 XZ 平面 B-mode 圖...');

% 準備網格
cos_th = cos(asin(sin_th));
Z_bmode = range.' * cos_th;
X_bmode = range.' * sin_th;

% 提取目標仰角層
QQ = squeeze(qqq(:, :, idx_el));

figure('Name', 'XZ 平面 B-mode 圖像');
pcolor(X_bmode, Z_bmode, QQ);
shading interp;
axis('ij'); 
axis equal; axis tight;
caxis([-DR 0]);
colormap hot;
colorbar;

xlabel('方位 Azimuth X [m]');
ylabel('深度 Depth Z [m]');
title(sprintf('XZ B-mode @ Elev Index %d (%.1f^o)', idx_el, asin(target_sin_phi)*180/pi));

% 標示目標
hold on;
plot(X_bmode(idx_r, idx_az), Z_bmode(idx_r, idx_az), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

%% 5. XY 平面 C-mode 圖像 (自動切片)
disp('繪製 XY 平面 C-mode 圖...');

% 取目標深度前後 10 層做 MIP (最大強度投影)
thickness = 10;
r_min = max(1, idx_r - thickness);
r_max = min(Nr, idx_r + thickness);
c_data_raw = squeeze(max(qqq(r_min:r_max, :, :), [], 1)); % (Azimuth x Elevation)

% 準備網格 (注意 meshgrid 順序)
% X軸: sin_th (Azimuth), Y軸: sin_phi (Elevation)
[Grid_X_sin, Grid_Y_sin] = meshgrid(sin_th, sin_phi);

% 轉換為物理米
XX_meter = target_range * Grid_X_sin;
YY_meter = target_range * Grid_Y_sin;

% 為了配合 pcolor, 需要將 c_data_raw 轉置為 (Elevation x Azimuth)
c_data_plot = c_data_raw.'; 

% 插值平滑
scale = 4; % 解析度倍數
[Gx_fine, Gy_fine] = meshgrid( ...
    linspace(min(sin_th), max(sin_th), length(sin_th)*scale), ...
    linspace(min(sin_phi), max(sin_phi), length(sin_phi)*scale));
C_fine = interp2(sin_th, sin_phi, c_data_plot, Gx_fine, Gy_fine, 'spline');
XX_fine = target_range * Gx_fine;
YY_fine = target_range * Gy_fine;

figure('Name', 'XY 平面 C-mode 圖像');
pcolor(XX_fine, YY_fine, C_fine);
shading interp;
axis equal; axis tight;
colormap("hot");
caxis([-DR 0]);
colorbar;

xlabel('方位 Azimuth X (m)');
ylabel('仰角 Elevation Y (m)');
title(sprintf('C-mode Top View @ Range %.2f m', target_range));

% 標示目標
hold on;
plot(target_x, target_y, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

%% 6. 一維波束模式 (MIP 強健版)
disp('繪製 Beampattern (MIP mode)...');

% --- 關鍵修正：使用 MIP 避免切到空訊號 ---
r_margin_bp = 5; 
r_start_bp = max(1, idx_r - r_margin_bp);
r_stop_bp = min(Nr, idx_r + r_margin_bp);

% 取出小塊數據 (Depth_Slab, All_Azimuth, Target_Elev)
slab = qqq(r_start_bp:r_stop_bp, :, idx_el);
% 沿深度取最大值 -> 確保抓到波束主瓣
bp_raw = max(slab, [], 1); 
bp_raw = bp_raw(:);

% X 軸 (物理寬度)
x_axis_raw = target_range * sin_th(:);

% 插值
x_axis_interp = linspace(min(x_axis_raw), max(x_axis_raw), length(x_axis_raw)*100);
bp_interp = interp1(x_axis_raw, bp_raw, x_axis_interp, 'spline');

% 清洗 NaN
bp_interp(isnan(bp_interp)) = -DR - 10;

figure('Name', '一維波束模式');
plot(x_axis_interp, bp_interp, 'LineWidth', 2, 'DisplayName', 'Beampattern');
hold on;

% 找峰值與繪製 -6dB 線
peak_dB = max(bp_interp);
if peak_dB <= -DR, peak_dB = -DR+5; end % 防呆
threshold_dB = peak_dB - 6;

plot([min(x_axis_interp) max(x_axis_interp)], [threshold_dB threshold_dB], ...
     'r--', 'LineWidth', 1.5, 'DisplayName', '-6 dB Threshold');

% 計算寬度
above_th = bp_interp >= threshold_dB;
[~, loc_max] = max(bp_interp);
left = loc_max; while left>1 && above_th(left), left=left-1; end
right = loc_max; while right<length(bp_interp) && above_th(right), right=right+1; end

if right > left
    width = x_axis_interp(right) - x_axis_interp(left);
    text(x_axis_interp(loc_max), threshold_dB+2, sprintf('Width: %.4f m', width), ...
        'Color','r', 'FontWeight','bold', 'HorizontalAlignment','center');
    fprintf('Beampattern -6dB 寬度: %.4f m\n', width);
else
    fprintf('無法計算寬度 (未低於 -6dB)\n');
end

% --- 安全的 ylim 設定 ---
y_up = peak_dB + 5;
y_low = max(peak_dB - 50, -DR - 5);
if y_low >= y_up, y_low = y_up - 40; end
ylim([y_low, y_up]);

xlabel('方位 Azimuth X (m)');
ylabel('強度 (dB)');
title(sprintf('Beampattern (MIP) @ %.2f m', target_range));
grid on; legend('Location','best');
hold off;

%% 7. 3D 體積切片 (自動聚焦 + Meshgrid 版)
disp('繪製 3D 體積切片 (Zoom In)...');

% 設定聚焦範圍 (目標前後 20 cm)
zoom = 0.2; 
vec_azim = linspace(target_x - zoom, target_x + zoom, 80);
vec_depth = linspace(target_z_depth - zoom, target_z_depth + zoom, 80);
vec_elev = linspace(target_y - zoom, target_y + zoom, 80);

% 生成網格 (meshgrid 確保 slice 兼容性)
% 順序: [Xgrid, Ygrid, Zgrid] = meshgrid(x, y, z)
% 這裡對應物理意義: [Azim, Depth, Elev]
[Azim_grid, Depth_grid, Elev_grid] = meshgrid(vec_azim, vec_depth, vec_elev);

% 逆向轉換求查詢點 (R, sin_th, sin_phi)
R_q = sqrt(Azim_grid.^2 + Depth_grid.^2 + Elev_grid.^2);
Sin_th_q = Azim_grid ./ R_q;
Sin_phi_q = Elev_grid ./ R_q;

% 3D 插值
% qqq 維度: (Range, Azim, Elev)
% interp3 輸入順序: (X=Azim, Y=Range, Z=Elev)
bb_mode = interp3(sin_th, range, sin_phi, qqq, ...
                  Sin_th_q, R_q, Sin_phi_q, 'linear', -100);

% 過濾背景 (設為 NaN 使其透明，或 -100 顯示黑色)
bb_mode(bb_mode <= -DR) = NaN; 

figure('Name', '3D Volume Slice');
% slice 的輸入順序直接對應 meshgrid 的輸出
h_slice = slice(Azim_grid, Depth_grid, Elev_grid, bb_mode, ...
                target_x, target_z_depth, target_y);

set(h_slice, 'EdgeColor', 'none', 'FaceColor', 'interp');
caxis([-DR 0]);
colormap("hot");
colorbar;

xlabel('Azimuth X (m)');
ylabel('Depth Z (m)'); % 這裡顯示為 Z (深度)
zlabel('Elevation Y (m)');

axis equal; axis tight;
grid on; box on; view(3);
title(sprintf('3D Slice Locked @ (%.2f, %.2f, %.2f)', target_x, target_z_depth, target_y));

fprintf('所有繪圖完成。\n');