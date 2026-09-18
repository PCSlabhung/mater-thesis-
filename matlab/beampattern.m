%% 疊加比較 5個 C++ 合成孔徑輸出的 Beampattern (含 Sum) - 公平比較版
% 修改重點：計算寬度時，每一條線使用自己的 Peak - 6dB 作為閾值
clc;
close all;
clear all;

%% 1. 設置參數
% 檔案名稱 (請確認資料夾名稱是否正確)
filename1 = '30_30_100/SA_image_tx01.txt';   
filename2 = '30_30_100/SA_image_tx02.txt';
filename3 = '30_30_100/SA_image_tx03.txt';
filename4 = '30_30_100/SA_image_tx04.txt';
filename5 = '30_30_100/SA_image_tx05.txt';
filename_sum = '30_30_100/SA_image_sum.txt';    

% 模擬常數
c = 340;                   
fs = 3.125e6/19;           
no_lines = 69;             
no_lines_phi = 69;         
d_th = 0.025574095887433;  

% 深度 (Range) 資訊
range = [300:6:2000]/fs*c/2; 
Nr = max(size(range));     
DR = 40;                   

% 計算角度向量
sin_th = -1 * (no_lines - 1) / 2 * d_th + ([1:1:no_lines] - 1) * d_th; 

%% 2. 定義讀取與鎖定目標的函數
get_bp_data = @(fname) extract_beampattern(fname, Nr, no_lines, no_lines_phi, range, sin_th, DR);

function [bp_interp, x_axis_interp, peak_info] = extract_beampattern(fname, Nr, no_lines, no_lines_phi_target, range_vec, sin_th_vec, DR)
    % 讀取檔案
    try
        % raw = dlmread(fname);
        raw = sqrt(max(dlmread(fname), 0));
    catch
        error('無法讀取 %s，請確認檔案路徑是否正確', fname);
    end
    
    [total_rows, total_cols] = size(raw);
    
    % 1. 檢查深度
    if total_cols ~= Nr
        warning('[%s] 警告: 檔案列數 (Depth) 為 %d，預期為 %d。', fname, total_cols, Nr);
        current_Nr = total_cols;
        range_vec = range_vec(1:min(length(range_vec), current_Nr));
    else
        current_Nr = Nr;
    end
    
    % 2. 檢查角度線數量
    expected_rows = no_lines * no_lines_phi_target;
    if total_rows < expected_rows
        valid_phi_lines = floor(total_rows / no_lines);
        if valid_phi_lines < 1
            error('檔案 %s 數據過少，無法構成單一完整掃描面。', fname);
        end
        raw = raw(1 : valid_phi_lines * no_lines, :);
        current_no_lines_phi = valid_phi_lines;
    else
        current_no_lines_phi = no_lines_phi_target;
        raw = raw(1 : expected_rows, :);
    end
    
    % 3. Reshape
    bf = reshape(raw.', [current_Nr, no_lines, current_no_lines_phi]);
    
    % 正規化與轉 dB
    bfn = abs(bf / max(abs(bf(:))));
    bfn(isnan(bfn)) = 1e-20;
    qqq = 20 * log10(bfn);
    % qqq = 10 * log10(bfn);
    % --- 自動鎖定最強點 ---
    [max_val, max_idx] = max(qqq(:));
    [idx_r, idx_az, idx_el] = ind2sub(size(qqq), max_idx);
    
    % --- [MIP] 提取波束 ---
    r_margin = 5;
    r_start = max(1, idx_r - r_margin);
    r_stop = min(current_Nr, idx_r + r_margin);
    
    slab = qqq(r_start:r_stop, :, idx_el);
    bp_raw = max(slab, [], 1); 
    bp_raw = bp_raw(:);       
    
    % 計算對應的物理 X 軸
    target_depth = range_vec(idx_r);
    x_axis_raw = target_depth * sin_th_vec;
    
    % --- 高解析度插值 ---
    x_axis_interp = linspace(min(x_axis_raw), max(x_axis_raw), length(x_axis_raw)*100);
    bp_interp = interp1(x_axis_raw, bp_raw, x_axis_interp, 'spline');
    
    % --- 清洗 NaN ---
    bp_interp(isnan(bp_interp)) = -DR - 10;
    
    % 紀錄峰值資訊
    peak_info.depth = target_depth;
    peak_info.elev_idx = idx_el;
    peak_info.max_dB = max(bp_interp);
end

%% 3. 處理數據
fprintf('正在處理數據 Tx01 (%s)...\n', filename1);
[bp1, x1, info1] = get_bp_data(filename1);
fprintf('正在處理數據 Tx02 (%s)...\n', filename2);
[bp2, x2, info2] = get_bp_data(filename2);
fprintf('正在處理數據 Tx03 (%s)...\n', filename3);
[bp3, x3, info3] = get_bp_data(filename3);
fprintf('正在處理數據 Tx04 (%s)...\n', filename4);
[bp4, x4, info4] = get_bp_data(filename4);
fprintf('正在處理數據 Tx05 (%s)...\n', filename5);
[bp5, x5, info5] = get_bp_data(filename5);
fprintf('正在處理數據 Sum (%s)...\n', filename_sum);
[bp_sum, x_sum, info_sum] = get_bp_data(filename_sum);

%% 4. 繪製疊加圖 (原始強度比較)
figure('Name', 'Beampattern 絕對強度比較');

% 繪製線條
plot(x1, bp1, 'LineWidth', 2, 'DisplayName', 'Tx01'); hold on;
plot(x2, bp2, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Tx02');
plot(x3, bp3, 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', 'Tx03');
plot(x4, bp4, 'LineWidth', 2, 'LineStyle', '-.', 'DisplayName', 'Tx04');
plot(x5, bp5, 'LineWidth', 1, 'DisplayName', 'Tx05'); 
plot(x_sum, bp_sum, 'LineWidth', 3, 'LineStyle', '-', 'Color', 'k', 'DisplayName', 'Sum');

grid on;
ylabel('強度 (dB)');
xlabel('方位角距離 Azimuth X (m)');
title(sprintf('絕對強度 Beampattern @ Depth %.2f m', info_sum.depth));
legend('Location', 'best');

% 設定 Y 軸範圍
peak_ref = max(bp_sum);
y_upper = peak_ref + 5;
y_lower = max(peak_ref - 50, -DR - 5);
if y_lower >= y_upper, y_lower = y_upper - 40; end
ylim([y_lower, y_upper]);
hold off;

%% 5. 公平比較：計算各自的波束寬度 (動態閾值)
all_bp = {bp1, bp2, bp3, bp4, bp5, bp_sum};
all_x  = {x1, x2, x3, x4, x5, x_sum};
names  = {'Tx01', 'Tx02', 'Tx03', 'Tx04', 'Tx05', 'Sum'};
widths = zeros(6,1);
peaks  = zeros(6,1);

fprintf('\n=== 公平比較：基於各自峰值下降 6dB 的寬度 ===\n');

for i = 1:6
    current_bp = all_bp{i};
    current_x  = all_x{i};
    
    % 1. 找出該線條自己的峰值
    my_peak = max(current_bp);
    peaks(i) = my_peak;
    
    % 2. 設定專屬於該線條的閾值 (Peak - 6dB)
    my_threshold = my_peak - 6;
    %my_threshold = my_peak - 3;
    % 3. 計算寬度
    widths(i) = calculate_width_robust(current_x, current_bp, my_threshold);
    
    fprintf('%s \tPeak: %6.2f dB \tWidth: %.4f m\n', names{i}, my_peak, widths(i));
end

%% 6. 繪製正規化視圖 (Normalized View) - 讓大家都在 0dB 比較
figure('Name', 'Beampattern 正規化比較 (形狀對比)');

hold on;
styles = {'-', '--', ':', '-.', '-', '-'};
colors = lines(6); colors(6,:) = [0 0 0]; % Sum 為黑色
widths_plot = [1.5, 1.5, 1.5, 1.5, 1.5, 2.5];

for i = 1:6
    % 將每個波束的最大值平移至 0 dB
    bp_norm = all_bp{i} - peaks(i); 
    plot(all_x{i}, bp_norm, 'LineStyle', styles{i}, 'Color', colors(i,:), ...
         'LineWidth', widths_plot(i), 'DisplayName', names{i});
end

% 畫統一的 -6 dB 參考線
yline(-6, 'r--', 'LineWidth', 1.5, 'DisplayName', '-6 dB Reference');
%yline(-3, 'r--', 'LineWidth', 1.5, 'DisplayName', '-3 dB Reference');
grid on;
ylabel('正規化強度 (dB)');
xlabel('方位角距離 Azimuth X (m)');
title('正規化 Beampattern (形狀與寬度公平比較)');
legend('Location', 'best');
ylim([-30 2]); % 專注看 Mainlobe 形狀
hold off;

% 在正規化圖上顯示計算結果
str_info = cell(6,1);
for i=1:6
    str_info{i} = sprintf('%s Width: %.4f m', names{i}, widths(i));
end
annotation('textbox', [0.15, 0.15, 0.2, 0.2], 'String', str_info, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FaceAlpha', 0.8);

%% 8. 輔助函數
function width = calculate_width_robust(x_axis, bp_data, threshold)
    % 找出大於閾值的區域
    above_threshold = bp_data >= threshold;
    
    % 找出最大值位置
    [~, max_idx] = max(bp_data);
    
    % 向左搜尋邊界
    left_idx = max_idx;
    while left_idx > 1 && above_threshold(left_idx)
        left_idx = left_idx - 1;
    end
    
    % 向右搜尋邊界
    right_idx = max_idx;
    while right_idx < length(bp_data) && above_threshold(right_idx)
        right_idx = right_idx + 1;
    end
    
    % 線性插值求更精確的交點 X
    if right_idx > left_idx && left_idx > 1 && right_idx < length(bp_data)
        % 左邊插值 (y = threshold)
        x1 = x_axis(left_idx); y1 = bp_data(left_idx);
        x2 = x_axis(left_idx+1); y2 = bp_data(left_idx+1);
        x_left_exact = x1 + (threshold - y1) * (x2 - x1) / (y2 - y1);
        
        % 右邊插值
        x3 = x_axis(right_idx-1); y3 = bp_data(right_idx-1);
        x4 = x_axis(right_idx); y4 = bp_data(right_idx);
        x_right_exact = x3 + (threshold - y3) * (x4 - x3) / (y4 - y3);
        
        width = x_right_exact - x_left_exact;
    elseif right_idx > left_idx
        % 如果無法插值 (例如邊界)，直接用索引差
        width = x_axis(right_idx) - x_axis(left_idx);
    else
        width = NaN; 
    end
end