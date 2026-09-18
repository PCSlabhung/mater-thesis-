%% 檢查不同 Tx 的「深度對齊」狀況 (A-scan Alignment Check)
clc; close all; clear all;

%% 1. 參數設置
folder_path = '0_0_150/'; 
c = 340; 
fs = 3.125e6/19;
range = [300:6:2000]/fs*c/2; % 深度向量
Nr = length(range);
no_lines = 69;             
no_lines_phi = 69;

% 要檢查的通道列表
check_list = [1, 2, 3, 4, 5]; 
% 必須補上第 5 個顏色 (這裡加一個 'k' 黑色給 Tx01)
colors = {'k', 'r', 'b', 'm', 'g'}; 

% 必須補上第 5 個線寬 (由粗到細排列)
line_widths = [6, 5, 4, 3, 2];
%% 2. 先讀取 Sum 檔案以鎖定目標角度
fname_sum = fullfile(folder_path, 'SA_image_sum.txt');
fprintf('正在讀取 Sum 檔案以鎖定目標位置...\n');

if isfile(fname_sum)
    raw_sum = dlmread(fname_sum);
    % 必須 Reshape 才能確保維度正確 [深度 x 方位 x 仰角]
    data_sum = reshape(raw_sum.', [Nr, no_lines, no_lines_phi]); 
    data_sum = abs(data_sum);
    
    % 找出最強點的索引 (Index)
    [max_val, max_idx] = max(data_sum(:));
    [idx_r, idx_az, idx_el] = ind2sub(size(data_sum), max_idx);
    
    fprintf('  -> 鎖定目標位置: 方位 Index %d, 仰角 Index %d\n', idx_az, idx_el);
    fprintf('  -> 峰值深度約: %.4f m (Index %d)\n', range(idx_r), idx_r);
else
    error('找不到 Sum 檔案，無法鎖定目標！');
end

%% 3. 繪圖迴圈
figure('Name', '深度對齊檢查 (Fixed Angle)');
hold on;

% 畫出各個 Tx 通道
for i = 1:length(check_list)
    tx_idx = check_list(i);
    fname = fullfile(folder_path, sprintf('SA_image_tx%02d.txt', tx_idx)); 
    
    if isfile(fname)
        raw = dlmread(fname);
        % Reshape
        data = reshape(raw.', [Nr, no_lines, no_lines_phi]);
        
        % --- 關鍵：取出與 Sum 相同的 (方位, 仰角) 的深度線 ---
        % 這樣我們比較的就是同一條視線上的深度訊號
        profile = squeeze(data(:, idx_az, idx_el));
        
        % 正規化 (Normalize)
        profile = abs(profile);
        profile = profile / max(profile); % 歸一化到 0~1
        profile_dB = 20*log10(profile + 1e-12);
        
        % 繪圖
        % 技巧：前面的 Tx 畫很粗，後面的 Tx 畫很細
        % 如果它們重疊，你會看到中間有一條細線，外圍有粗線的顏色
        plot(range, profile_dB, ...
             'Color', colors{i}, ...
             'LineWidth', line_widths(i), ...
             'DisplayName', sprintf('Tx %02d', tx_idx));
    else
        fprintf('警告: 找不到檔案 %s\n', fname);
    end
end

% 最後畫上 Sum (黑色虛線)
profile_sum = squeeze(data_sum(:, idx_az, idx_el));
profile_sum = profile_sum / max(profile_sum);
plot(range, 20*log10(profile_sum + 1e-12), 'k--', 'LineWidth', 2, 'DisplayName', 'SA Sum');

%% 4. 設定顯示範圍與標籤
grid on;
grid minor;
xlabel('深度 Range (m)');
ylabel('強度 (dB)');
title(['深度對齊檢查 (固定角度 Az:' num2str(idx_az) ', El:' num2str(idx_el) ')']);
legend('Location', 'best');

% 自動 Zoom In 到峰值附近 (前後 10 cm)
center_depth = range(idx_r);
xlim([center_depth - 0.1, center_depth + 0.1]); 
ylim([-30 1]); % 顯示前 30dB

hold off;

%% 5. (新增) 特別診斷: Tx01 vs Tx05 深度差異比較
% 目的: 確認 Tx01 是否有偏移，導致總和 (Sum) 的波束變寬或位置跑掉

figure('Name', '深度診斷: Tx01 vs Tx05');
hold on; grid on; grid minor;

% --- 讀取 Tx01 數據 ---
fname_tx01 = fullfile(folder_path, 'SA_image_tx01.txt');
if isfile(fname_tx01)
    raw01 = dlmread(fname_tx01);
    data01 = reshape(raw01.', [Nr, no_lines, no_lines_phi]);
    % 鎖定與 Sum 相同的角度
    prof01 = squeeze(data01(:, idx_az, idx_el));
    prof01 = abs(prof01) / max(abs(prof01)); % 正規化
    prof01_dB = 20*log10(prof01 + 1e-12);
    
    % 找出 Tx01 的峰值深度
    [~, loc01] = max(prof01);
    depth01 = range(loc01);
    
    plot(range, prof01_dB, 'r', 'LineWidth', 2.5, 'DisplayName', 'Tx 02 (嫌疑目標)');
else
    depth01 = NaN;
    fprintf('找不到 Tx01 檔案\n');
end

% --- 讀取 Tx05 數據 (作為標準參考) ---
fname_tx05 = fullfile(folder_path, 'SA_image_tx05.txt');
if isfile(fname_tx05)
    raw05 = dlmread(fname_tx05);
    data05 = reshape(raw05.', [Nr, no_lines, no_lines_phi]);
    % 鎖定與 Sum 相同的角度
    prof05 = squeeze(data05(:, idx_az, idx_el));
    prof05 = abs(prof05) / max(abs(prof05)); % 正規化
    prof05_dB = 20*log10(prof05 + 1e-12);
    
    % 找出 Tx05 的峰值深度
    [~, loc05] = max(prof05);
    depth05 = range(loc05);
    
    plot(range, prof05_dB, 'g--', 'LineWidth', 2.5, 'DisplayName', 'Tx 05 (對照組)');
else
    depth05 = NaN;
    fprintf('找不到 Tx05 檔案\n');
end

% --- 畫出 SA Sum 作為參考 ---
% (使用上一段程式已經讀好的 data_sum)
prof_sum = squeeze(data_sum(:, idx_az, idx_el));
prof_sum = abs(prof_sum)/max(abs(prof_sum));
plot(range, 20*log10(prof_sum+1e-12), 'k:', 'LineWidth', 1.5, 'DisplayName', 'SA Sum');

% --- 計算並顯示誤差 ---
if ~isnan(depth01) && ~isnan(depth05)
    diff_m = depth01 - depth05;
    diff_cm = diff_m * 100;
    
    title_str = sprintf('診斷結果: Tx01 與 Tx05 深度差 = %.4f m (%.2f cm)', diff_m, diff_cm);
    title(title_str, 'Color', 'r', 'FontSize', 12);
    fprintf('\n=== 深度誤差診斷 ===\n');
    fprintf('Tx01 峰值位置: %.5f m\n', depth01);
    fprintf('Tx05 峰值位置: %.5f m\n', depth05);
    fprintf('兩者誤差: %.5f m\n', diff_m);
    
    if abs(diff_m) < 1e-4
        fprintf('-> 判定: 深度完全對齊。\n');
    else
        fprintf('-> 判定: 深度不對齊！這可能是波束變寬的原因。\n');
    end
else
    title('無法計算誤差 (檔案缺失)');
end

% --- 圖表設定 ---
xlabel('深度 (m)');
ylabel('強度 (dB)');
legend('Location', 'best');

% 自動 Zoom In 到 Tx05 的峰值附近
if ~isnan(depth05)
    xlim([depth05 - 0.05, depth05 + 0.05]); % 前後 5 公分
    ylim([-20 1]); 
end

hold off;