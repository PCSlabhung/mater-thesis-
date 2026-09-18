%% 批量診斷：SA Sum 與 所有單次發射 (Tx) 的波束寬度比較
clc; close all; clear all;

%% 1. 設定基礎參數
folder_path = '30_30_150_float/';       % 資料夾路徑
sum_filename = 'SA_image_sum.txt';
c = 340; fs = 3.125e6/19;
no_lines = 69; no_lines_phi = 69; d_th = 0.025574095887433;
range = [300:6:2000]/fs*c/2;
Nr = length(range);
sin_th = -1 * (no_lines - 1) / 2 * d_th + ([1:1:no_lines] - 1) * d_th;

% 準備計算寬度的匿名函數
calc_width = @(x, y, th) (x(find(y>th,1,'last')) - x(find(y>th,1,'first')));

%% 2. 先處理 Sum 檔案 (作為基準)
fprintf('正在處理 Sum 檔案...\n');
full_path_sum = fullfile(folder_path, sum_filename);
[sum_width, sum_peak_depth] = analyze_beam_width(full_path_sum, Nr, no_lines, no_lines_phi, range, sin_th);
fprintf('>>> SA Sum -6dB 寬度: %.4f m (於深度 %.2f m)\n', sum_width, sum_peak_depth);

%% 3. 批量處理所有 Tx 檔案
% 假設檔案命名格式為 SA_image_tx0.txt, SA_image_tx1.txt ... 或 tx01.txt
% 這裡嘗試搜尋所有 tx 開頭的 txt
files = dir(fullfile(folder_path, 'SA_image_tx*.txt'));

tx_indices = [];
tx_widths = [];

fprintf('正在分析單次發射檔案 (共 %d 個)...\n', length(files));

for i = 1:length(files)
    fname = files(i).name;
    % 嘗試從檔名提取編號 (例如 tx05 -> 5)
    tokens = regexp(fname, 'tx(\d+)', 'tokens');
    if ~isempty(tokens)
        idx = str2double(tokens{1}{1});
        
        % 分析該檔案
        [w, d] = analyze_beam_width(fullfile(folder_path, fname), Nr, no_lines, no_lines_phi, range, sin_th);
        
        % 只有當偵測到的深度跟 Sum 差不多時才納入比較 (避免鎖定到雜訊)
        if abs(d - sum_peak_depth) < 0.2 % 容許 0.2m 的深度誤差
            tx_indices(end+1) = idx;
            tx_widths(end+1) = w;
            fprintf('  Tx %d: 寬度 %.4f m (深度 %.2f m)\n', idx, w, d);
        else
            fprintf('  Tx %d: 略過 (峰值深度 %.2f m 與 Sum 不符)\n', idx, d);
        end
    end
end

%% 4. 繪圖比較
figure('Name', '波束寬度診斷');
% 畫出所有單次發射的寬度
plot(tx_indices, tx_widths, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Single Tx Widths');
hold on;

% 畫出 SA Sum 的寬度 (紅線)
yline(sum_width, 'r--', 'LineWidth', 2, 'DisplayName', 'SA Sum Width');

xlabel('Tx Index (發射編號)');
ylabel('-6 dB Beam Width (m)');
title(sprintf('為什麼 SA 沒比較窄？ (目標深度 %.2f m)', sum_peak_depth));
legend('Location', 'best');
grid on;

% 統計分析
avg_tx_width = mean(tx_widths);
min_tx_width = min(tx_widths);
fprintf('\n=== 統計結果 ===\n');
fprintf('SA Sum 寬度: %.4f m\n', sum_width);
fprintf('Single Tx 平均寬度: %.4f m\n', avg_tx_width);
fprintf('Single Tx 最小寬度: %.4f m\n', min_tx_width);

if sum_width > min_tx_width
    fprintf('結論: SA Sum 比最窄的單次發射還要寬 -> 存在散焦或相位誤差。\n');
else
    fprintf('結論: SA Sum 是最窄的 -> 合成效果良好。\n');
end

%% --- 輔助函數 ---
function [width, peak_depth] = analyze_beam_width(filepath, Nr, no_lines, no_lines_phi, range_vec, sin_th_vec)
    try
        raw = dlmread(filepath);
        bf = reshape(raw.', [Nr, no_lines, no_lines_phi]);
        bfn = abs(bf / max(abs(bf(:))));
        bfn(isnan(bfn)) = 1e-20;
        qqq = 20 * log10(bfn);
        
        % 自動鎖定最強點
        [~, max_idx] = max(qqq(:));
        [idx_r, idx_az, idx_el] = ind2sub(size(qqq), max_idx);
        
        peak_depth = range_vec(idx_r);
        
        % 提取波束
        bp_raw = squeeze(qqq(idx_r, :, idx_el));
        x_axis = peak_depth * sin_th_vec;
        
        % 插值
        x_interp = linspace(min(x_axis), max(x_axis), length(x_axis)*100);
        bp_interp = interp1(x_axis, bp_raw, x_interp, 'spline');
        
        % 計算寬度
        threshold = -6;
        above = bp_interp >= (max(bp_interp) + threshold); % 相對峰值 -6dB
        if any(above)
            first = x_interp(find(above, 1, 'first'));
            last = x_interp(find(above, 1, 'last'));
            width = last - first;
        else
            width = NaN;
        end
    catch
        width = NaN; peak_depth = NaN;
    end
end
