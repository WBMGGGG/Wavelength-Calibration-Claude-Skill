%% ==============================================================
% Wavelength Calibration - Generate Excel
% Function: Read ENVI peak results, 3rd-order calibration, formatted Excel with chart
% Dependencies: Windows + Excel + MATLAB R2020b+
% Input: offset value (interactive), pixel data from ENVI
% Output: YYYYMMDD.xlsx on desktop
%% ==============================================================

clear; clc; close all;

%% ===================== 1. Interactive Input =====================
% offset value from user input, no hard-code
offset_str = input('Enter offset integer value: ', 's');
offsetValue = str2double(offset_str);
while isnan(offsetValue) || mod(offsetValue, 1) ~= 0
    offset_str = input('Invalid input. Enter integer offset value: ', 's');
    offsetValue = str2double(offset_str);
end

% ENVI peak result sample format (replace with actual ENVI output)
% dir: sub-folder name; pixels: peak array (NaN = no peak found)
pixelData = [
    struct('dir','HG-1bin', 'pixels',[25.0, 105.0, 377.0]), ...
    struct('dir','NE-1bin', 'pixels',[493.0, 570.0, NaN, 594.0, 734.0, 799.0, 844.0]), ...
    struct('dir','AR-1bin', 'pixels',[696.0, 750.0, 763.5, 772.0, 794.5, 800.0, 811.5, 826.0, 842.0, 852.0, 912.0]) ...
];

%% ===================== 2. Built-in Constants =====================
approx_a = 0.4;
approx_b = 400;

% 22 standard wavelengths (row 1~22)
lambda_std = [
    404.66; 435.83; 546.07; 576.96; 579.07; 585.25; 594.48; 614.31; ...
    696.54; 703.24; 724.52; 743.89; 750.39; 763.51; 772.42; 794.82; ...
    800.62; 811.53; 826.45; 842.46; 852.14; 912.30 ...
];

% Lamp segment row mapping (1-based)
lampMap = struct( ...
    'HG', [1,2,3], ...
    'NE', [4,5,6,7,8,10,11,12], ...
    'AR', [9,13,14,15,16,17,18,19,20,21,22] ...
);

% B column background color (decimal)
lampColor = struct( ...
    'HG', 15790320, ...
    'NE', 15921894, ...
    'AR', 16764095 ...
);

% Output path
desktopPath = fullfile(getenv('USERPROFILE'), 'Desktop');
dt = datetime('now', 'Format', 'yyyyMMdd');
excelFile = [char(dt), '.xlsx'];
filePath = fullfile(desktopPath, excelFile);

%% ===================== 3. Identify Bin Modes =====================
% Pre-allocate bin list (max 3 modes)
binList = cell(1, 3);
binCount = 0;

for f = 1:length(pixelData)
    bin_token = regexp(pixelData(f).dir, '\d+bin', 'match', 'once');
    if ~isempty(bin_token)
        if ~ismember(bin_token{1}, binList(1:binCount))
            binCount = binCount + 1;
            binList{binCount} = bin_token{1};
        end
    end
end
binList = binList(1:binCount);

if isempty(binList)
    error('No bin mode detected from folder names. Check pixelData naming.');
end

%% ===================== 4. Create Multi-Sheet Template =====================
for b = 1:length(binList)
    binMode = binList{b};
    sheetName = ['offset', num2str(offsetValue), '-', binMode];
    
    col_A = cell(22, 1);
    col_B = num2cell(lambda_std);
    col_C = cell(22, 1);
    col_D = cell(22, 1);
    
    T = table(col_A, col_B, col_C, col_D, ...
        'VariableNames', {'A_pixel', 'B_wavelength', 'C_calculated', 'D_residual'});
    
    if b == 1
        writetable(T, filePath, 'Sheet', sheetName);
    else
        writetable(T, filePath, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');
    end
    
    % Initialize E1:E4 as empty (replace xlswrite)
    writecell(cell(4, 1), filePath, 'Sheet', sheetName, 'Range', 'E1:E4');
end

%% ===================== 5. Fill Column A =====================
% Pre-allocate warning list
debugWarn = cell(1, length(pixelData));
warnCount = 0;

for f = 1:length(pixelData)
    item = pixelData(f);
    dir_name = item.dir;
    
    % Identify lamp type
    if contains(dir_name, 'HG')
        lampType = 'HG';
    elseif contains(dir_name, 'NE')
        lampType = 'NE';
    elseif contains(dir_name, 'AR')
        lampType = 'AR';
    else
        continue;
    end
    
    % Identify bin mode
    bin_token = regexp(dir_name, '\d+bin', 'match', 'once');
    if isempty(bin_token), continue; end
    binMode = bin_token{1};
    
    sheetName = ['offset', num2str(offsetValue), '-', binMode];
    rowIdx = lampMap.(lampType);
    
    % Read current sheet
    sheetData = readcell(filePath, 'Sheet', sheetName);
    
    % Pre-allocate missing row index
    missRow = zeros(1, length(rowIdx));
    missCount = 0;
    
    for k = 1:length(rowIdx)
        r = rowIdx(k);
        pix_val = item.pixels(k);
        
        if isnan(pix_val)
            missCount = missCount + 1;
            missRow(missCount) = r;
            sheetData{r, 1} = [];
        else
            sheetData{r, 1} = round(pix_val * 100) / 100;
        end
    end
    
    writecell(sheetData, filePath, 'Sheet', sheetName);
    
    % Record missing warning
    if missCount > 0
        missRow = missRow(1:missCount);
        warnMsg = sprintf('DEBUG[Data Error] Sheet:%s Lamp:%s Missing rows:%s', ...
            sheetName, lampType, num2str(missRow));
        warnCount = warnCount + 1;
        debugWarn{warnCount} = warnMsg;
        warning('%s', warnMsg);
    end
end

% Print all warnings
if warnCount > 0
    fprintf('\n===== Missing Data Warnings =====\n');
    for w = 1:warnCount
        fprintf('%s\n', debugWarn{w});
    end
end

%% ===================== 6. Fit, Fill, Format, Chart =====================
for b = 1:length(binList)
    binMode = binList{b};
    sheetName = ['offset', num2str(offsetValue), '-', binMode];
    fprintf('\n===== Processing Sheet: %s =====\n', sheetName);
    
    tbl = readtable(filePath, 'Sheet', sheetName, 'ReadVariableNames', true);
    rawA = tbl.A_pixel;
    B = tbl.B_wavelength;
    
    % Pre-allocate index arrays (max 22 rows)
    validIdx = zeros(1, 22);
    ovlpIdx  = zeros(1, 22);
    p_all = nan(22, 1);
    vCount = 0;
    oCount = 0;
    
    for i = 1:22
        val = rawA{i};
        if ischar(val) && strcmp(val, '不可用')
            oCount = oCount + 1;
            ovlpIdx(oCount) = i;
        elseif ~ismissing(val) && ~isempty(val)
            vCount = vCount + 1;
            validIdx(vCount) = i;
            p_all(i) = str2double(val);
        end
    end
    validIdx = validIdx(1:vCount);
    ovlpIdx  = ovlpIdx(1:oCount);
    
    if vCount < 5
        warning('%s: Less than 5 valid points. Result for demo only.', sheetName);
    end
    
    % 3rd-order polynomial fit
    pV = p_all(validIdx);
    lV = B(validIdx);
    coeff3 = polyfit(pV, lV, 3);
    a3 = coeff3(1);
    a2 = coeff3(2);
    a1 = coeff3(3);
    a0 = coeff3(4);
    
    % Calculate column C and D (valid + overlap rows only)
    lPredAll = cell(22, 1);
    resAbsAll = cell(22, 1);
    calcIdx = [validIdx, ovlpIdx];
    
    for ii = 1:length(calcIdx)
        idx = calcIdx(ii);
        lPredAll{idx} = polyval(coeff3, p_all(idx));
        resAbsAll{idx} = abs(B(idx) - lPredAll{idx});
    end
    
    % Write back C and D columns
    tbl.C_calculated = lPredAll;
    tbl.D_residual = resAbsAll;
    writetable(tbl, filePath, 'Sheet', sheetName, 'WriteMode', 'overwritesheet');
    
    % Write E1:E4 coefficients (replace xlswrite)
    coeff_cell = {a3; a2; a1; a0};
    writecell(coeff_cell, filePath, 'Sheet', sheetName, 'Range', 'E1:E4');
    
    % Linear fit and quality judgement
    coeffLin = polyfit(pV, lV, 1);
    a_lin = coeffLin(1);
    
    resV = zeros(1, vCount);
    for v = 1:vCount
        resV(v) = resAbsAll{validIdx(v)};
    end
    maxResV = max(resV);
    
    if all(resV < a_lin)
        judge_result = 'Pass';
    else
        judge_result = 'Fail';
    end
    
    fprintf('Linear slope a_lin = %.6f nm/pixel\n', a_lin);
    fprintf('Max absolute residual = %.2f nm\n', maxResV);
    fprintf('Quality judgement: %s\n', judge_result);
    
    %% ----- Excel Format + Linear Chart (formula only, no R^2) -----
    try
        excel = actxserver('Excel.Application');
        excel.Visible = false;
        wb = excel.Workbooks.Open(filePath);
        sh = wb.Sheets.Item(sheetName);
        
        % Full table center alignment
        sh.Range('A1:D22').HorizontalAlignment = -4108;
        sh.Range('A1:D22').VerticalAlignment = -4108;
        
        % B column color by lamp segment
        lamp_names = fieldnames(lampMap);
        for lm = 1:length(lamp_names)
            rs = lampMap.(lamp_names{lm});
            rng_str = sprintf('B%d:B%d', min(rs), max(rs));
            sh.Range(rng_str).Interior.Color = lampColor.(lamp_names{lm});
        end
        
        % Insert scatter chart + linear trendline (formula only, R^2 off)
        chart = wb.Charts.Add;
        chart.ChartType = 73; % xlXYScatter
        chart.SetSourceData(sh.Range(sprintf('A%d:B%d', min(validIdx), max(validIdx))));
        series = chart.SeriesCollection(1);
        tr = series.Trendlines.Add;
        tr.Type = 1; % xlLinear
        tr.DisplayEquation = true;
        tr.DisplayRSquared = false;
        chart.Name = [sheetName, '_LinearFit'];
        
        wb.Save;
        wb.Close;
        excel.Quit;
        fprintf('Excel format and chart set completed\n');
        
    catch ME
        fprintf('Auto format failed. Please do manually:\n');
        fprintf('1. Center align all cells\n');
        fprintf('2. Set B column background by lamp segment\n');
        fprintf('3. Insert scatter chart with linear trendline, show formula only\n');
        fprintf('Error: %s\n', ME.message);
    end
    
    
end

fprintf('\nAll processing completed. Output Excel saved to desktop\n');
fprintf('File path: %s\n', filePath);
