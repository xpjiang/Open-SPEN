clear all
close all
clc



% Calculate SPEN-conditions: 
% sweepBw*rf_dur = gexc.amplitude*rf_dur*fov = gacq.amplitude*Ny*mr.calcDuration(gacq)*fov




rf_dur = 4e-3; % 4ms
delay = 1.1000e-04; %0.11ms
sweepBw=25000; % 25KHz
ang = 90;
n_fac = 40;
sys = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,'B0',3);

% Create chirped
% RF and excitation gradient
rf = makeChirpedRfPulse('duration',rf_dur,'delay',delay,'bandwidth',sweepBw, ...
    'ang',90,'n_fac',40,'system',sys); %+gexc.riseTime nachträglich geändert. Messungen waren ohne-> leichte Verzerrungen

% 此射频对应的FOV和编码梯度强度关系：sweepBw = gexc.amplitude*fov;



rfsignal = rf.signal; %1us一个点

% normalize
maxS = max(abs(rfsignal));

rfsignal = rfsignal/maxS * 117960;

% 创建画布
figure('Position', [100 100 800 600])

% 上子图：幅度
subplot(2,1,1)
plot(abs(rfsignal), 'LineWidth', 1.5)
title('RF Signal Amplitude')
ylabel('Amplitude (a.u.)')
grid on

% 下子图：相位（弧度转角度）
subplot(2,1,2)
plot(rad2deg(angle(rfsignal)), 'LineWidth', 1.5)
title('RF Signal Phase')
xlabel('Time (us)')
ylabel('Phase (degrees)')
grid on

complex2txt(rfsignal, 'test.txt');


function complex2txt(complexArray, filename)
    % 提取实部和虚部并转换为int32
    realPart = int32(bankers_round(real(complexArray)));
    imagPart = int32(bankers_round(imag(complexArray)));
    
    % 组合数据矩阵
    data = [realPart(:), imagPart(:)];
    
    % 写入文件
    fid = fopen(filename, 'w');
    fprintf(fid, '%d,%d\n', data');  % 写入数据
    fclose(fid);
end

function y = bankers_round(x, n)
    % 输入: x为数值或数组, n为保留小数位数(默认为0)
    % 输出: 银行家舍入结果
    
    if nargin < 2, n = 0; end
    scale = 10^n;
    x_scaled = x * scale;
    
    % 分离整数和小数部分
    integer_part = fix(x_scaled);
    decimal_part = abs(x_scaled - integer_part);
    
    % 判断舍入条件
    mask_gt5 = (decimal_part > 0.5);
    mask_eq5 = (abs(decimal_part - 0.5) < eps);
    mask_odd = (mod(integer_part, 2) ~= 0);
    
    % 应用银行家规则
    y = integer_part + mask_gt5 + (mask_eq5 & mask_odd);
    y = sign(x_scaled) .* abs(y) / scale;
end