clear; clc; close all;

s = tf("s");

W_init = 1500 / ( s * (0.05 * s + 1) * (0.02 * s + 1));

W_closed = feedback(W_init, 1);
disp(M(W_closed)) % 1.33

d = cell2mat(W_closed.Denominator);
disp(mat2str(d));       % [0.001 0.07 1 1500]
disp(mat2str(gurv(d))); % [0.07 -1.43 -2145] - неустойчива
figure("Name", "Исходная система");
step(W_closed);

w_range = logspace(-2, 4, 1000);
detailed_bode(W_init, "Исходные ЛАЧХ и ЛФЧХ");

% По рисунку
w = [0.015 1.9 89 1500];
W_cor = (1 / (s / w(1) + 1)) ... % апериодическое -20 дБ
      * (s / w(2) + 1) ...       % форсирующее +20 дБ
      * (s / w(3) + 1) ...       % форсирующее +20 дБ
      * (1 / (s / w(4) + 1));    % апериодическое -20 дБ
disp(W_cor); % num: 22.5 2045 3804; den: 169.1 253650 3804.7
detailed_bode(W_cor, "Регулятор");
W_fin = W_cor * W_init;
detailed_bode(W_fin, "Отрегулированные ЛАЧХ и ЛФЧХ");
W_fin_closed = feedback(W_fin, 1);
figure("Name", "Отрегулированная система");
step(W_fin_closed);
disp(M(W_fin_closed)); % 1.2886
disp(stepinfo(W_fin_closed));
%      RiseTime: 0.1079
% TransientTime: 1.0820
%  SettlingTime: 1.0820
%   SettlingMin: 0.9180
%   SettlingMax: 1.2847
%     Overshoot: 28.4703
%    Undershoot: 0
%          Peak: 1.2847
%      PeakTime: 0.2726

function g = gurv(a)
    gurvMatrix = [
        a(2), a(4),    0
        a(1), a(3),    0
           0, a(2), a(4)
    ];
    first = a(2);
    second = det(gurvMatrix(1:2, 1:2));
    third = second * a(4);
    g = [first, second, third];
end

function m = M(W)
    [mag, ~, ~] = bode(W);
    m = max(mag);
end

function detailed_bode(W, name)
    [mag, phase, w] = bode(W, logspace(-2, 4, 1000));

    figure("Name", name);
    subplot(2,1,1);
    mag_dB = 20*log10(squeeze(mag));
    semilogx(w, mag_dB);
    grid on;
    ylabel('Magnitude (dB)');
    
    mag_min = floor(min(mag_dB)/20)*20;
    mag_max = ceil(max(mag_dB)/20)*20;
    yticks(mag_min:20:mag_max);
    ylim([mag_min mag_max]);
    
    title('Bode Plot');
    
    subplot(2,1,2);
    phase_deg = squeeze(phase);
    semilogx(w, phase_deg);
    grid on;
    ylabel('Phase (deg)');
    xlabel('Frequency (rad/s)');
    
    ph_min = floor(min(phase_deg)/10)*10;
    ph_max = ceil(max(phase_deg)/10)*10;
    yticks(ph_min:10:ph_max);
    ylim([ph_min ph_max]);
end
