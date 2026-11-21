clc;
clear;
close all;

% 1. Исходная система
s = tf('s');
K_open = 282;
T_open = 0.1;
M_limit = 1.5;
overshoot_limit = 37;
peak_tolerance = 1e-8;

W_open = K_open / (s*(T_open*s + 1));
W_closed = feedback(W_open, 1);
info_open = stepinfo(W_closed);
[~, ~, ~, wc_initial] = margin(W_open);

if isstable(W_closed)
    fprintf('Исходная система устойчива\n');
else
    fprintf('Исходная система неустойчива\n');
end

fprintf('Tпп = %.4f с, перерегулирование = %.2f %%, M = %.4f\n', ...
    info_open.SettlingTime, info_open.Overshoot, getPeakGain(W_closed, peak_tolerance));
fprintf('Частота среза = %.2f рад/с', wc_initial);

% 2. Расчёт ПКУ и подбор частоты среза
phase_margin_required = asind(1/M_limit);

% Tпп в варианте не задано, поэтому начинаем с исходной частоты среза.
% Остальные частоты выбраны вручную.
wc_trials = [wc_initial 100 200 203 210];
best_settling_time = Inf;
wc = NaN;
fprintf('\nПодбор частоты среза (M <= %.2f, перерегулирование <= %.0f %%)\n', ...
    M_limit, overshoot_limit);

for wc_trial = wc_trials
    [magnitude, phase] = bode(W_open, wc_trial);
    phase_lead = phase_margin_required - (180 + phase);
    alpha = (1 + sind(phase_lead)) / (1 - sind(phase_lead));

    T2_trial = 1 / (wc_trial*sqrt(alpha));
    T1_trial = alpha*T2_trial;

    L_open_wc = 20*log10(magnitude);
    L_lead_wc = 20*log10(sqrt(alpha));
    K_trial = 10^((-L_open_wc - L_lead_wc)/20);

    W_ku_trial = K_trial*(T1_trial*s + 1)/(T2_trial*s + 1);
    W_closed_trial = feedback(W_open*W_ku_trial, 1);
    info_trial = stepinfo(W_closed_trial);
    M_trial = getPeakGain(W_closed_trial, peak_tolerance);
    fits = isstable(W_closed_trial) && M_trial <= M_limit && ...
        info_trial.Overshoot <= overshoot_limit;

    fprintf('wc = %7.2f: Tпп = %.4f с, sigma = %.2f %%, M = %.4f, подходит: %d\n', ...
        wc_trial, info_trial.SettlingTime, info_trial.Overshoot, M_trial, fits);
    if fits && info_trial.SettlingTime < best_settling_time
        wc = wc_trial;
        best_settling_time = info_trial.SettlingTime;
        T1 = T1_trial;
        T2 = T2_trial;
        K = K_trial;
        W_ku = W_ku_trial;
    end
end

if isnan(wc)
    error('Среди проверенных частот нет подходящей.');
end

fprintf('\nВыбрано wc = %.2f рад/с, K = %.4f, T1 = %.6f с, T2 = %.6f с\n', ...
    wc, K, T1, T2);
fprintf('Частоты сопряжения ПКУ: %.2f и %.2f рад/с\n', 1/T1, 1/T2);
disp('Передаточная функция ПКУ:');
disp(W_ku);

% 3. Исходные и желаемые ЛАЧХ и ЛФЧХ
W_desired = W_open*W_ku;
W_closed_corrected = feedback(W_desired, 1);

figure('Theme', 'light');
bode(W_open, W_desired);
grid on;
legend('Исходная', 'Желаемая', 'Location', 'best');
title('Частотные характеристики разомкнутых систем');

% 4. Геометрическое вычитание асимптотических ЛАЧХ
w = logspace(0, 4, 2000);
L_open = 20*log10(K_open./w) - 20*max(0, log10(T_open*w));
L_desired = 20*log10(K*K_open./w) + 20*max(0, log10(T1*w)) ...
    - 20*max(0, log10(T_open*w)) - 20*max(0, log10(T2*w));
L_ku = L_desired - L_open;

figure('Theme', 'light');
subplot(2, 1, 1);
semilogx(w, L_open, w, L_desired, 'LineWidth', 1.3);
grid on;
ylabel('L, дБ');
title('Асимптотические ЛАЧХ');
legend('Исходная', 'Желаемая', 'Location', 'best');

subplot(2, 1, 2);
semilogx(w, L_ku, 'LineWidth', 1.3);
xline(1/T1, ':');
xline(1/T2, ':');
grid on;
xlabel('\omega, рад/с');
ylabel('L_{ку}, дБ');
title('ЛАЧХ ПКУ как разность желаемой и исходной');

% 5. Проверка скорректированной замкнутой системы
[~, phase_margin, ~, wc_actual] = margin(W_desired);
fprintf('Фактическая частота среза = %.2f рад/с, запас по фазе = %.2f град\n', ...
    wc_actual, phase_margin);

figure('Theme', 'light');
step(W_closed, W_closed_corrected);
grid on;
legend('Исходная', 'Скорректированная', 'Location', 'best');
title('Переходные характеристики замкнутых систем');

figure('Theme', 'light');
bode(W_closed_corrected);
grid on;
title('ЛАЧХ и ЛФЧХ скорректированной замкнутой системы');
