clc;
clear;
close all;

num = 1500;
Ta = 0.05;
Tb = 0.02;
den = [Ta*Tb Ta+Tb 1 0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

D1 = den(2);
D2 = den(2)*den(3) - den(1)*num;
D3 = num*D2;
[~, Pm_open, ~, wc_open] = margin(W_open);
fprintf('Вариант 7. Исходная замкнутая система\n');
fprintf('Определители Гурвица: %.3f, %.3f, %.3f.\n', D1, D2, D3);
fprintf('Разомкнутый контур до коррекции: частота среза %.2f рад/с, запас по фазе %.2f град.\n', wc_open, Pm_open);
if ~isstable(W_closed)
    fprintf('Система неустойчива. Время установления и перерегулирование не определены.\n');
end

M_req = 1.5;
sigma_req = 37;
band = 0.02;
gamma_req = asind(1/M_req);
best_time = Inf;

for w_try = 1:0.25:90
    gamma_base = 90 - atand(Ta*w_try) - atand(Tb*w_try);
    delta_phi = gamma_req - gamma_base;
    if delta_phi <= 0 || delta_phi >= 90
        continue;
    end
    q = sqrt((1 + sind(delta_phi))/(1 - sind(delta_phi)));
    Tz = q/w_try;
    Tp = 1/(q*w_try);
    mag_open = abs(evalfr(W_open, 1i*w_try));
    gain_des = num/(mag_open*q);
    W_des_try = tf(gain_des*[Tz 1], conv(den, [Tp 1]));
    W_try = feedback(W_des_try, 1);
    decay = -max(real(pole(W_try)));
    if decay <= 0 || getPeakGain(W_try, 1e-4) > M_req
        continue;
    end
    t_try = (0:0.00025:8/decay)';
    y_try = step(W_try, t_try);
    info_try = stepinfo(y_try, t_try, 1, 'SettlingTimeThreshold', band);
    if info_try.Overshoot <= sigma_req && info_try.SettlingTime < best_time
        best_time = info_try.SettlingTime;
        wc = w_try;
        T1 = Tz;
        T2 = Tp;
        K_des = gain_des;
    end
end
if isinf(best_time)
    error('Для заданных ограничений корректор на выбранной сетке не найден.');
end

w = unique([logspace(-1, 3.5, 1000) 1/Ta 1/Tb 1/T1 1/T2 wc])';
bend = @(T) 20*log10(max(1, T*w));
L_open_as = 20*log10(num./w) - bend(Ta) - bend(Tb);
L_des_as = 20*log10(K_des./w) + bend(T1) - bend(Ta) - bend(Tb) - bend(T2);
L_ku_as = L_des_as - L_open_as;
K = 10^(L_ku_as(1)/20);
W_ku = K*tf([T1 1], [T2 1]);
W_final = W_ku*W_open;
W_closed_final = feedback(W_final, 1);

[~, Pm, ~, wcp] = margin(W_final);
[M_fact, w_peak] = getPeakGain(W_closed_final, 1e-4);
t = (0:0.00025:8/min(-real(pole(W_closed_final))))';
y = step(W_closed_final, t);
info = stepinfo(y, t, 1, 'SettlingTimeThreshold', band);
fprintf('\nСкорректированная система\n');
fprintf('Частота среза = %.2f рад/с. Запас по фазе = %.2f град.\n', wcp, Pm);
fprintf('Частоты сопряжения исходной системы: %.2f и %.2f рад/с.\n', 1/Ta, 1/Tb);
fprintf('Частоты сопряжения ПКУ: %.2f и %.2f рад/с.\n', 1/T1, 1/T2);
fprintf('Wку(s) = %.7f * (%.7f s + 1) / (%.7f s + 1)\n', K, T1, T2);
fprintf('M = %.3f (допустимо до %.2f). Перерегулирование = %.2f %% (до %.0f %%).\n', M_fact, M_req, info.Overshoot, sigma_req);
fprintf('Время установления в полосе 2 %% = %.4f с. Частота выбрана по минимуму времени с шагом 0.25 рад/с.\n', info.SettlingTime);
if isstable(W_closed_final) && M_fact <= M_req && info.Overshoot <= sigma_req
    fprintf('Замкнутая система устойчива. Оба ограничения выполнены.\n');
end

newFigure('Переходные характеристики');
t_before = (0:0.00025:0.20)';
subplot(2,1,1);
plot(t_before, step(W_closed, t_before), 'Color', [0.63 0.20 0.10], 'LineWidth', 1.6);
grid on; title('Исходная замкнутая система'); xlabel('t, с'); ylabel('y(t)');
subplot(2,1,2);
plot(t, y, 'Color', [0 0.32 0.65], 'LineWidth', 1.6);
xlim([0 2*info.SettlingTime]); yline(1-band, ':k'); yline(1+band, ':k');
xline(info.SettlingTime, '--k', 'T_{пп}');
grid on; title('Скорректированная замкнутая система'); xlabel('t, с'); ylabel('y(t)');

[mag_open, phase_open] = bode(W_open, w);
[mag_des, phase_des] = bode(W_final, w);
newFigure('Частотный синтез ПКУ');
subplot(2,2,1);
semilogx(w, 20*log10(squeeze(mag_open)), '--', w, 20*log10(squeeze(mag_des)), '-', 'LineWidth', 1.5);
grid on; yline(0, ':k'); xline(wc, ':k'); title('ЛАЧХ разомкнутого контура');
xlabel('\omega, рад/с'); ylabel('L, дБ'); legend('Исходная', 'Желаемая', 'Location', 'southwest');
subplot(2,2,3);
semilogx(w, squeeze(phase_open), '--', w, squeeze(phase_des), '-', 'LineWidth', 1.5);
grid on; xline(wc, ':k'); title('ЛФЧХ разомкнутого контура');
xlabel('\omega, рад/с'); ylabel('\phi, град'); legend('Исходная', 'Желаемая', 'Location', 'southwest');
subplot(2,2,2);
semilogx(w, L_open_as, '--', w, L_des_as, '-', 'LineWidth', 1.5);
grid on; title('Асимптотические ЛАЧХ'); xlabel('\omega, рад/с'); ylabel('L, дБ');
legend('Исходная', 'Желаемая', 'Location', 'southwest');
subplot(2,2,4);
semilogx(w, L_ku_as, 'Color', [0 0.32 0.65], 'LineWidth', 1.6);
grid on; xline(1/T1, ':k'); xline(1/T2, ':k'); title('ЛАЧХ ПКУ по разности асимптот');
xlabel('\omega, рад/с'); ylabel('L_{ку}, дБ');

[mag_closed, phase_closed] = bode(W_closed_final, w);
newFigure('Частотные характеристики замкнутой системы');
subplot(2,1,1);
semilogx(w, 20*log10(squeeze(mag_closed)), 'LineWidth', 1.6); hold on;
plot(w_peak, 20*log10(M_fact), 'ko', 'MarkerFaceColor', 'k');
yline(20*log10(M_req), '--k', 'Предел M');
grid on; title('ЛАЧХ после замыкания обратной связи'); xlabel('\omega, рад/с'); ylabel('L, дБ');
subplot(2,1,2);
semilogx(w, squeeze(phase_closed), 'LineWidth', 1.6);
grid on; title('ЛФЧХ после замыкания обратной связи'); xlabel('\omega, рад/с'); ylabel('\phi, град');

function newFigure(name)
    f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'Position', [0.10 0.12 0.80 0.72]);
    if isprop(f, 'Theme')
        f.Theme = 'light';
    end
    set(f, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
    set(f, 'DefaultLegendColor', 'w', 'DefaultLegendTextColor', 'k', 'DefaultAxesColorOrder', [0 0.32 0.65; 0.65 0.18 0.04]);
end
