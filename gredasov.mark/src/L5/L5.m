clear;
clc;
close all;

num = 32;
den = [0.001 0.035 0.35 1];

W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = den + [0 0 0 num];

a3 = char_eq(1);
a2 = char_eq(2);
a1 = char_eq(3);
a0 = char_eq(4);

D1 = a2;
D2 = a2*a1 - a3*a0;
D3 = a0*D2;

fprintf('Определители Гурвица: D1=%.3f  D2=%.3f  D3=%.3f\n', D1, D2, D3);

[Am, Pm, Wcg, Wcp] = margin(W_open);

figure('Name','Исходная САУ');
subplot(2,1,1);
step(W_closed);
title('Переходная функция исходной САУ');
grid on;

subplot(2,1,2);
bode(W_open);
grid on;
title('ЛАЧХ/ЛФЧХ исходной разомкнутой САУ');

M_req = 1.19;
Tpp_req = 0.6;
sigma_req = 20;

gamma_req = asind(1/M_req);
wc_req = 8;
phi_req = -180 + gamma_req;

[T1, T2] = calculate_T(num, M_req, wc_req);

fprintf('\nПараметры желаемой системы:\n');
fprintf('gamma = %.2f град\n', gamma_req);
fprintf('wc = %.2f рад/с\n', wc_req);
fprintf('T1 = %.6f с\n', T1);
fprintf('T2 = %.6f с\n', T2);

W_zh = tf(num, conv([T1 1], [T2 1]));

disp('Передаточная функция желаемой системы:');
W_zh

w = logspace(-2, 3, 1000);

[mag_i, phase_i] = bode(W_open, w);
mag_i = squeeze(mag_i);
phase_i = squeeze(phase_i);

[mag_zh, phase_zh] = bode(W_zh, w);
mag_zh = squeeze(mag_zh);
phase_zh = squeeze(phase_zh);

L_i = 20*log10(mag_i);
L_zh = 20*log10(mag_zh);

figure('Name','Исходная и желаемая частотные характеристики');

subplot(2,1,1);
semilogx(w, L_i, 'b', w, L_zh, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('L(\omega), дБ');
title('ЛАЧХ исходной и желаемой систем');
legend('Исходная система', 'Желаемая система', 'Location', 'best');

subplot(2,1,2);
semilogx(w, phase_i, 'b', w, phase_zh, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('\phi(\omega), град');
title('ЛФЧХ исходной и желаемой систем');
legend('Исходная система', 'Желаемая система', 'Location', 'best');

w_i = [5 10 20];
w_zh = [1/T1 1/T2];

L0 = 20*log10(num);

L_i_asym = L0 - 20*max(0, log10(w/5)) - 20*max(0, log10(w/10)) - 20*max(0, log10(w/20));
L_zh_asym = L0 - 20*max(0, log10(w/(1/T1))) - 20*max(0, log10(w/(1/T2)));
L_ku_asym = L_zh_asym - L_i_asym;

figure('Name','Асимптотические ЛАЧХ исходной системы, желаемой системы и ПКУ');

semilogx(w, L_i_asym, 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, L_zh_asym, 'r--', 'LineWidth', 1.5);
semilogx(w, L_ku_asym, 'g-.', 'LineWidth', 1.5);

grid on;
xlabel('\omega, рад/с');
ylabel('L(\omega), дБ');
title('Асимптотические ЛАЧХ исходной системы, желаемой системы и ПКУ');
legend('Исходная система', 'Желаемая система', ...
    'Корректирующее устройство (ПКУ)', 'Location', 'best');
hold off;

W_ku = minreal(W_zh / W_open);
disp('Передаточная функция корректирующего устройства ПКУ:');
W_ku

W_final = W_ku * W_open;
W_closed_final = feedback(W_final, 1);

[Am2, Pm2, Wcg2, Wcp2] = margin(W_final);
M_fact = 1 / sind(Pm2);
fprintf('Фактическое M = 1/sin(Pm) = %.2f (требуется <= %.2f)\n', M_fact, M_req);

info = stepinfo(W_closed_final);
Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;
yesNo = @(b) subsref({'НЕТ','ДА'}, struct('type','{}','subs',{{b+1}}));
fprintf('Соответствие: Tпп %.3f <= %.1f ? %s    sigma %.1f <= %.0f ? %s\n', ...
    Tpp_fact, Tpp_req, yesNo(Tpp_fact<=Tpp_req), sigma_fact, sigma_req, yesNo(sigma_fact<=sigma_req));

figure('Name','ЛАЧХ/ЛФЧХ систем');

subplot(2,1,1);
[mag_open, phase_open, wout] = bode(W_open, w);
[mag_ku, phase_ku] = bode(W_ku, w);
[mag_final, phase_final] = bode(W_final, w);

semilogx(w, 20*log10(squeeze(mag_open)), 'b', ...
         w, 20*log10(squeeze(mag_ku)), 'r--', ...
         w, 20*log10(squeeze(mag_final)), 'g-.', ...
         'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('L(\omega), дБ');
legend('Исходная','ПКУ', ...
       'Скорректированная','Location','best');
title('ЛАЧХ — исходная, ПКУ, скорректированная');

subplot(2,1,2);
semilogx(w, squeeze(phase_open), 'b', ...
         w, squeeze(phase_ku), 'r--', ...
         w, squeeze(phase_final), 'g-.', ...
         'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('\phi(\omega), град');
legend('Исходная','ПКУ', ...
       'Скорректированная','Location','best');
title('ЛФЧХ — исходная, ПКУ, скорректированная');

t = 0:0.01:1;

figure('Name','Переходные характеристики');
step(W_closed, W_closed_final, t);
grid on;
legend('Исходная САУ','Скорректированная САУ','Location','best');
title('Переходные характеристики замкнутых САУ');

function [T1, T2] = calculate_T(K, M_req, wc_req)
    gamma_req = asind(1 / M_req);

    alpha = 180 - gamma_req;

    S = K * sind(alpha);
    P = 1 - K * cosd(alpha);

    D = S^2 - 4 * P;

    if D < 0
        error('Нет действительных решений для T1 и T2');
    end

    x1 = (S + sqrt(D)) / 2;
    x2 = (S - sqrt(D)) / 2;

    T1 = x1 / wc_req;
    T2 = x2 / wc_req;
end
