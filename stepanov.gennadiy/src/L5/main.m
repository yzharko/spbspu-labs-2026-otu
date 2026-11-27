
clear;
clc;
close all;

num = 100;
den = [0.00002 0.0032 0.13 1 0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);
char_eq = den + [0 0 0 0 num];
a4 = char_eq(1);
a3 = char_eq(2);
a2 = char_eq(3);
a1 = char_eq(4);
a0 = char_eq(5);
D1 = a3;
D2 = a3*a2 - a4*a1;
D3 = a3*a2*a1 - a3^2*a0 - a4*a1^2;
D4 = a0 * D3;
fprintf('Определители Гурвица:\n');
fprintf('D1 = %.6f\nD2 = %.6f\nD3 = %.6f\nD4 = %.6f\n', ...
        D1, D2, D3, D4);
figure(1);
step(W_closed);
grid on;
title('Переходная характеристика исходной САУ');

M_req = 1.3;
Tpp_req = 0.8;
sigma_req = 25;
gamma_req = asind(1 / M_req);
wc = 10;
[mag_wc, phase_wc] = bode(W_open, wc);
mag_wc = squeeze(mag_wc);
phase_wc = squeeze(phase_wc);
gamma_base = 180 + phase_wc;
delta_phi = gamma_req - gamma_base;
alpha = (1 + sind(delta_phi)) / (1 - sind(delta_phi));

T2 = 1 / (wc * sqrt(alpha));
T1 = alpha * T2;
K = 1 / (mag_wc * sqrt(alpha));
W_ku = K * tf([T1 1], [T2 1]);
W_desired = W_open * W_ku;
w = logspace(-1, 4, 2000);
[mag_open, phase_open] = bode(W_open, w);
mag_open = squeeze(mag_open);
phase_open = squeeze(phase_open);
[mag_des, phase_des] = bode(W_desired, w);
mag_des = squeeze(mag_des);
phase_des = squeeze(phase_des);
L_open = 20 * log10(mag_open);
L_des = 20 * log10(mag_des);
figure(2);
semilogx(w, L_open, 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, L_des, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Частота, рад/с');
ylabel('Амплитуда, дБ');
title('ЛАЧХ исходной и желаемой систем');
legend('Исходная', 'Желаемая', 'Location', 'southwest');
hold off;
figure(3);
semilogx(w, phase_open, 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, phase_des, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Частота, рад/с');
ylabel('Фаза, град');
title('ЛФЧХ исходной и желаемой систем');
legend('Исходная', 'Желаемая', 'Location', 'southwest');
hold off;

L_ku = L_des - L_open;
L_ku_asym = 20 * log10(K) ...
    + 20 * log10(max(1, T1*w)) ...
    - 20 * log10(max(1, T2*w));
figure(4);
semilogx(w, L_ku, 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, L_ku_asym, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Частота, рад/с');
ylabel('L_{ку}(\omega), дБ');
title('ЛАЧХ корректирующего устройства');
legend('Точная ЛАЧХ ПКУ', 'Асимптотическая ЛАЧХ', ...
       'Location', 'southwest');
hold off;
fprintf('\nПараметры ПКУ:\n');
fprintf('K = %.6f\n', K);
fprintf('T1 = %.6f с\n', T1);
fprintf('T2 = %.6f с\n', T2);
fprintf('Частота сопряжения 1 = %.4f рад/с\n', 1/T1);
fprintf('Частота сопряжения 2 = %.4f рад/с\n', 1/T2);
disp('Передаточная функция ПКУ:');
W_ku

W_closed_final = feedback(W_desired, 1);
p_final = pole(W_closed_final);
disp(p_final);
t = 0:0.001:1.5;
h = step(W_closed_final, t);
y_ust = dcgain(W_closed_final);
figure(5);
plot(t, h, 'b', 'LineWidth', 1.5);
hold on;
plot(t, 1.02*y_ust*ones(size(t)), 'r--');
plot(t, 0.98*y_ust*ones(size(t)), 'r--');
plot(t, y_ust*ones(size(t)), 'k:');
grid on;
xlabel('Время, с');
ylabel('h(t)');
title('Переходная характеристика скорректированной САУ');
legend('h(t)', '+2%', '-2%', ...
       'Установившееся значение', ...
       'Location', 'southeast');
hold off;
[mag_closed, phase_closed] = bode(W_closed, w);
mag_closed = squeeze(mag_closed);
phase_closed = squeeze(phase_closed);
[mag_final, phase_final] = bode(W_closed_final, w);
mag_final = squeeze(mag_final);
phase_final = squeeze(phase_final);
figure(6);
subplot(2,1,1);
semilogx(w, 20*log10(mag_closed), 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, 20*log10(mag_final), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Частота, рад/с');
ylabel('Амплитуда, дБ');
title('ЛАЧХ замкнутых систем');
legend('Исходная', 'Скорректированная', ...
       'Location', 'southwest');
hold off;
subplot(2,1,2);
semilogx(w, phase_closed, 'b', 'LineWidth', 1.5);
hold on;
semilogx(w, phase_final, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Частота, рад/с');
ylabel('Фаза, град');
title('ЛФЧХ замкнутых систем');
legend('Исходная', 'Скорректированная', ...
       'Location', 'southwest');
hold off;
