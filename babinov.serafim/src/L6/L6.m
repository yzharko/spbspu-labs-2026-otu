clc;
clear;
close all;

s = tf('s');
a = 1;
b = 1;
Go = 1/(s^2 + a*s + b);
sigma_req = 5;
t_p_req = 2.5;

% 1. Исходная замкнутая система
G_closed = feedback(Go, 1);
S0 = stepinfo(G_closed);
fprintf('Исходная система. Полюса:\n');
disp(pole(G_closed));
fprintf('Устойчива: %d, tп = %.4f с, sigma = %.2f %%, ошибка = %.4f\n\n', ...
    isstable(G_closed), S0.SettlingTime, S0.Overshoot, 1 - dcgain(G_closed));

% 2. Первый метод Зиглера-Николса
w = sqrt(b - (a/2)^2);
t_inflection = atan(2*w/a)/w;
slope = exp(-a*t_inflection/2)*sin(w*t_inflection)/w;
K_object = dcgain(Go);

y_inflection = K_object*(1 - a*slope);
L = t_inflection - y_inflection/slope;
T = K_object/slope;

Kp = 1.2*T/L;
Ti = 2*L;
Td = 0.5*L;
Ki = Kp/Ti;
Kd = Kp*Td;
C_zn = pid(Kp, Ki, Kd);
fprintf('Первый метод Зиглера-Николса: L = %.4f с, T = %.4f с\n', L, T);
fprintf('Kp = %.4f, Ti = %.4f с, Td = %.4f с, Ki = %.4f, Kd = %.4f\n', ...
    Kp, Ti, Td, Ki, Kd);

[y, t] = step(Go);
figure('Theme', 'light');
plot(t, y, [L L+T], [0 K_object], '--', 'LineWidth', 1.2);
grid on;
yline(K_object, ':', 'HandleVisibility', 'off');
xlabel('Время, с');
ylabel('Выход объекта');
title('Переходная характеристика объекта и касательная');
legend('Объект', 'Касательная в точке перегиба', 'Location', 'best');

% 3. Замкнутая система с регулятором Зиглера-Николса
G_zn = feedback(C_zn*Go, 1);
S_zn = stepinfo(G_zn);
fprintf('Устойчива: %d, ошибка = %.4f\n', isstable(G_zn), 1 - dcgain(G_zn));
fprintf('tп = %.4f с (<= %.1f с: %d), sigma = %.2f %% (<= %.0f %%: %d)\n\n', ...
    S_zn.SettlingTime, t_p_req, S_zn.SettlingTime <= t_p_req, ...
    S_zn.Overshoot, sigma_req, S_zn.Overshoot <= sigma_req);

% 4. Автоматическая настройка ПИД-регулятора
C_auto = pidtune(Go, 'PID');
G_auto = feedback(C_auto*Go, 1);
S_auto = stepinfo(G_auto);
fprintf('Автонастройка: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', ...
    C_auto.Kp, C_auto.Ki, C_auto.Kd);
fprintf('Устойчива: %d, ошибка = %.4f\n', isstable(G_auto), 1 - dcgain(G_auto));
fprintf('tп = %.4f с (<= %.1f с: %d), sigma = %.2f %% (<= %.0f %%: %d)\n\n', ...
    S_auto.SettlingTime, t_p_req, S_auto.SettlingTime <= t_p_req, ...
    S_auto.Overshoot, sigma_req, S_auto.Overshoot <= sigma_req);

% 5. Сравнение настроек
[GM_zn, PM_zn] = margin(C_zn*Go);
[GM_auto, PM_auto] = margin(C_auto*Go);
fprintf('Запасы устойчивости по усилению (дБ) и фазе (град):\n');
fprintf('Зиглер-Николс: %.2f дБ, %.2f град.\n', 20*log10(GM_zn), PM_zn);
fprintf('pidtune: %.2f дБ, %.2f град.\n', 20*log10(GM_auto), PM_auto);

figure('Theme', 'light');
step(G_zn, G_auto);
grid on;
title('Переходные характеристики замкнутых систем с ПИД');
legend('Зиглер-Николс', 'Автонастройка pidtune', 'Location', 'best');
