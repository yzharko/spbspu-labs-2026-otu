clear;
clc;
close all;

s = tf('s');

Go = 1 / ((s + 0.1) * (s + 0.5) * (s + 2));

sigma_req = 10;
t_p_req = 2;

G_closed = feedback(Go, 1);
S0 = stepinfo(G_closed);

fprintf('Исходная система\n');
fprintf('Полюса:\n');
disp(pole(G_closed));

fprintf('Устойчивость: %d\n', isstable(G_closed));
fprintf('t_p = %.4f с\n', S0.SettlingTime);
fprintf('sigma = %.2f %%\n', S0.Overshoot);
fprintf('Ошибка = %.4f\n\n', 1 - dcgain(G_closed));

figure;
step(G_closed);
grid on;
title('Переходный процесс исходной системы');

[Kcr, ~, Wcr, ~] = margin(Go);

Pcr = 2*pi / Wcr;

Kp = 0.6 * Kcr;
Ti = 0.5 * Pcr;
Td = 0.125 * Pcr;

Ki = Kp / Ti;
Kd = Kp * Td;

fprintf('Метод Зиглера-Николса\n');
fprintf('Kcr = %.4f\n', Kcr);
fprintf('Pcr = %.4f с\n', Pcr);
fprintf('Kp = %.4f\n', Kp);
fprintf('Ti = %.4f с\n', Ti);
fprintf('Td = %.4f с\n\n', Td);

C = pid(Kp, Ki, Kd);

G_pid = feedback(C * Go, 1);
Sp = stepinfo(G_pid);

fprintf('Система с ПИД-регулятором\n');
fprintf('Устойчивость: %d\n', isstable(G_pid));
fprintf('t_p = %.4f с, требование выполнено: %d\n', ...
    Sp.SettlingTime, Sp.SettlingTime <= t_p_req);
fprintf('sigma = %.2f %%, требование выполнено: %d\n', ...
    Sp.Overshoot, Sp.Overshoot <= sigma_req);
fprintf('Ошибка = %.4f\n', 1 - dcgain(G_pid));

figure;
step(G_pid);
grid on;
title('Переходный процесс системы с ПИД-регулятором');

Ca = pidtune(Go, 'PID');

Ga = feedback(Ca * Go, 1);
Sa = stepinfo(Ga);

fprintf('\nАвтоматическая настройка pidtune\n');
fprintf('Kp = %.4f\n', Ca.Kp);
fprintf('Ki = %.4f\n', Ca.Ki);
fprintf('Kd = %.4f\n', Ca.Kd);
fprintf('t_p = %.4f с\n', Sa.SettlingTime);
fprintf('sigma = %.2f %%\n', Sa.Overshoot);
fprintf('Ошибка = %.4f\n', 1 - dcgain(Ga));
fprintf('Устойчивость: %d\n', isstable(Ga));

figure;
step(G_pid, Ga);
grid on;

title('Сравнение методов настройки ПИД-регулятора');
xlabel('Время, с');
ylabel('Амплитуда');

legend('Зиглер-Николс', 'Автонастройка pidtune', ...
    'Location', 'best');

[GM_zn, PM_zn] = margin(C * Go);
[GM_auto, PM_auto] = margin(Ca * Go);

fprintf('\nЗапасы устойчивости\n');
fprintf('Зиглер-Николс: GM = %.4f, PM = %.2f град.\n', ...
    GM_zn, PM_zn);
fprintf('pidtune: GM = %.4f, PM = %.2f град.\n', ...
    GM_auto, PM_auto);
