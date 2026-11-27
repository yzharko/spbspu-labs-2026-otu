clc;
clear;
close all;

s = tf('s');
G = 4/(s*(s+0.5)*(s+8));
G0 = feedback(G, 1);
fprintf('Полюса исходной замкнутой системы:\n');
disp(pole(G0));
fprintf('Устойчивость: %d\n', isstable(G0));
[Gm, Pm, Wcg, Wcp] = margin(G);
Ku = Gm;
Tu = 2*pi/Wcg;
Kp = 0.6*Ku;
Ti = 0.5*Tu;
Td = 0.125*Tu;
Czn = pid(Kp, Kp/Ti, Kp*Td);
Gzn = feedback(Czn*G, 1);
Cauto = pidtune(G, 'PID', 5);
Gauto = feedback(Cauto*G, 1);
[Gm_zn, Pm_zn] = margin(Czn*G);
[Gm_auto, Pm_auto] = margin(Cauto*G);
if Gm_auto == 0
    Gm_auto = Inf;
end
t = 0:0.01:20;

fprintf('Ku = %.4f, Tu = %.4f\n', Ku, Tu);
fprintf('Зиглер-Николс: Kp = %.4f, Ti = %.4f, Td = %.4f\n', Kp, Ti, Td);
fprintf('Автонастройка: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Cauto.Kp, Cauto.Ki, Cauto.Kd);
fprintf('Исходная ошибка: %.4f\n', 1 - dcgain(G0));
fprintf('Запасы: ЗН Pm = %.2f, Gm = %.4f; авто Pm = %.2f, Gm = %.4f\n', Pm_zn, Gm_zn, Pm_auto, Gm_auto);
fprintf('Ошибки: ЗН %.4f, авто %.4f\n', 1 - dcgain(Gzn), 1 - dcgain(Gauto));
disp('Исходная система:');
disp(stepinfo(G0));
disp('Система с Зиглером-Николсом:');
disp(stepinfo(Gzn));
disp('Система с автонастройкой:');
disp(stepinfo(Gauto));

figure;
step(G0, t);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Исходная система');

figure;
step(feedback(Ku*G, 1), 0:0.01:10);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Колебания при K = K_{кр}');

figure;
subplot(1,2,1);
[y0, t0] = step(G0, t);
[yz, tz] = step(Gzn, t);
[ya, ta] = step(Gauto, t);
plot(t0, y0, 'b', tz, yz, 'r', ta, ya, 'g', 'LineWidth', 1.5);
grid on;
xlabel('t, с');
ylabel('y(t)');
legend('исходная', 'Зиглер-Николс', 'автонастройка');
subplot(1,2,2);
plot(t0, y0, 'b', tz, yz, 'r', ta, ya, 'g', 'LineWidth', 1.5);
grid on;
xlim([0 3]);
ylim([0 1.8]);
xlabel('t, с');
ylabel('y(t)');
