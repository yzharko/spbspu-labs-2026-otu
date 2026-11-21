clear; 
clc; 
close all;

s = tf('s');
Go = 2 / ((s + 0.2) * (s + 0.8));
sigma_req = 8;
t_p_req = 5.0;

G_closed = feedback(Go,1);
S = stepinfo(G_closed);
fprintf('Полюса: '); fprintf('%.3f  ', pole(G_closed)); fprintf('\n');
fprintf('Устойчивость: %d\n', isstable(G_closed));
fprintf('t_п = %.4f c, sigma = %.2f %%, ошибка = %.4f\n\n', ...
        S.SettlingTime, S.Overshoot, dcgain(feedback(1,Go)));

tau = 0.1;
Go_delay = Go * exp(-tau * s);
[Kcr, ~, Wcg, ~] = margin(Go_delay);
Tcr = 2*pi/Wcg;

Kp_second_method = 0.6*Kcr;  
Ti_second_method = 0.5*Tcr;
Td_second_method = 0.125*Tcr;

fprintf('K_кр = %.4f, T_кр = %.4f c\n', Kcr, Tcr);
fprintf('Kp = %.4f, Ti = %.4f c, Td = %.4f c\n\n', ...
    Kp_second_method, Ti_second_method, Td_second_method);

Kp = Kp_second_method;
Ti = Ti_second_method;
Td = Td_second_method;

C = pid(Kp, Kp/Ti, Kp*Td);
Sp = stepinfo(feedback(C*Go,1));
fprintf('C(s) = %.4f + %.4f/s + %.4f*s\n', Kp, Kp/Ti, Kp*Td);
fprintf('t_п = %.4f c (Соблюдение требования: %d), sigma = %.2f %% (Соблюдение требования: %d), ошибка = %.4f\n\n', ...
        Sp.SettlingTime, Sp.SettlingTime <= t_p_req, ...
        Sp.Overshoot,    Sp.Overshoot <= sigma_req, ...
        dcgain(feedback(1,C*Go)));

opts = pidtuneOptions('PhaseMargin', 60);
Ca = pidtune(Go, 'PID', 2, opts);
Sa = stepinfo(feedback(Ca*Go,1));
fprintf('Автоматическая настройка: Kp=%.4f Ki=%.4f Kd=%.4f | t_п=%.4f c, sigma=%.2f %%\n', ...
        Ca.Kp, Ca.Ki, Ca.Kd, Sa.SettlingTime, Sa.Overshoot);

figure('Name', 'Переходные характеристики');
step(G_closed, 'k--', feedback(C*Go,1), 'm-.', feedback(Ca*Go,1), 'b-', 15);
grid on;
legend('Исходная без ПИД', 'Зиглер-Николс', 'pidtune', 'Location', 'best');
title('Сравнение переходных процессов');
xlabel('Время, с');
ylabel('Выходной сигнал');
saveas(gcf, 'step_response.png');

