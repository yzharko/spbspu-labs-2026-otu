clc;
clear; 
close all;
s = tf('s');
Go = 1 / (s + 1)^2;
W_cl_init = feedback(Go, 1);

info_init = stepinfo(W_cl_init);
e_ss_init = abs(1 - dcgain(W_cl_init)); %находит коэффициент передачи, т.е. установившуюся ошибку

disp('Полюса исходной замкнутой системы:');
disp(pole(W_cl_init));
fprintf('Время переходного процесса t_n = %.3f с\n', info_init.SettlingTime);
fprintf('Перерегулирование sigma = %.2f %%\n', info_init.Overshoot);
fprintf('Установившаяся ошибка e_ss = %.3f\n\n', e_ss_init);

%1-й метод З-Н (по S-образной кривой разгона в точке перегиба t=1)
L = 1 - (1 - 2*exp(-1)) / exp(-1); 
T = 1 / exp(-1);                  

Kp_zn1 = 1.2 * (T / L);
Ti_zn1 = 2 * L;
Td_zn1 = 0.5 * L;
W_pid_zn1 = pidstd(Kp_zn1, Ti_zn1, Td_zn1); %создаёт объект ПИД-регулятора
W_cl_zn1 = feedback(W_pid_zn1 * Go, 1);
info_zn1 = stepinfo(W_cl_zn1);

%2-й метод З-Н с малым вспомогательным запаздыванием
tau = 0.1;
Go_delay = Go * exp(-tau * s);
[Gm, ~, Wcg, ~] = margin(Go_delay); %находит запас устойчивости по амплитуде,фазе,частоте среза по фазе -180, частоте
K_cr = Gm;
T_cr = 2 * pi / Wcg; %переводит критическую частоту в период незатухающих колебаний

Kp_zn2 = 0.6 * K_cr;
Ti_zn2 = 0.5 * T_cr;
Td_zn2 = 0.125 * T_cr;
W_pid_zn2 = pidstd(Kp_zn2, Ti_zn2, Td_zn2, 10);
W_cl_zn2 = feedback(W_pid_zn2 * Go, 1);
info_zn2 = stepinfo(W_cl_zn2);

disp('Метод Зиглера - Николса');
fprintf('1-й метод (кривая разгона): Kp=%.3f, Ti=%.3f, Td=%.3f -> t_n=%.3f с, sigma=%.2f %%\n', ...
    Kp_zn1, Ti_zn1, Td_zn1, info_zn1.SettlingTime, info_zn1.Overshoot);
fprintf('2-й метод (tau=0.1 с): K_cr=%.3f, T_cr=%.3f -> Kp=%.3f, Ti=%.3f, Td=%.3f -> t_n=%.3f с, sigma=%.2f %%\n\n', ...
    K_cr, T_cr, Kp_zn2, Ti_zn2, Td_zn2, info_zn2.SettlingTime, info_zn2.Overshoot);

opts = pidtuneOptions('PhaseMargin', 85);
[W_pid_auto, info_tune] = pidtune(Go, 'PID', 20, opts);
W_cl_auto = feedback(W_pid_auto * Go, 1);
info_auto = stepinfo(W_cl_auto);

%расчет запасов устойчивости по фазе для всех методов
[~, Pm_zn1] = margin(W_pid_zn1 * Go);
[~, Pm_zn2] = margin(W_pid_zn2 * Go);

disp('Автоматическая настройка pidtune');
disp(W_pid_auto);
fprintf('t_n = %.3f с (треб. <= 0.5 с), sigma = %.2f %% (треб. <= 5 %%)\n', ...
    info_auto.SettlingTime, info_auto.Overshoot);
fprintf('Запасы по фазе: З-Н (1-й) = %.1f град, З-Н (2-й) = %.1f град, pidtune = %.1f град\n', ...
    Pm_zn1, Pm_zn2, info_tune.PhaseMargin);

figure('Name', 'Переходные характеристики');
step(W_cl_init, 'k--', W_cl_zn1, 'r-', W_cl_zn2, 'm-.', W_cl_auto, 'b-', 4);
grid on;
legend('Исходная без ПИД', 'Зиглер-Николс (1-й метод)', 'Зиглер-Николс (2-й метод)', 'pidtune', 'Location', 'best');
title('Сравнение переходных процессов');
xlabel('Время, с');
ylabel('Выходной сигнал');
