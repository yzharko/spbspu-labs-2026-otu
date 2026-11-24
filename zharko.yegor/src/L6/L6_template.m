clear; 
clc; 
close all;

Go = ???;
sigma_req = 7;
t_p_req = 7;

%% 1. Оценка исходной ЗАМКНУТОЙ системы
G_closed = feedback(Go,1);
S = stepinfo(G_closed);
fprintf('Полюса: '); fprintf('%.3f  ', pole(G_closed)); fprintf('\n');
fprintf('Устойчивость: %d\n', isstable(G_closed));
fprintf('t_п = %.4f c, sigma = %.2f %%, ошибка = %.4f\n\n', ...
        S.SettlingTime, S.Overshoot, dcgain(feedback(1,Go)));

%% 2. Настройка ПИД методом Зиглера-Николса
% Метод 1 (1, 4, 6, 8, 10, 12 варианты)
L = ???;
T = ???;
Kp_first_method = 1.2*T/L;   
Ti_first_method = 2*L;   
Td_first_method = 0.5*L;
fprintf('Метод 1: L = %.4f c, T = %.4f c\n', L, T);
fprintf('         Kp = %.4f, Ti = %.4f c, Td = %.4f c\n\n', Kp1, Ti1, Td1);

% Метод 2 (все остальные варианты)
[Kcr, ~, Wcg, ~] = margin(???);
Tcr = 2*pi/Wcg;
Kp_second_method = 0.6*Kcr;  
Ti_second_method = 0.5*Tcr;
Td_second_method = 0.125*Tcr;
fprintf('Метод 2: K_кр = %.4f, T_кр = %.4f c\n', Kcr, Tcr);
fprintf('         Kp = %.4f, Ti = %.4f c, Td = %.4f c\n', Kp2, Ti2, Td2);

%% 3. Оценка системы с ПИД регулятором
% т.к. в рамках работы будет использован один из методов выше
% под Kp, Ti и Td ниже мысленно принимаем соответсвующие переменные
% с припиской "..._first_method" / "..._second_method"
C = pid(Kp, Kp/Ti, Kp*Td);
Sp = stepinfo(feedback(C*Go,1));
fprintf('C(s) = %.4f + %.4f/s + %.4f*s\n', Kp, Kp/Ti, Kp*Td);
fprintf('t_п = %.4f c (Соблюдение требования: %d), sigma = %.2f %% (Соблюдение требования: %d), ошибка = %.4f\n\n', ...
        Sp.SettlingTime, Sp.SettlingTime <= sigma_req, ...
        Sp.Overshoot,    Sp.Overshoot <= t_p_req, ...
        dcgain(feedback(1,C*Go)));

%% 4. Автоматическая настройка с помощью pidtune
Ca = pidtune(???);
Sa = stepinfo(feedback(???));
fprintf('Автоматическая настройка: Kp=%.4f Ki=%.4f Kd=%.4f | t_п=%.4f c, sigma=%.2f %%\n', ...
        Ca.Kp, Ca.Ki, Ca.Kd, Sa.SettlingTime, Sa.Overshoot);
