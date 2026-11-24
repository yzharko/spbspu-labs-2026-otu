clear;
clc;
close all;

s  = tf('s');
Go = 1/(s*(s+2)*(s+10));
sigma_req = 10;
t_p_req   = 1;

%% 1. Оценка исходной замкнутой системы
G_closed = feedback(Go, 1);
S0 = stepinfo(G_closed);

fprintf('\nИсходная замкнутая система\n');
fprintf('Полюса:       %s\n', sprintf('%.3f  ', pole(G_closed)));
fprintf('Устойчивость: %d\n', isstable(G_closed));
print_quality(S0, dcgain(feedback(1, Go)), t_p_req, sigma_req);

%% 2. Настройка ПИД методом Зиглера-Николса
[Kcr, ~, Wcg, ~] = margin(Go);
Tcr = 2*pi/Wcg;
Kp = 0.6*Kcr;
Ti = 0.5*Tcr;
Td = 0.125*Tcr;
C_zn = pid(Kp, Kp/Ti, Kp*Td);

fprintf('\nРасчёт по методу Зиглера-Николса\n');
fprintf('K_кр = %.4f, T_кр = %.4f с\n', Kcr, Tcr);

r_zn = evaluate('Зиглер-Николс', C_zn, Go);
print_result(r_zn, t_p_req, sigma_req);

figure;
step(r_zn.T);
grid on;
title('Переходный процесс с ПИД-регулятором (Зиглер-Николс)');

%% 3. Автоматическая настройка pidtune по умолчанию
C_auto = pidtune(Go, 'pid');

r_auto = evaluate('pidtune по умолчанию', C_auto, Go);
print_result(r_auto, t_p_req, sigma_req);

%% 4. Автоматическая настройка pidtune с параметрами
pidtune_wc = 4;
pidtune_pm = 70;
pidtune_focus = 'reference-tracking';

opt = pidtuneOptions('PhaseMargin', pidtune_pm, 'DesignFocus', pidtune_focus);
C_tuned = pidtune(Go, 'pid', pidtune_wc, opt);

name_tuned = sprintf('pidtune: wc = %g рад/с, PM = %g град, %s', ...
                     pidtune_wc, pidtune_pm, pidtune_focus);
r_tuned = evaluate(name_tuned, C_tuned, Go);
print_result(r_tuned, t_p_req, sigma_req);

%% 5. Сравнение переходных процессов
figure;
step(r_zn.T, 'b', r_auto.T, 'r--', r_tuned.T, 'm');
grid on;
legend(r_zn.name, r_auto.name, r_tuned.name, 'Location', 'northeast');
title('Переходные процессы: сравнение настроек');

function r = evaluate(name, C, Go)
    L = C * Go;
    r.name = name;
    r.C = C;
    r.T = feedback(L, 1);
    r.S = stepinfo(r.T);
    r.err = dcgain(feedback(1, L));
    [r.Gm, r.Pm, ~, r.Wc] = margin(L);
end

function print_pid(C)
    fprintf('Регулятор: Kp = %.4f, Ki = %.4f, Kd = %.4f (Ti = %.4f с, Td = %.4f с)\n', ...
            C.Kp, C.Ki, C.Kd, C.Kp/C.Ki, C.Kd/C.Kp);
end

function print_quality(S, err, t_req, sigma_req)
    ok = {'нет', 'да'};
    fprintf('t_п   = %.4f с (требование выполнено: %s)\n', ...
            S.SettlingTime, ok{(S.SettlingTime <= t_req) + 1});
    fprintf('sigma = %.2f %% (требование выполнено: %s)\n', ...
            S.Overshoot, ok{(S.Overshoot <= sigma_req) + 1});
    fprintf('Ошибка = %.4f\n', err);
end

function print_margins(Gm, Pm, Wc)
    fprintf('Запас устойчивости: Gm = %.2f дБ, Pm = %.2f град, Wc = %.3f рад/с\n', ...
            20*log10(Gm), Pm, Wc);
end

function print_result(r, t_req, sigma_req)
    fprintf('\n%s\n', r.name);
    print_pid(r.C);
    print_quality(r.S, r.err, t_req, sigma_req);
    print_margins(r.Gm, r.Pm, r.Wc);
end
