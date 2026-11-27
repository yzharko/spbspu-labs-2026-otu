clear;
clc;
close all;

den = conv(conv([1 0.2], [1 0.2]), [1 4]);
Go = tf(1, den);
sigma_req = 8;
t_p_req = 5;
band = 0.02;

%% Исходная замкнутая система
G_closed = feedback(Go, 1);
S = stepinfo(G_closed, 'SettlingTimeThreshold', band);
y0_inf = dcgain(G_closed);
e0 = dcgain(feedback(1, Go));
fprintf(['\nИсходная система: устойчива = %d, t_п = %.4f с, ', ...
    'sigma = %.2f %%, e(∞) = %.6f. Полюса:\n'], ...
    isstable(G_closed), S.SettlingTime, S.Overshoot, e0);
disp(pole(G_closed));

t0 = linspace(0, 1.2*max(t_p_req, S.SettlingTime), 4001).';
figure('Color', 'w', 'Name', 'Исходная замкнутая система');
plot(t0, step(G_closed, t0), 'Color', [0.20 0.40 0.67], 'LineWidth', 1.5);
hold on;
yline(1, 'k:', 'Задание');
yline(y0_inf, 'k--', sprintf('y(∞) = %.3f', y0_inf));
yline(y0_inf*(1 + band), ':', 'Color', [0.55 0.55 0.55]);
yline(y0_inf*(1 - band), ':', 'Color', [0.55 0.55 0.55]);
xline(S.SettlingTime, 'k:', sprintf('t_п = %.2f с', S.SettlingTime));
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Переходная характеристика исходной замкнутой системы');

%% Второй метод Зиглера - Николса
a2 = den(2);
a1 = den(3);
a0 = den(4);
omega_cr = sqrt(a1);
Kcr = a2*a1 - a0;
Tcr = 2*pi/omega_cr;

Kp = 0.6*Kcr;
Ti = 0.5*Tcr;
Td = 0.125*Tcr;
Ki = Kp/Ti;
Kd = Kp*Td;
C = pid(Kp, Ki, Kd);

fprintf('\nЗиглер - Николс: K_кр = %.6f, omega_кр = %.6f рад/с, T_кр = %.6f с\n', Kcr, omega_cr, Tcr);
fprintf('Kp = %.6f, Ti = %.6f с, Td = %.6f с, Ki = %.6f, Kd = %.6f\n', Kp, Ti, Td, Ki, Kd);

%% Система с ПИД-регулятором Зиглера - Николса
L_zn = tf(C)*Go;
G_zn = feedback(L_zn, 1);
Sp = stepinfo(G_zn, 'SettlingTimeThreshold', band);
e_zn = dcgain(feedback(1, L_zn));
pass_time_zn = Sp.SettlingTime <= t_p_req;
pass_os_zn = Sp.Overshoot <= sigma_req;

t_zn = linspace(0, 1.2*max(t_p_req, Sp.SettlingTime), 4001).';
figure('Color', 'w', 'Name', 'ПИД-регулятор Зиглера — Николса');
plot(t_zn, step(G_zn, t_zn), 'Color', [0.73 0.26 0.18], 'LineWidth', 1.6);
plotPidLimits(t_p_req, sigma_req, band);
title('Переходная характеристика с ПИД-регулятором Зиглера - Николса');

%% Автоматическая настройка ПИД и сравнение
Ca = pidtune(Go, 'PID');
L_auto = tf(Ca)*Go;
G_auto = feedback(L_auto, 1);
Sa = stepinfo(G_auto, 'SettlingTimeThreshold', band);
e_auto = dcgain(feedback(1, L_auto));
pass_time_auto = Sa.SettlingTime <= t_p_req;
pass_os_auto = Sa.Overshoot <= sigma_req;

fprintf('\npidtune: Kp = %.6f, Ki = %.6f, Kd = %.6f\n', Ca.Kp, Ca.Ki, Ca.Kd);

[gm_zn, pm_zn, ~, wc_zn] = margin(L_zn);
[gm_auto, pm_auto, ~, wc_auto] = margin(L_auto);
comparison = table({'Зиглер - Николс'; 'pidtune'}, ...
    [isstable(G_zn); isstable(G_auto)], ...
    [Sp.SettlingTime; Sa.SettlingTime], [Sp.Overshoot; Sa.Overshoot], ...
    [e_zn; e_auto], [pass_time_zn; pass_time_auto], [pass_os_zn; pass_os_auto], ...
    [gm_zn; gm_auto], 20*log10([gm_zn; gm_auto]), [pm_zn; pm_auto], [wc_zn; wc_auto], ...
    'VariableNames', {'Method', 'Stable', 'Tp_s', 'Overshoot_pct', 'Ess', ...
    'TimeOK', 'OvershootOK', 'GainMargin', 'GainMargin_dB', 'PhaseMargin_deg', 'Wc_rad_s'});
disp(comparison);

t = linspace(0, 1.2*max([t_p_req, Sp.SettlingTime, Sa.SettlingTime]), 4001).';
figure('Color', 'w', 'Name', 'Сравнение ПИД-регуляторов');
curves = plot(t, step(G_zn, t), '-', t, step(G_auto, t), '--', 'LineWidth', 1.6);
curves(1).Color = [0.73 0.26 0.18];
curves(2).Color = [0.20 0.40 0.67];
plotPidLimits(t_p_req, sigma_req, band);
title('Сравнение переходных характеристик с ПИД-регуляторами');
legend(curves, 'Зиглер - Николс', 'Автонастройка pidtune', 'Location', 'best');

function plotPidLimits(t_req, sigma_req, band)
    hold on;
    yline(1, 'k:', 'Задание', 'HandleVisibility', 'off');
    yline(1 + band, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    yline(1 - band, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
    yline(1 + sigma_req/100, 'k--', 'Предел перерегулирования', 'HandleVisibility', 'off');
    xline(t_req, 'k:', 'Предел t_п', 'HandleVisibility', 'off');
    grid on;
    xlabel('t, с');
    ylabel('y(t)');
end
