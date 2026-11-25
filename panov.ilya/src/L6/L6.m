clear;
clc;
close all;

Go = tf(2, [1 0.2 0]);
sigma_req = 5;
t_p_req = 6;
band = 0.02;

G_closed = feedback(Go, 1);
S = stepinfo(G_closed, 'SettlingTimeThreshold', band);
error_initial = dcgain(feedback(1, Go));
fprintf('Исходная замкнутая система\n');
fprintf('Полюса:\n');
disp(pole(G_closed));
fprintf('Устойчивость: %d\n', isstable(G_closed));
fprintf('Время установления = %.4f с, перерегулирование = %.2f %%, ошибка = %.6f\n', S.SettlingTime, S.Overshoot, error_initial);

tau = 0.1;
Go_aux = Go;
Go_aux.InputDelay = tau;
[Kcr, ~, Wcg, ~] = margin(Go_aux);
if ~isfinite(Kcr) || ~isfinite(Wcg) || Kcr <= 0 || Wcg <= 0
    error('Не получены конечные положительные критические параметры.');
end
Tcr = 2*pi/Wcg;
Kp = 0.6*Kcr;
Ti = 0.5*Tcr;
Td = 0.125*Tcr;
fprintf('\nМетод Зиглера-Николса\n');
fprintf('Вспомогательное запаздывание для расчёта критических параметров = %.3f с\n', tau);
fprintf('Критическая частота = %.8f рад/с, Kкр = %.8f, Tкр = %.8f с\n', Wcg, Kcr, Tcr);
fprintf('Kp = %.8f, Ti = %.8f с, Td = %.8f с\n', Kp, Ti, Td);

C = pid(Kp, Kp/Ti, Kp*Td);
G_pid = feedback(C*Go, 1);
Sp = stepinfo(G_pid, 'SettlingTimeThreshold', band);
error_pid = dcgain(feedback(1, C*Go));
[Gm_pid, Pm_pid, ~, wc_pid] = margin(C*Go);
fprintf('Проверка регулятора на исходном объекте без вспомогательного запаздывания\n');
fprintf('C(s) = %.8f + %.8f/s + %.8f*s\n', C.Kp, C.Ki, C.Kd);
fprintf('Полюса замкнутой системы:\n');
disp(pole(G_pid));
fprintf('Устойчивость: %d\n', isstable(G_pid));
fprintf('Время установления = %.4f с, требование <= %.2f с: %d\n', Sp.SettlingTime, t_p_req, Sp.SettlingTime <= t_p_req);
fprintf('Перерегулирование = %.2f %%, требование <= %.2f %%: %d\n', Sp.Overshoot, sigma_req, Sp.Overshoot <= sigma_req);
fprintf('Установившаяся ошибка = %.6f\n', error_pid);
fprintf('Запас по фазе = %.4f град, запас по усилению = %.4f дБ, частота среза = %.4f рад/с\n', Pm_pid, 20*log10(Gm_pid), wc_pid);

Ca = pidtune(Go, 'PID');
G_auto = feedback(Ca*Go, 1);
Sa = stepinfo(G_auto, 'SettlingTimeThreshold', band);
error_auto = dcgain(feedback(1, Ca*Go));
[Gm_auto, Pm_auto, ~, wc_auto] = margin(Ca*Go);
fprintf('\nАвтоматическая настройка pidtune\n');
fprintf('Cauto(s) = %.8f + %.8f/s + %.8f*s\n', Ca.Kp, Ca.Ki, Ca.Kd);
fprintf('Полюса замкнутой системы:\n');
disp(pole(G_auto));
fprintf('Устойчивость: %d\n', isstable(G_auto));
fprintf('Время установления = %.4f с, требование <= %.2f с: %d\n', Sa.SettlingTime, t_p_req, Sa.SettlingTime <= t_p_req);
fprintf('Перерегулирование = %.2f %%, требование <= %.2f %%: %d\n', Sa.Overshoot, sigma_req, Sa.Overshoot <= sigma_req);
fprintf('Установившаяся ошибка = %.6f\n', error_auto);
fprintf('Запас по фазе = %.4f град, запас по усилению = %.4f дБ, частота среза = %.4f рад/с\n', Pm_auto, 20*log10(Gm_auto), wc_auto);

t = linspace(0, 1.2*max([S.SettlingTime Sp.SettlingTime Sa.SettlingTime t_p_req]), 10000);
y_initial = step(G_closed, t);
y_pid = step(G_pid, t);
y_auto = step(G_auto, t);

newFigure('Оценка исходной системы и настройка ПИД');
subplot(2,1,1);
plot(t, y_initial, 'Color', [0.65 0.18 0.04], 'LineWidth', 1.6);
hold on;
yline(1, ':k', 'HandleVisibility', 'off');
yline(1-band, '--k', 'HandleVisibility', 'off');
yline(1+band, '--k', 'HandleVisibility', 'off');
grid on;
title('Исходная замкнутая система');
xlabel('t, с');
ylabel('y(t)');
subplot(2,1,2);
plot(t, y_initial, '--', t, y_pid, '-', 'LineWidth', 1.6);
hold on;
yline(1, ':k', 'HandleVisibility', 'off');
yline(1-band, '--k', 'HandleVisibility', 'off');
yline(1+band, '--k', 'HandleVisibility', 'off');
grid on;
title('Сравнение с настройкой Зиглера-Николса');
xlabel('t, с');
ylabel('y(t)');
legend('Без регулятора', 'Зиглер-Николс', 'Location', 'best');

newFigure('Сравнение способов настройки ПИД');
subplot(2,1,1);
plot(t, y_pid, '-', t, y_auto, '--', 'LineWidth', 1.6);
hold on;
yline(1, ':k', 'HandleVisibility', 'off');
yline(1-band, ':k', 'HandleVisibility', 'off');
yline(1+band, ':k', 'HandleVisibility', 'off');
grid on;
title('Зиглер-Николс и автоматическая настройка');
xlabel('t, с');
ylabel('y(t)');
legend('Зиглер-Николс', 'pidtune', 'Location', 'best');
subplot(2,1,2);
plot(t, y_pid, '-', t, y_auto, '--', 'LineWidth', 1.6);
hold on;
yline(1, ':k', 'HandleVisibility', 'off');
yline(1+sigma_req/100, '--k', 'Предел перерегулирования', 'HandleVisibility', 'off');
xline(t_p_req, ':k', 'Требуемое время', 'HandleVisibility', 'off');
xlim([0 2*t_p_req]);
grid on;
title('Начальный участок переходных процессов');
xlabel('t, с');
ylabel('y(t)');
legend('Зиглер-Николс', 'pidtune', 'Location', 'best');

function newFigure(name)
    f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Units', 'normalized', 'Position', [0.10 0.08 0.80 0.80]);
    if isprop(f, 'Theme')
        f.Theme = 'light';
    end
    set(f, 'DefaultAxesColor', 'w', 'DefaultAxesXColor', 'k', 'DefaultAxesYColor', 'k', 'DefaultTextColor', 'k');
    set(f, 'DefaultLegendColor', 'w', 'DefaultLegendTextColor', 'k', 'DefaultAxesColorOrder', [0 0.32 0.65; 0.65 0.18 0.04]);
end
