clc;
clear;
close all;

%% 1. Исходная замкнутая система, вариант 2
s = tf('s');
Go = 1/((s + 0.2)^2*(s + 4));
timeLimit = 5;
overshootLimit = 8;
band = 0.02;

%% 2. Второй метод Зиглера - Николса
[~, a] = tfdata(Go, 'v');
wcr = sqrt(a(3));
Kcr = a(2)*a(3) - a(4);
Tcr = 2*pi/wcr;
criticalClosed = feedback(Kcr*Go, 1);

Kp = 0.6*Kcr;
Ti = 0.5*Tcr;
Td = 0.125*Tcr;
Czn = pid(Kp, Kp/Ti, Kp*Td);

%% 3. Автоматическая настройка
Cauto = pidtune(Go, 'PID');
wcTarget = 1.6;
pmTarget = 70;
Cfinal = pidtune(Go, 'PID', wcTarget, pidtuneOptions('PhaseMargin', pmTarget));

%% 4. Сравнение систем
openSystems = {Go, Czn*Go, Cauto*Go, Cfinal*Go};
names = {'Исходная', 'Зиглер - Николс', 'pidtune', 'pidtune, 1.6 рад/с и 70°'};
t = (0:0.002:45).';
y = zeros(numel(t), numel(openSystems));
steady = zeros(1, numel(openSystems));
errorFinal = steady;
settlingTime = steady;
overshoot = steady;
gainMargin = steady;
phaseMargin = steady;
stable = false(size(steady));

fprintf('Kкр = %.4f, ωкр = %.4f рад/с, Tкр = %.4f с.\n', Kcr, wcr, Tcr);
fprintf('Зиглер - Николс: Kp = %.4f, Ti = %.4f с, Td = %.4f с.\n', Kp, Ti, Td);
fprintf('Итоговый ПИД: Kp = %.4f, Ki = %.4f, Kd = %.4f.\n', Cfinal.Kp, Cfinal.Ki, Cfinal.Kd);
for k = 1:numel(openSystems)
    L = openSystems{k};
    closed = feedback(L, 1);
    stable(k) = isstable(closed);
    steady(k) = dcgain(closed);
    errorFinal(k) = dcgain(feedback(1, L));
    y(:,k) = step(closed, t);
    stats = stepinfo(y(:,k), t, steady(k), 'SettlingTimeThreshold', band);
    settlingTime(k) = stats.SettlingTime;
    overshoot(k) = stats.Overshoot;
    [gainMargin(k), phaseMargin(k)] = margin(L);
    fprintf(['%s: уст. %d, tп %.3f с, σ %.2f %%, e∞ %.4f, ' ...
         'запас %.2f дБ / %.2f°.\n'], names{k}, stable(k), settlingTime(k), ...
         overshoot(k), errorFinal(k), 20*log10(gainMargin(k)), phaseMargin(k));
end
fprintf('Итоговое соответствие: tп ≤ %.1f с - %d, σ ≤ %.0f %% - %d.\n', ...
    timeLimit, settlingTime(4) <= timeLimit, overshootLimit, overshoot(4) <= overshootLimit);

%% 5. Переходные характеристики
figures = gobjects(3,1);
figures(1) = figure('Theme','light','Name','Исходная система', 'Position',[70 70 780 435]);
plot(t, y(:,1), 'LineWidth',1.5);
yline(steady(1), ':', 'y_\infty');
grid on;
xlim([0 30]);
xlabel('t, с'); ylabel('y(t)');
title('Переходная характеристика исходной замкнутой системы');

figures(2) = figure('Theme','light','Name','Критический режим', 'Position',[90 90 780 435]);
tCritical = (0:0.01:25).';
plot(tCritical, step(criticalClosed, tCritical), 'LineWidth',1.5);
grid on;
xlabel('t, с'); ylabel('y(t)');
title(sprintf('Пропорциональный регулятор при K_{кр} = %.3f', Kcr));

figures(3) = figure('Theme','light','Name','Сравнение ПИД-регуляторов', 'Position',[110 110 780 455]);
plot(t, y(:,2:4), 'LineWidth',1.5);
hold on;
yline(1, '-.', 'HandleVisibility','off');
yline(1+overshootLimit/100, ':', 'Предел 8 %', 'HandleVisibility','off');
xline(timeLimit, ':', 'Предел 5 с', 'HandleVisibility','off');
grid on;
xlim([0 30]);
xlabel('t, с'); ylabel('y(t)');
title('Переходные характеристики с ПИД-регуляторами');
legend(names(2:4), 'Location','northeast');
