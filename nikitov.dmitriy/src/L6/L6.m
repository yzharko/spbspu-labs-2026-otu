clear;
clc;

s = tf('s');
G = 1 / ((s + 0.2)^2 * (s + 4));
H0 = feedback(G, 1);
timeLimit = 5;
overshootLimit = 8;

[Kcr, ~, wcr] = margin(G);
Tcr = 2 * pi / wcr;
Kp = 0.6 * Kcr;
Ti = 0.5 * Tcr;
Td = 0.125 * Tcr;
Czn = pid(Kp, Kp / Ti, Kp * Td);
Lzn = Czn * G;
Hzn = feedback(Lzn, 1);

wc = 1.8;
Cauto = pidtune(G, 'PID', wc);
Lauto = Cauto * G;
Hauto = feedback(Lauto, 1);

fprintf('Kкр = %.4f, wкр = %.4f рад/с, Tкр = %.4f с\n', Kcr, wcr, Tcr);
fprintf('Зиглер–Николс: Kp = %.4f, Ti = %.4f с, Td = %.4f с\n', Kp, Ti, Td);
fprintf('Зиглер–Николс: Ki = %.4f, Kd = %.4f\n', Czn.Ki, Czn.Kd);
fprintf('pidtune: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Cauto.Kp, Cauto.Ki, Cauto.Kd);

names = {'Без регулятора', 'Зиглер–Николс', 'pidtune'};
loops = {G, Lzn, Lauto};
systems = {H0, Hzn, Hauto};
for k = 1:numel(systems)
    H = systems{k};
    info = stepinfo(H, 'SettlingTimeThreshold', 0.02);
    [gm, pm] = margin(loops{k});
    steadyError = 1 - dcgain(H);
    meetsLimits = info.SettlingTime <= timeLimit && info.Overshoot <= overshootLimit;
    fprintf('\n%s\n', names{k});
    disp(pole(H));
    fprintf('Устойчива: %d. Условия по tп и sigma выполнены: %d\n', isstable(H), meetsLimits);
    fprintf('tп = %.4f с, sigma = %.2f %%, ошибка = %.4f\n', ...
        info.SettlingTime, info.Overshoot, steadyError);
    fprintf('Запас по усилению: %.2f дБ, по фазе: %.2f град.\n', 20 * log10(gm), pm);
end

t = 0:0.02:30;
Hcr = feedback(Kcr * G, 1);

figure('Name', 'Исходная система и граница устойчивости', 'Position', [100 150 900 350]);
tiledlayout(1, 2);
nexttile;
plot(t, step(H0, t), 'LineWidth', 1.5);
grid on;
title('Без регулятора');
xlabel('Время, с');
ylabel('Выход');
nexttile;
plot(t, step(Hcr, t), 'LineWidth', 1.5);
grid on;
title('П-регулятор при критическом усилении');
xlabel('Время, с');
ylabel('Выход');

figure('Name', 'Сравнение ПИД-регуляторов', 'Position', [130 170 900 350]);
plot(t, step(Hzn, t), t, step(Hauto, t), 'LineWidth', 1.5);
hold on;
yline(1, ':', 'HandleVisibility', 'off');
grid on;
title('Переходные характеристики с ПИД-регуляторами');
xlabel('Время, с');
ylabel('Выход');
legend(names(2:3), 'Location', 'northeast');
