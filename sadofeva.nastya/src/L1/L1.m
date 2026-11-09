clear; clc; close all;

K0 = 1.9;
T  = 0.7;
T1 = 1.1;
T2 = 1.3;

T1_5 = T1;
if T1_5 < 2*T2
    T1_5 = 3.0;
    fprintf('Звено 5: T1 = %.1f < 2*T2 = %.1f, принято T1 = %.1f с\n', T1, 2*T2, T1_5);
end
T1_6 = T1;
if T1_6 >= 2*T2
    T1_6 = T2;
    fprintf('Звено 6: T1 = %.1f >= 2*T2 = %.1f, принято T1 = %.1f с\n', T1, 2*T2, T1_6);
end

K_set = K0 * [0.5 1 1.5];
time  = 0:0.01:15;

syms s t
syms K positive

U_step    = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

links = {
    'Интегрирующее',                K/s
    'Идеальное дифференцирующее',   K*s
    'Усилительное',                 K
    'Апериодическое 1-го порядка',  K/(T*s + 1)
    'Апериодическое 2-го порядка',  K/(T2^2*s^2 + T1_5*s + 1)
    'Колебательное',                K/(T2^2*s^2 + T1_6*s + 1)
    'Консервативное',               K/(T*s^2 + 1)
    'Инерционное интегрирующее',    K/(s*(T*s + 1))
    'Реальное дифференцирующее',    K*s/(T*s + 1)
    'Форсирующее 1-го порядка',     K*(T*s + 1)
    'Изодромное',                   K*(T*s + 1)/s
};
names = links(:, 1);
n = numel(names);

h = sym(zeros(n, 1));
w = sym(zeros(n, 1));

for i = 1:n
    W = links{i, 2};
    Y_step    = W * U_step;
    h(i)      = simplify(ilaplace(Y_step, s, t));
    Y_impulse = W * U_impulse;
    w(i)      = simplify(ilaplace(Y_impulse, s, t));

    fprintf('\n%2d. %s\n    W(s) = %s\n    h(t) = %s\n    w(t) = %s\n', ...
        i, names{i}, char(W), char(h(i)), char(w(i)));
end

plotAll(names, h, K, t, K_set, time, 'h(t)', ...
    'Переходные характеристики звеньев');
plotAll(names, w, K, t, K_set, time, 'w(t)', ...
    'Импульсные переходные характеристики звеньев');

p = tf('s');
W_check = cell(n, 1);
W_check{1}  = K0/p;
W_check{3}  = tf(K0);
W_check{4}  = K0/(T*p + 1);
W_check{5}  = K0/(T2^2*p^2 + T1_5*p + 1);
W_check{6}  = K0/(T2^2*p^2 + T1_6*p + 1);
W_check{7}  = K0/(T*p^2 + 1);
W_check{8}  = K0/(p*(T*p + 1));
W_check{9}  = K0*p/(T*p + 1);
W_check{11} = K0*(T*p + 1)/p;

eps_d   = 1e-3;
approx  = [2 10];
W_check{2}  = K0*p/(eps_d*p + 1)^2;
W_check{10} = K0*(T*p + 1)/(eps_d*p + 1)^2;

checkAll(names, h, W_check, approx, eps_d, K, K0, t, time, 'step', ...
    'Проверка переходных характеристик (K = 1.9)');
checkAll(names, w, W_check, approx, eps_d, K, K0, t, time, 'impulse', ...
    'Проверка импульсных характеристик (K = 1.9)');

T_list = T * [1/5 1/2 1 2 5];
time4  = 0:0.01:20;

figure('Name', 'Задание *: звено 4', 'Color', 'w');
hold on; grid on; box on;
fprintf('\nЗвено 4 (зона 5%%):\n   T, с    t_п, с    3T, с\n');
for Tv = T_list
    W4 = K0/(Tv*p + 1);
    y  = step(W4, time4);
    si = stepinfo(W4, 'SettlingTimeThreshold', 0.05);
    plot(time4, y, 'LineWidth', 1.4, 'DisplayName', ...
        sprintf('T = %.2f с, t_п = %.2f с', Tv, si.SettlingTime));
    plot(Tv, K0*(1 - exp(-1)), 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
    fprintf('  %5.2f   %6.2f    %5.2f\n', Tv, si.SettlingTime, 3*Tv);
end
yline(K0, 'k-', 'HandleVisibility', 'off');
yline(0.95*K0, 'k:', 'зона 5%', 'HandleVisibility', 'off');
title('Апериодическое звено 1-го порядка: влияние постоянной времени T');
xlabel('t, с'); ylabel('h(t)');
legend('Location', 'southeast');

T1_list = T1_6 * [1 1/2 1/4 1/10 0];
time6   = 0:0.02:100;
w0      = 1/T2;

figure('Name', 'Задание *: звено 6', 'Color', 'w');
hold on; grid on; box on;
fprintf(['\nЗвено 6 (T2 = %.1f с, зона 5%%):\n' ...
    '   T1, с    zeta    sigma_т, %%  sigma, %%   t_п_т, с   t_п, с\n'], T2);
for T1v = T1_list
    W6   = K0/(T2^2*p^2 + T1v*p + 1);
    y    = step(W6, time6);
    zeta = T1v/(2*T2);
    if T1v > 0
        si    = stepinfo(W6, 'SettlingTimeThreshold', 0.05);
        sig_t = 100*exp(-pi*zeta/sqrt(1 - zeta^2));
        ts_t  = 3/(zeta*w0);
        lbl = sprintf('T_1 = %.3f с, \\sigma = %.1f%%, t_п = %.1f с', ...
            T1v, si.Overshoot, si.SettlingTime);
        fprintf('  %6.3f   %5.3f   %8.1f   %7.1f   %8.1f   %6.1f\n', ...
            T1v, zeta, sig_t, si.Overshoot, ts_t, si.SettlingTime);
    else
        lbl = 'T_1 = 0 (консервативное), колебания не затухают';
        fprintf('  %6.3f   %5.3f   колебания не затухают (sigma = 100%%)\n', T1v, zeta);
    end
    plot(time6, y, 'LineWidth', 1.2, 'DisplayName', lbl);
end
yline(K0, 'k-', 'HandleVisibility', 'off');
yline(0.95*K0, 'k:', 'HandleVisibility', 'off');
yline(1.05*K0, 'k:', 'HandleVisibility', 'off');
ylim([0 5.2]);
title(sprintf('Колебательное звено: T_2 = %.1f с, уменьшение T_1', T2));
xlabel('t, с'); ylabel('h(t)');
legend('Location', 'northeast');

function plotAll(names, F, K, t, K_set, time, yname, figName)
figure('Name', figName, 'Color', 'w', 'Position', [40 40 1400 900]);
tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(names)
    nexttile; hold on; grid on; box on;
    [reg, sing] = splitSingular(F(i), t);
    for k = K_set
        y = evalOnGrid(subs(reg, K, k), t, time);
        plot(time, y, 'LineWidth', 1.3, 'DisplayName', sprintf('K = %.2f', k));
    end
    if ~isequal(sing, sym(0))
        markSingular(sing);
    end
    title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
    xlabel('t, с'); ylabel(yname);
    if i == 1
        legend('Location', 'northwest', 'FontSize', 8);
    end
end
sgtitle(figName, 'FontWeight', 'bold');
end

function checkAll(names, F, W_check, approx, eps_d, K, K0, t, time, mode, figName)
figure('Name', figName, 'Color', 'w', 'Position', [60 60 1400 900]);
tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(names)
    nexttile; hold on; grid on; box on;
    [reg, sing] = splitSingular(F(i), t);
    isApprox = ismember(i, approx);
    if isApprox
        tt = 0:eps_d/10:time(end);
    else
        tt = time;
    end
    y_an = evalOnGrid(subs(reg, K, K0), t, tt);
    if strcmp(mode, 'step')
        y_num = step(W_check{i}, tt);
    else
        y_num = impulse(W_check{i}, tt);
    end
    y_num = y_num(:)';
    plot(tt, y_an, 'LineWidth', 2, 'DisplayName', 'аналитически');
    plot(tt, y_num, '--', 'LineWidth', 1.5, 'DisplayName', mode);

    if isApprox
        t0   = 50*eps_d;
        far  = tt >= t0;
        err  = max(abs(y_num(far) - y_an(far)));
        ts = 0:eps_d/1000:t0;
        if strcmp(mode, 'step')
            y_s = step(W_check{i}, ts);
        else
            y_s = impulse(W_check{i}, ts);
        end
        S    = trapz(ts, y_s(:)' - evalOnGrid(subs(reg, K, K0), t, ts));
        c    = deltaCoef(sing, K, K0, t);
        ylim([min(y_an) - 3, max(y_an) + 3]);
        title({sprintf('%d. %s (приближение, \\epsilon = %g с)', i, names{i}, eps_d), ...
            sprintf('t \\geq %.2f: max|\\Delta| = %.1e;  S = %.3f, коэф. при \\delta = %.3f', ...
            t0, err, S, c)}, 'FontSize', 8);
    else
        err = max(abs(y_num - y_an));
        title(sprintf('%d. %s, max|\\Delta| = %.1e', i, names{i}, err), 'FontSize', 9);
    end
    xlabel('t, с'); ylabel('Выходной сигнал');
    if i == 1
        legend('Location', 'northwest', 'FontSize', 8);
    end
end
sgtitle(figName, 'FontWeight', 'bold');
end

function c = deltaCoef(sing, K, K0, t)
c = subs(sing, dirac(1, t), 0);
c = subs(c, dirac(t), 1);
c = double(subs(c, K, K0));
end

function [reg, sing] = splitSingular(expr, t)
e   = subs(expr, heaviside(t), 1);
reg = subs(e, dirac(1, t), 0);
reg = subs(reg, dirac(t), 0);
sing = simplify(e - reg);
end

function y = evalOnGrid(expr, t, time)
f = matlabFunction(expr, 'Vars', t);
y = f(time);
if isscalar(y)
    y = repmat(y, size(time));
end
end

function markSingular(sing)
yl = ylim;
ylim([min(yl(1), 0), max(yl(2), 1)*1.25]);
yl  = ylim;
top = 0.8*yl(2);
plot([0 0], [0 top], 'k-', 'LineWidth', 1.8, 'HandleVisibility', 'off');
plot(0, top, 'k^', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
lbl = char(sing);
lbl = regexprep(lbl, 'dirac\(1,\s*t\)|dirac\(t,\s*1\)', 'δ''(t)');
lbl = regexprep(lbl, 'dirac\(t\)', 'δ(t)');
text(0.08, 0.9, ['+ ' lbl], 'Units', 'normalized', 'Interpreter', 'none', ...
    'FontSize', 8, 'BackgroundColor', 'w');
end
