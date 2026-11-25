clc;
clear;
close all;

%% Исходная система
s = tf('s');

num = 10;
den = [2 3 3 0.5 0.5 1];
W   = tf(num, den);

Dcl = den + [zeros(1, numel(den) - 1) num];
Wcl = tf(num, Dcl);

%% Критерий Гурвица
[Hopen,  dOpen ] = hurwitz(den);
[Hclose, dClose] = hurwitz(Dcl);
disp('Миноры Гурвица разомкнутой системы:');  disp(dOpen);
disp('Миноры Гурвица замкнутой системы:');    disp(dClose);
fprintf('Разомкнутая система устойчива: %d\n', all(den > 0) && all(dOpen > 0));
fprintf('Замкнутая система устойчива:   %d\n', all(Dcl > 0) && all(dClose > 0));

disp('Полюса разомкнутой системы:'); disp(roots(den));
disp('Полюса замкнутой системы:');   disp(roots(Dcl));

%% Построение годографа Найквиста
figure('Color', 'w');
nyquist(W); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 3]); ylim([-2 2]);
title('Область критической точки');

%% Построение кривой Михайлова
[num_cl, den_cl] = tfdata(Wcl, 'v');
D = round(den_cl, 10);

w   = 0:0.001:3;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)),   imag(Djw(1)),   'go', 'MarkerFaceColor', 'g', ...
     'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', ...
     'MarkerSize', 8);

text(real(Djw(1)),   imag(Djw(1)),   '  \omega = 0',   'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = max', 'Color', 'r');

plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);

grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega) замкнутой системы');
axis equal;

%% Переходные процессы
t = 0:0.01:25;
figure('Color', 'w'); step(W, t);   grid on;
title('Переходный процесс разомкнутой системы');
figure('Color', 'w'); step(Wcl, t); grid on;
title('Переходный процесс замкнутой системы');

%% Прямая проверка устойчивости
p_cl = roots(Dcl);
is_stable_closed = all(real(p_cl) < 0);
fprintf('Прямая проверка: замкнутая система устойчива = %d\n', is_stable_closed);

p_open = roots(den);
is_stable_open = all(real(p_open) < 0);
fprintf('Прямая проверка: разомкнутая система устойчива = %d\n', is_stable_open);

fprintf('Критерий Рауса для замкнутой системы:\n');
raus_table = raus(Dcl);
disp(raus_table);

%% Вспомогательные функции
function [H, minors] = hurwitz(c)
    n = numel(c) - 1;
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            k = 2*j - i;
            if k >= 0 && k <= n
                H(i, j) = c(k + 1);
            end
        end
    end
    minors = arrayfun(@(m) det(H(1:m, 1:m)), 1:n);
end

function T = raus(c)
    c = c(:).';
    n = numel(c) - 1;
    m = ceil((n + 1) / 2);
    T = zeros(n + 1, m);

    row1 = c(1:2:end);
    T(1, 1:numel(row1)) = row1;

    row2 = c(2:2:end);
    T(2, 1:numel(row2)) = row2;

    for i = 3:n + 1
        for j = 1:m - 1
            T(i, j) = (T(i-1, 1) * T(i-2, j+1) - T(i-2, 1) * T(i-1, j+1)) / T(i-1, 1);
        end
    end
end
