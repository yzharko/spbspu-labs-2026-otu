clc;
clear;
close all;

%% Вариант 3: W(s) = 1 / (0.1*s^3 + 0.1*s^2 + s + 1)
s = tf('s');
W_open = 1 / (0.1*s^3 + 0.1*s^2 + s + 1);
W_closed = feedback(W_open, 1);

%% Характеристические уравнения
[~, den_open]   = tfdata(W_open, 'v'); % D_откр(s) = 0.1s^3 + 0.1s^2 + s + 1
[~, den_closed] = tfdata(W_closed, 'v'); % D_закр(s) = 0.1s^3 + 0.1s^2 + s + 2

disp('Характеристический полином разомкнутой системы:');
disp(poly2str(den_open, 's'));
disp('Характеристический полином замкнутой системы:');
disp(poly2str(den_closed, 's'));

%% Критерий Гурвица (система 3-го порядка)
names_h = {'разомкнутой', 'замкнутой'};
dens = {den_open, den_closed};
for k = 1:2
    a = dens{k};
    a0 = a(1); a1 = a(2); a2 = a(3); a3 = a(4);
    H = [a1 a3 0;
        a0 a2 0;
        0  a1 a3];
    d1 = H(1,1);
    d2 = det(H(1:2,1:2));
    d3 = det(H);
    fprintf('\nГурвиц, %s система: a_i = [%g %g %g %g]\n', names_h{k}, a0, a1, a2, a3);
    fprintf('  Delta1 = %g, Delta2 = %g, Delta3 = %g\n', d1, d2, d3);
    if all(a > 0) && d1 > 0 && abs(d2) < 1e-9
        disp('  Вывод: система на границе устойчивости (Delta2 = 0)');
    elseif all(a > 0) && d1 > 0 && d2 > 0 && d3 > 0
        disp('  Вывод: система устойчива');
    else
        disp('  Вывод: система неустойчива');
    end
end

%% Построение годографа Найквиста разомкнутой системы
figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_open); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

%% Построение кривой Михайлова замкнутой системы
D = round(den_closed, 10);

w = 0:0.01:10;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', ...
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
title('Кривая Михайлова D(j\omega)');
axis equal;

%% Прямая проверка устойчивости
disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутойы системы:');
disp(pole(W_closed));

t = 0:0.01:20;

figure('Color', 'w');
subplot(2, 1, 1);
step(W_open, t);
grid on;
title('Переходный процесс разомкнутой системы');
xlabel('t, c'); ylabel('y(t)');

subplot(2, 1, 2);
step(W_closed, t);
grid on;
title('Переходный процесс замкнутой системы');
xlabel('t, c'); ylabel('y(t)');