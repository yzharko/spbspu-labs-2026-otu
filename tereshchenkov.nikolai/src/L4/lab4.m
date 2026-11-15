clc;
clear;
close all;

s = tf('s');
W = 1 / (8*s^3 + 4*s^2 + 2*s + 1);
Wz = feedback(W, 1);

disp('Разомкнутая W(s)');
W
disp('Замкнутая Wz(s)');
Wz

fprintf('Полюса разомкнутой системы:\n');
disp(pole(W));
fprintf('Полюса замкнутой системы:\n');
disp(pole(Wz));

a = [8 4 2 1];
b = [8 4 2 2];

Ho = [a(2) a(4) 0;
      a(1) a(3) 0;
      0    a(2) a(4)];
Hc = [b(2) b(4) 0;
      b(1) b(3) 0;
      0    b(2) b(4)];

fprintf('Матрица Гурвица разомкнутой системы:\n');
disp(Ho);
d1o = a(2);
d2o = a(2)*a(3) - a(1)*a(4);
d3o = a(4)*d2o;
fprintf('Миноры разомкнутой: d1 = %g, d2 = %g, d3 = %g\n', d1o, d2o, d3o);

fprintf('Матрица Гурвица замкнутой системы:\n');
disp(Hc);
d1c = b(2);
d2c = b(2)*b(3) - b(1)*b(4);
d3c = b(4)*d2c;
fprintf('Миноры замкнутой: d1 = %g, d2 = %g, d3 = %g\n', d1c, d2c, d3c);

%% Годограф Найквиста разомкнутой системы
wLo = logspace(-3, log10(0.49), 500);
wHi = logspace(log10(0.51), 2, 500);

figure('Color', 'w');
nyquist(W, wLo); hold on;
nyquist(W, wHi);
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-6 4]); ylim([-5 5]);
grid on;
title('Годограф Найквиста разомкнутой системы');

figure('Color', 'w');
nyquist(W, wLo); hold on;
nyquist(W, wHi);
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
grid on;
title('Область критической точки');

%% Кривая Михайлова замкнутой системы
[~, den] = tfdata(Wz, 'v');
D = round(den, 10);

w = 0:0.01:2;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', 'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = 2', 'Color', 'r');
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(2, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(1, 0, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
plot(0, -sqrt(2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(2, 0.06, '  \omega = 0', 'Color', [0 0.5 0]);
text(1, 0.12, '  \omega = 0.5', 'Color', [0.5 0 0.5]);
xlim([-2.2 2.6]); ylim([-2 0.8]);
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Начальный участок кривой Михайлова');

%% Прямая проверка устойчивости
t = 0:0.01:40;

figure('Color', 'w');
step(W, t);
grid on;
title('Переходный процесс разомкнутой системы');

figure('Color', 'w');
step(Wz, t);
grid on;
title('Переходный процесс замкнутой системы');

figure('Color', 'w');
step(W, Wz, t);
grid on;
legend('разомкнутая', 'замкнутая', 'Location', 'best');
title('Сравнение переходных процессов');
