
clc;
clear;
close all;

s = tf('s');

W = 3 / (0.1*s^3 + 0.01*s^2 + 0.1*s + 1);
Wz = feedback(W, 1);

disp('Полюса разомкнутой системы:');
disp(pole(W));

disp('Полюса замкнутой системы:');
disp(pole(Wz));

a = [0.1 0.01 0.1 1];
b = [0.1 0.01 0.1 4];

Hr = [a(2) a(4) 0;
      a(1) a(3) 0;
      0 a(2) a(4)];

Hz = [b(2) b(4) 0;
      b(1) b(3) 0;
      0 b(2) b(4)];

disp('Матрица Гурвица разомкнутой системы:');
disp(Hr);

disp('Миноры разомкнутой системы:');
disp([Hr(1,1), det(Hr(1:2,1:2)), det(Hr)]);

disp('Матрица Гурвица замкнутой системы:');
disp(Hz);

disp('Миноры замкнутой системы:');
disp([Hz(1,1), det(Hz(1:2,1:2)), det(Hz)]);

figure('Color', 'w');
nyquist(W);
grid on;
hold on;
plot(-1, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
xlim([-1.2 3.3]);
ylim([-1.8 1.8]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста');


w = linspace(0, 30, 3000);

D = 0.1*(1i*w).^3 + 0.01*(1i*w).^2 ...
    + 0.1*(1i*w) + 4;

figure('Color', 'w');

subplot(1, 2, 1);
plot(real(D), imag(D), 'b', 'LineWidth', 2);
hold on;
plot(4, 0, 'go', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Полная кривая Михайлова');

subplot(1, 2, 2);
plot(real(D), imag(D), 'b', 'LineWidth', 2);
hold on;
plot(4, 0, 'go', 'MarkerSize', 8, 'LineWidth', 2);
xlim([3.95 4.01]);
ylim([-0.1 0.1]);
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Начальный участок');

sgtitle('Кривая Михайлова замкнутой системы');

t = linspace(0, 5, 2000);

figure('Color', 'w');
step(W, t);
grid on;
xlabel('Время, с');
ylabel('Выходной сигнал');
title('Разомкнутая система');

figure('Color', 'w');
step(Wz, t);
grid on;
xlabel('Время, с');
ylabel('Выходной сигнал');
title('Замкнутая система');

figure('Color', 'w');
step(W, Wz, t);
grid on;
legend('Разомкнутая', 'Замкнутая', 'Location', 'best');
xlabel('Время, с');
ylabel('Выходной сигнал');
title('Сравнение переходных процессов');
