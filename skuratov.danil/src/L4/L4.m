clc;
clear;
close all;
s = tf('s');
W_open = 1 / (0.1*s^3 + 0.1*s^2 + s + 1);
W_closed = feedback(W_open, 1);
output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

[~, D_open] = tfdata(W_open, 'v');
[~, D_closed] = tfdata(W_closed, 'v');
D_open = round(D_open, 10);
D_closed = round(D_closed, 10);
H_open = [D_open(2) D_open(4) 0;
          D_open(1) D_open(3) 0;
          0 D_open(2) D_open(4)];
H_closed = [D_closed(2) D_closed(4) 0;
            D_closed(1) D_closed(3) 0;
            0 D_closed(2) D_closed(4)];
for j = 1:3
    delta_open(j) = det(H_open(1:j, 1:j));
    delta_closed(j) = det(H_closed(1:j, 1:j));
end
disp('Миноры Гурвица разомкнутой системы:');
disp(delta_open);
disp('Миноры Гурвица замкнутой системы:');
disp(delta_closed);
disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));

omega0 = sqrt(10);
omega1 = linspace(0, omega0 - 0.001, 1500);
omega2 = logspace(log10(omega0 + 0.001), 2, 1500);
figure('Color', 'w', 'Position', [100 100 900 650]);
nyquist(W_open, omega1);
hold on;
nyquist(W_open, omega2);
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-4 3]);
ylim([-5 5]);
grid on;
title('Годограф Найквиста разомкнутой системы');
print(gcf, fullfile(output_dir, 'nyquist.png'), '-dpng', '-r180');

omega = 0:0.01:5;
Djw = polyval(D_closed, 1i*omega);
figure('Color', 'w', 'Position', [100 100 900 650]);
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(2, 0, 'go', 'MarkerFaceColor', 'g');
plot(1, 0, 'mo', 'MarkerFaceColor', 'm');
plot(0, -sqrt(20), 'rs', 'MarkerFaceColor', 'r');
text(2, 0.3, '  \omega = 0');
text(1, -0.5, '  \omega = \surd10');
text(0, -sqrt(20), '  \omega = \surd20');
xline(0, 'k:');
yline(0, 'k:');
xlim([-0.6 2.4]);
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');
print(gcf, fullfile(output_dir, 'mikhailov.png'), '-dpng', '-r180');

t = 0:0.01:20;
figure('Color', 'w', 'Position', [100 100 900 550]);
step(W_open, t);
grid on;
title('Переходный процесс разомкнутой системы');
xlabel('Время');
ylabel('h(t)');
print(gcf, fullfile(output_dir, 'step_open.png'), '-dpng', '-r180');

figure('Color', 'w', 'Position', [100 100 900 550]);
step(W_closed, t);
grid on;
title('Переходный процесс замкнутой системы');
xlabel('Время');
ylabel('h(t)');
print(gcf, fullfile(output_dir, 'step_closed.png'), '-dpng', '-r180');
