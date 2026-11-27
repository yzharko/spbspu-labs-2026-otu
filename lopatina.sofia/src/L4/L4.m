clc; clear; close all;

s = tf('s');
W_open  = 10 / (s^5 + 3*s^4 + 2*s^3 + 2*s^2 + s + 1);
W_close = 10 / (s^5 + 3*s^4 + 2*s^3 + 2*s^2 + s + 11);

%% Годограф Найквиста разомкнутой системы
figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы');

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

%% Кривая Михайлова замкнутой системы
[num, den] = tfdata(W_close, 'v');
D = round(den, 10);

w = 0:0.02:20;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

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
time = 0:0.02:200;
[y_op, t_op] = step(W_open, time);
[y_cl, t_cl] = step(W_close, time);

figure('Color', 'w');
subplot(1,2,1);
plot(t_op, y_op, 'b', 'LineWidth', 2); grid on;
title('Разомкнутая система'); xlabel('t, с'); ylabel('y(t)');
ylim([-max(abs(y_op))*1.1 max(abs(y_op))*1.1]);

subplot(1,2,2);
plot(t_cl, y_cl, 'r', 'LineWidth', 2); grid on;
title('Замкнутая система'); xlabel('t, с'); ylabel('y(t)');
ylim([-max(abs(y_cl))*1.1 max(abs(y_cl))*1.1]);

sgtitle('Прямая проверка устойчивости', 'FontSize', 14);
