clc;
clear;
close all;

s = tf('s');
W_open = 10 / (2*s^5 + 3*s^4 + 3*s^3 + 0.5*s^2 + 0.5*s + 1);
W_closed = feedback(W_open, 1);

fig1 = figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');
saveas(fig1, 'Nyquist_full.png');

fig2 = figure('Color', 'w');
nyquist(W_open); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-5 5]); ylim([-5 5]);
title('Область критической точки');
saveas(fig2, 'Nyquist_zoom.png');

[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);
w = 0:0.02:20;
Djw = polyval(D, 1i*w);

fig3 = figure('Color', 'w');

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
saveas(fig3, 'Mikhailov.png');

disp('Полюса ЗАМКНУТОЙ системы:');
disp(pole(W_closed));
disp('Полюса РАЗОМКНУТОЙ системы:');
disp(pole(W_open));

fig4 = figure('Color', 'w');
step(W_open, 'b', W_closed, 'r');
grid on;
legend('Разомкнутая система', 'Замкнутая система');
title('Переходные процессы при ступенчатом воздействии');
xlabel('Время (с)');
ylabel('Амплитуда');
saveas(fig4, 'Step_response.png');


