clc;
clear;
close all;

s = tf('s');
Q = 2*s^5 + 3*s^4 + 3*s^3 + 0.5*s^2 + 0.5*s + 1;
W_open = 10 / Q;
W_closed = feedback(W_open, 1);

w_nyquist = 0:0.005:20;
[re, im] = nyquist(W_open, w_nyquist);
re = squeeze(re);
im = squeeze(im);
re_full = [flipud(re(2:end)); re];
im_full = [-flipud(im(2:end)); im];
figure('Color', 'w');
plot(re_full, im_full, 'b-', 'LineWidth', 1.5);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
plot(re_full, im_full, 'b-', 'LineWidth', 1.5); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

[~, den] = tfdata(W_closed, 'v');
D = round(den, 10);

w = 0:0.005:1.8;
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

format short;
disp('Полюсы разомкнутой системы:');
disp(pole(W_open));
disp('Полюсы замкнутой системы:');
disp(pole(W_closed));

t = 0:0.02:14;
figure('Color', 'w');
subplot(2, 1, 1);
step(W_open, t);
grid on;
title('Переходный процесс разомкнутой системы');

subplot(2, 1, 2);
step(W_closed, t);
grid on;
title('Переходный процесс замкнутой системы');

