clc;
clear;
close all;

s = tf('s');
W_open   = 1 / (0.1*s^3 + 0.1*s^2 + s + 1);
W_closed = feedback(W_open, 1);

[~, D_open]   = tfdata(W_open, 'v');
[~, D_closed] = tfdata(W_closed, 'v');
D_open   = round(D_open, 10);
D_closed = round(D_closed, 10);


disp('Миноры Гурвица разомкнутой системы:');
disp(hurwitz_minors(D_open));
disp('Миноры Гурвица замкнутой системы:');
disp(hurwitz_minors(D_closed));

disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));


figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-4 3]); ylim([-4 4]);
title('Годограф Найквиста разомкнутой системы W(s)');


w = 0:0.01:10;
Djw = polyval(D_closed, 1i*w);

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


t = 0:0.01:20;

figure('Color', 'w');
subplot(2, 1, 1);
step(W_open, t); grid on;
title('Переходный процесс разомкнутой системы');
subplot(2, 1, 2);
step(W_closed, t); grid on;
title('Переходный процесс замкнутой системы');

function minors = hurwitz_minors(a)
    H = [a(2), a(4), 0;
         a(1), a(3), 0;
         0,    a(2), a(4)];
    minors = [H(1, 1), det(H(1:2, 1:2)), det(H)];
end
