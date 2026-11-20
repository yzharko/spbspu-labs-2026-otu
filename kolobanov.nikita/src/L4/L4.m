clc;
clear;
close all;

s = tf('s');
FEEDBACK_GAIN = 1;
W_open = 10 / (s^5 + 3*s^4 + 2*s^3 + 2*s^2 + s + 1);
W_closed = feedback(W_open, FEEDBACK_GAIN);
CRITICAL_POINT = -1 / FEEDBACK_GAIN;

FREQUENCY_STEP = 0.01;
FREQUENCY_END = 20;
TIME_STEP = 0.02;
TIME_END = 25;
NYQUIST_ZOOM_LIMITS = [-2 2];
MARKER_SIZE = 8;
CRITICAL_MARKER_SIZE = 14;
LINE_WIDTH = 1.8;

[~, D_open] = tfdata(W_open, 'v');
[~, D_closed] = tfdata(W_closed, 'v');

disp('Характеристический полином разомкнутой системы:');
disp(D_open);
disp('Характеристический полином замкнутой системы:');
disp(D_closed);

H_open = [
    D_open(2) D_open(4) D_open(6) 0 0;
    D_open(1) D_open(3) D_open(5) 0 0;
    0 D_open(2) D_open(4) D_open(6) 0;
    0 D_open(1) D_open(3) D_open(5) 0;
    0 0 D_open(2) D_open(4) D_open(6)
    ];
H_closed = [
    D_closed(2) D_closed(4) D_closed(6) 0 0;
    D_closed(1) D_closed(3) D_closed(5) 0 0;
    0 D_closed(2) D_closed(4) D_closed(6) 0;
    0 D_closed(1) D_closed(3) D_closed(5) 0;
    0 0 D_closed(2) D_closed(4) D_closed(6)
    ];

minors_open = [H_open(1, 1), ...
    det(H_open(1:2, 1:2)), det(H_open(1:3, 1:3)), ...
    det(H_open(1:4, 1:4)), det(H_open)];
minors_closed = [H_closed(1, 1), ...
    det(H_closed(1:2, 1:2)), det(H_closed(1:3, 1:3)), ...
    det(H_closed(1:4, 1:4)), det(H_closed)];
disp('Миноры Гурвица разомкнутой системы:');
disp(minors_open);
disp('Миноры Гурвица замкнутой системы:');
disp(minors_closed);

figure;
nyquist(W_open);
grid on;
hold on;
plot(CRITICAL_POINT, 0, 'r+', 'MarkerSize', CRITICAL_MARKER_SIZE);
title('Годограф Найквиста разомкнутой системы');

figure;
nyquist(W_open);
grid on;
hold on;
plot(CRITICAL_POINT, 0, 'r+', 'MarkerSize', CRITICAL_MARKER_SIZE);
xlim(NYQUIST_ZOOM_LIMITS); ylim(NYQUIST_ZOOM_LIMITS);
title('Область критической точки');

w = 0:FREQUENCY_STEP:FREQUENCY_END;
Djw = polyval(D_closed, 1i*w);
figure;
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', LINE_WIDTH);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', ...
    'MarkerFaceColor', 'g', 'MarkerSize', MARKER_SIZE);
plot(real(Djw(end)), imag(Djw(end)), 'rs', ...
    'MarkerFaceColor', 'r', 'MarkerSize', MARKER_SIZE);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', ...
    'VerticalAlignment', 'bottom');
text(real(Djw(end)), imag(Djw(end)), ...
    sprintf('  \\omega = %g рад/с', FREQUENCY_END), ...
    'Color', 'r', 'VerticalAlignment', 'top');
xline(0, 'k:'); yline(0, 'k:');
grid on;
box on;
xlabel('Re D(j\omega)'); ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');
axis equal;
axis padded;

disp('Полюсы разомкнутой системы:');
disp(pole(W_open));
disp('Полюсы замкнутой системы:');
disp(pole(W_closed));

t = 0:TIME_STEP:TIME_END;
figure;
subplot(1, 2, 1);
step(W_open, t);
grid on;
title('Разомкнутая система');
xlabel('t'); ylabel('h(t)');

subplot(1, 2, 2);
step(W_closed, t);
grid on;
title('Замкнутая система');
xlabel('t'); ylabel('h(t)');

figure;
step(W_open, W_closed, t);
grid on;
title('Сравнение переходных процессов');
xlabel('t'); ylabel('h(t)');
legend('Разомкнутая система', 'Замкнутая система', 'Location', 'best');
