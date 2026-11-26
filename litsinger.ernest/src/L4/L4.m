clc;
clear;
close all;

K   = 10;
den = [1 3 2 2 1 1];

fontSize  = 14;
lineWidth = 2;
figPos    = [100 100 1000 700];

W_open   = tf(K, den);
W_closed = feedback(W_open, 1);

[~, D_open]   = tfdata(W_open, 'v');
[~, D_closed] = tfdata(W_closed, 'v');
D_open   = round(D_open, 10);
D_closed = round(D_closed, 10);

[H_open, minors_open]     = hurwitz(D_open);
[H_closed, minors_closed] = hurwitz(D_closed);

disp('Матрица Гурвица разомкнутой системы:');
disp(H_open);
disp('Главные миноры разомкнутой системы:');
disp(minors_open);
disp('Матрица Гурвица замкнутой системы:');
disp(H_closed);
disp('Главные миноры замкнутой системы:');
disp(minors_closed);

disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));

critPoint = [-1 0];
zoomX = [-3 2];
zoomY = [-2 2];

figure('Color', 'w', 'Position', figPos);
nyquist(W_open);
grid on;
hold on;
plot(critPoint(1), critPoint(2), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');
beautify(gcf, fontSize, lineWidth);

figure('Color', 'w', 'Position', figPos);
nyquist(W_open);
grid on;
hold on;
plot(critPoint(1), critPoint(2), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim(zoomX);
ylim(zoomY);
title('Область критической точки');
beautify(gcf, fontSize, lineWidth);

w_mikh = 0:0.001:1.5;
Djw = polyval(D_closed, 1i*w_mikh);

figure('Color', 'w', 'Position', figPos);
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', 'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = max', 'Color', 'r');
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([0 max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
grid on;
box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы D(j\omega)');
beautify(gcf, fontSize, lineWidth);

t_open   = 0:0.01:40;
t_closed = 0:0.01:15;

figure('Color', 'w', 'Position', figPos);
step(W_open, t_open);
grid on;
title('Переходный процесс разомкнутой системы');
beautify(gcf, fontSize, lineWidth);

figure('Color', 'w', 'Position', figPos);
step(W_closed, t_closed);
grid on;
title('Переходный процесс замкнутой системы');
beautify(gcf, fontSize, lineWidth);

figure('Color', 'w', 'Position', figPos);
step(W_open, W_closed, t_closed);
grid on;
legend('разомкнутая', 'замкнутая', 'Location', 'northwest');
title('Сравнение переходных процессов');
beautify(gcf, fontSize, lineWidth);

function beautify(fig, fontSize, lineWidth)
    set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
    set(findall(fig, 'Type', 'line'), 'LineWidth', lineWidth);
end

function [H, minors] = hurwitz(a)
    n = numel(a) - 1;
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            k = 2*j - i;
            if k >= 0 && k <= n
                H(i, j) = a(k + 1);
            end
        end
    end
    minors = zeros(1, n);
    for k = 1:n
        minors(k) = det(H(1:k, 1:k));
    end
end
