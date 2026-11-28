clc;
clear;
close all;

s = tf('s');
W_opened = 2/(s^4 + 5*s^3 + 5*s^2 + 3*s + 1);
W_closed = feedback(W_opened, 1);

%% Исследование САУ, применяя критерий Гурвица.
[num_open, den_open] = tfdata(W_opened, 'v');
[num_closed, den_closed] = tfdata(W_closed, 'v');

function H = matrix(den)
    n = length(den) - 1;
    a = den;
    H = zeros(n);
    for i = 1:n
        for j = 1:n
            k = 2*i - j;  
            if k >= 0 && k <= n
                H(i,j) = a(k+1);
            end
        end
    end
end

H_open   = matrix(den_open);
H_closed = matrix(den_closed);

for k = 1:length(den_open)-1
    fprintf('D%d(open)   = %g\n', k, det(H_open(1:k,1:k)));
end
for k = 1:length(den_closed)-1
    fprintf('D%d(closed) = %g\n', k, det(H_closed(1:k,1:k)));
end

%% Построение годографа Найквиста
figure('Color', 'w');
nyquist(W_closed);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста замкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_closed); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

%% Построение кривой Михайлова
[num, den] = tfdata(W_closed, 'v');
D = den;
w_max = 5.5;  w_step = 0.005;
w = 0:w_step:w_max;
Djw = polyval(D, 1i*w);

figure('Color','w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)),   imag(Djw(1)),   'go', 'MarkerFaceColor','g', 'MarkerSize',8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor','r', 'MarkerSize',8);
text(real(Djw(1)),   imag(Djw(1)),   '  \omega = 0', 'Color',[0 0.6 0]);
text(real(Djw(end)), imag(Djw(end)), sprintf('  \\omega = %.1f', w_max), 'Color','r');

xline(0, 'k:', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.8);
grid on; box on;
xlabel('Re D(j\omega)');  ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы D(j\omega)');

%% Прямая проверка устойчивости
%% Полюса
fprintf('Полюса РАЗОМКНУТОЙ системы:\n'); disp(pole(W_opened));
fprintf('max Re (разомкнутая) = %g\n', max(real(pole(W_opened))));

fprintf('Полюса ЗАМКНУТОЙ системы:\n');   disp(pole(W_closed));
fprintf('max Re (замкнутая)  = %g\n', max(real(pole(W_closed))));

%% Переходные процессы
t_sim = 0:0.01:10;
figure('Color','w');

subplot(2,1,1);
step(W_opened, t_sim);
grid on;
title('Переходный процесс РАЗОМКНУТОЙ системы');
ylabel('y(t)');

subplot(2,1,2);
step(W_closed, t_sim);
grid on;
title('Переходный процесс ЗАМКНУТОЙ системы');
ylabel('y(t)'); xlabel('Время, с');
