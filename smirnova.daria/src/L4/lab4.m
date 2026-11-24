clc;
clear;
close all;

s = tf('s');
W = 1 / (0.05*s^4 + 0.1*s^3 + s^2 + s + 1);

% Передаточная функция замкнутой системы w_closed= w/(1+w)
W_closed = feedback(W, 1);

figure('Color', 'w');
nyquist(W);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W); 
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

[num, den] = tfdata(W_closed, 'v'); % извлекаем коэффициенты числителя и знаменателя и превращаем в вектор
D = round(den, 10);
w = 0:0.02:20;               
Djw = polyval(D, 1i*w);
 
figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', ...
     'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', ...
     'MarkerSize', 8);
 
text(real(Djw(1)),   imag(Djw(1)),   '   \omega = 0',   'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '   \omega = max', 'Color', 'r');
 
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
 
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega)');
xlim([-5 50]); 
ylim([-20 20]); % Ограничили график, поскольку было не видно поведение кривой из-за масштаба

figure('Color', 'w', 'Name', 'Переходные процессы');
step(W, 'b-', W_closed, 'r-');
grid on;
legend('Разомкнутая система', 'Замкнутая система', 'Location', 'best');
title('Сравнение переходных процессов при ступенчатом воздействии');
xlabel('Время (с)');
ylabel('Амплитуда');

disp('Полюса ЗАМКНУТОЙ системы:');
disp(pole(W_closed));
disp('Полюса РАЗОМКНУТОЙ системы:');
disp(pole(W));
