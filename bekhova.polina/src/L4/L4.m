clc;
clear;
close all;

s = tf('s');
W_open = 100 / (5*s^4 + 0.1*s^3 + 2*s^2 + 2*s + 1);
W_closed = 100 / (5*s^4 + 0.1*s^3 + 2*s^2 + 2*s + 101);

D_closed = [5, 0.1, 2, 2, 101];
disp('Характеристическое уравнение замкнутой системы:');
disp('D(s) = 5*s^4 + 0.1*s^3 + 2*s^2 + 2*s + 101 = 0');

a0=5; a1=0.1; a2=2; a3=2; a4=1;
H_open = [a1, a3, 0, 0;
    a0, a2, a4, 0;
    0,  a1, a3, 0;
    0,  a0, a2, a4];
D1_open = det(H_open(1:1,1:1));
D2_open = det(H_open(1:2,1:2));
disp('Критерий Гурвица (разомкнутая):');
disp(['D1 = ', num2str(D1_open), ', D2 = ', num2str(D2_open)]);
if D1_open>0 && D2_open>0, disp('Устойчива'); else disp('Неустойчива'); end

a4=101;
H_closed = [a1, a3, 0, 0;
    a0, a2, a4, 0;
    0,  a1, a3, 0;
    0,  a0, a2, a4];
D1_closed = det(H_closed(1:1,1:1));
D2_closed = det(H_closed(1:2,1:2));
disp('Критерий Гурвица (замкнутая):');
disp(['D1 = ', num2str(D1_closed), ', D2 = ', num2str(D2_closed)]);
if D1_closed>0 && D2_closed>0, disp('Устойчива'); else disp('Неустойчива'); end

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_open); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);
w = 0:0.001:10;              
Djw = polyval(D, 1i*w);

scale = max(abs(Djw));
Re_norm = real(Djw) / scale;
Im_norm = imag(Djw) / scale;

figure('Color', 'w');
plot(Re_norm, Im_norm, 'b-', 'LineWidth', 1.8); hold on;
plot(Re_norm(1), Im_norm(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(Re_norm(end), Im_norm(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(Re_norm(1), Im_norm(1), '  \omega = 0', 'Color', 'g');
text(Re_norm(end), Im_norm(end), '  \omega = max', 'Color', 'r');
plot([0 0], [min(Im_norm) max(Im_norm)], 'k:', 'LineWidth', 0.8);
plot([min(Re_norm) max(Re_norm)], [0 0], 'k:', 'LineWidth', 0.8);
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega)');
axis equal;

disp('Корни характеристического уравнения:');
disp(roots(D_closed));

t = 0:0.01:5;
figure('Color', 'w');
step(W_open, t); hold on;
step(W_closed, t);
grid on;
legend('Разомкнутая система', 'Замкнутая система');
title('Переходные процессы при ступенчатом воздействии');
xlabel('Время, с');
ylabel('Выход, y(t)');
