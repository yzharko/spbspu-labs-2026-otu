clc;
clear;
close all;

s = tf('s');
W_open = 1 / (0.1*s^3 + 0.1*s^2 + s + 1);
W_closed = feedback(W_open, 1);

[N_open, D_open] = tfdata(W_open, 'v');
[N_closed, D_closed] = tfdata(W_closed, 'v');
D_open = D_open / N_open(end);
D_closed = D_closed / N_closed(end);

disp('Передаточная функция замкнутой системы:');
disp(W_closed);
disp('Коэффициенты D_closed(s) от старшей степени к младшей:');
disp(D_closed);

Delta2_open = D_open(2)*D_open(3) - D_open(1)*D_open(4);
Delta2_closed = D_closed(2)*D_closed(3) - D_closed(1)*D_closed(4);

fprintf('Гурвиц: Delta_2 разомкнутой системы = %.4f\n', Delta2_open);
fprintf('Гурвиц: Delta_2 замкнутой системы   = %.4f\n', Delta2_closed);

w0 = sqrt(10);
radius = 0.05;
w_low = linspace(1e-4, w0 - radius, 1000);
w_high = logspace(log10(w0 + radius), 2, 1000);
[re_low, im_low] = nyquist(W_open, w_low);
[re_high, im_high] = nyquist(W_open, w_high);
z_low = [1, reshape(re_low + 1i*im_low, 1, [])];
z_high = reshape(re_high + 1i*im_high, 1, []);

s_arc = 1i*w0 + radius*exp(1i*linspace(-pi/2, pi/2, 300));
z_arc = 1 ./ polyval(D_open, s_arc);
z_freq = [z_low, NaN, z_high, NaN, conj(z_low), NaN, conj(z_high)];
z_bypass = [z_arc, NaN, conj(z_arc)];

z_contour = [conj(fliplr(z_high)), conj(fliplr(z_arc)), ...
    conj(fliplr(z_low)), z_low, z_arc, z_high, 0, conj(z_high(end))];
phase = unwrap(angle(1 + z_contour));
turns = round((phase(end) - phase(1)) / (2*pi));
fprintf('Найквист: охватов точки (-1, 0) по часовой стрелке = %d\n', -turns);

for k = 1:2
    figure('Color', 'w');
    plot(real(z_freq), imag(z_freq), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(real(z_bypass), imag(z_bypass), 'k--', 'LineWidth', 1.2);
    plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
    grid on;
    xlabel('Re W(j\omega)');
    ylabel('Im W(j\omega)');
    legend('W(j\omega)', 'Обход полюсов', 'Критическая точка', ...
        'Location', 'best');

    if k == 1
        title('Годограф Найквиста разомкнутой системы');
    else
        xlim([-3 2]);
        ylim([-2 2]);
        title('Область критической точки');
    end
end

w = 0:0.02:6;
Djw = polyval(D_closed, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g');
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r');
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0');
text(real(Djw(end)), imag(Djw(end)), '  \omega = 6');
xline(0, 'k:');
yline(0, 'k:');
grid on;
box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');

disp('Полюсы разомкнутой системы:');
disp(pole(W_open));
disp('Полюсы замкнутой системы:');
disp(pole(W_closed));

t = 0:0.01:15;
[y_open, t_open] = step(W_open, t);
[y_closed, t_closed] = step(W_closed, t);

figure('Color', 'w');
subplot(2, 1, 1);
plot(t_open, y_open, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Разомкнутая система: незатухающие колебания');

subplot(2, 1, 2);
plot(t_closed, y_closed, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Замкнутая система: расходящиеся колебания');
