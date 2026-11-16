clc;
clear;
close all;

s = tf('s');
W_open = 100 / (5*s^4 + 0.1*s^3 + 2*s^2 + 2*s + 1);
W_closed = feedback(W_open, 1);

% Характеристические полиномы
D_open = [5 0.1 2 2 1];
D_closed = [5 0.1 2 2 101];

%% Критерий Гурвица и корни характеристических полиномов
polynomials = {D_open, D_closed};
right_half_plane = zeros(1, 2);
hurwitz_matrices = cell(1, 2);
hurwitz_minors = cell(1, 2);
system_roots = cell(1, 2);
is_stable_hurwitz = false(1, 2);

for system_index = 1:2
    a = polynomials{system_index};
    n = numel(a) - 1;
    H = zeros(n);
    for row = 1:n
        for col = 1:n
            coefficient_index = 2*col - row + 1;
            if coefficient_index >= 1 && coefficient_index <= numel(a)
                H(row, col) = a(coefficient_index);
            end
        end
    end

    minors = zeros(1, n);
    for k = 1:n
        minors(k) = det(H(1:k, 1:k));
    end
    poles = roots(a);
    right_half_plane(system_index) = sum(real(poles) > 1e-8);
    hurwitz_matrices{system_index} = H;
    hurwitz_minors{system_index} = minors;
    system_roots{system_index} = poles;
    is_stable_hurwitz(system_index) = all(minors > 0);
end

%% Годограф Найквиста разомкнутой системы
w_nyquist = [0, logspace(-4, 4, 5000)];
[re_nyquist, im_nyquist] = nyquist(W_open, w_nyquist);
re_nyquist = squeeze(re_nyquist);
im_nyquist = squeeze(im_nyquist);
figure('Color', 'w', 'Name', 'Найквист');
plot(re_nyquist, im_nyquist, 'b-', 'LineWidth', 1.5); hold on;
plot(re_nyquist, -im_nyquist, 'b-', 'LineWidth', 1.5);
grid on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы при коэффициенте усиления 100');

figure('Color', 'w', 'Name', 'Найквист: область критической точки');
plot(re_nyquist, im_nyquist, 'b-', 'LineWidth', 1.5); hold on;
plot(re_nyquist, -im_nyquist, 'b-', 'LineWidth', 1.5);
grid on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста: область точки (-1, 0)');

P = right_half_plane(1);
Z = right_half_plane(2);
% Полный частотный контур: от отрицательных частот через 0 к положительным
w_positive = logspace(-6, 4, 10000);
w_full = [-fliplr(w_positive), 0, w_positive];
W_full = 100 ./ polyval(D_open, 1i*w_full);
nyquist_phase = unwrap(angle(1 + W_full));
N = round((nyquist_phase(end) - nyquist_phase(1)) / (2*pi));
Z_nyquist = P + N;

%% Кривая Михайлова для полинома замкнутой системы
w = linspace(0, 20, 4001);
Djw = polyval(D_closed, 1i*w);
phase_change = unwrap(angle(Djw));
mikhailov_phase_change = phase_change(end) - phase_change(1);

figure('Color', 'w', 'Name', 'Михайлов');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', ...
    'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:');
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:');
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');
legend('D(j\omega)', '\omega = 0', '\omega = 20', 'Location', 'best');

figure('Color', 'w', 'Name', 'Михайлов: начальный участок');
near_origin = w <= 5;
plot(real(Djw(near_origin)), imag(Djw(near_origin)), ...
    'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g');
plot([0 0], [min(imag(Djw(near_origin))) max(imag(Djw(near_origin)))], 'k:');
plot([min(real(Djw(near_origin))) max(real(Djw(near_origin)))], [0 0], 'k:');
grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова: 0 \leq \omega \leq 5');

%% Прямая проверка устойчивости: переходные процессы
t_open_grid = 0:0.01:24;
t_closed_grid = 0:0.01:12;
[y_open, t_open] = step(W_open, t_open_grid);
[y_closed, t_closed] = step(W_closed, t_closed_grid);

figure('Color', 'w', 'Name', 'Переходные процессы');
subplot(2, 1, 1);
plot(t_open, y_open, 'b-', 'LineWidth', 1.3);
grid on;
ylabel('y_{open}(t)');
title('Ступенчатый отклик разомкнутой системы');

subplot(2, 1, 2);
plot(t_closed, y_closed, 'r-', 'LineWidth', 1.3);
grid on;
xlabel('Время, с');
ylabel('y_{closed}(t)');
title('Ступенчатый отклик замкнутой системы');
