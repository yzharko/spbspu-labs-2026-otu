clear; clc; close all;

s = tf('s');
W_open = 10/(2*s^5 + 3*s^4 + 3*s^3 + 0.5*s^2 + 0.5*s + 1);
W_closed = feedback(W_open, 1);

[~, A] = tfdata(W_open, 'v');
[~, D] = tfdata(W_closed, 'v');
fprintf('A(s): %s\n', num2str(A));
fprintf('D(s): %s\n', num2str(D));

system_names = {'разомкнутой', 'замкнутой'};
systems = {W_open, W_closed};
char_polys = {A, D};
for k = 1:2
    a = char_polys{k};
    hurwitz = [a(2) a(4) a(6) 0 0;
        a(1) a(3) a(5) 0 0;
        0 a(2) a(4) a(6) 0;
        0 a(1) a(3) a(5) 0;
        0 0 a(2) a(4) a(6)];
    minors = zeros(1, 5);
    for m = 1:5
        minors(m) = det(hurwitz(1:m, 1:m));
    end
    fprintf('\nМатрица Гурвица %s системы:\n', system_names{k});
    disp(hurwitz);
    fprintf('Главные миноры: %s\n', num2str(minors));
    fprintf('Полюсы %s системы:\n', system_names{k});
    disp(pole(systems{k}));
end

nyquist_zoom = [-2 2];
nyquist_titles = {'Весь годограф', 'Окрестность точки (-1; j0)'};
figure('Name', 'Найквист', 'Position', [50 50 800 420]);
layout_nyquist = tiledlayout(1, 2, 'TileSpacing', 'compact');
title(layout_nyquist, 'Годограф Найквиста разомкнутой системы');
for k = 1:2
    nexttile;
    nyquist(W_open);
    hold on;
    plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    title(nyquist_titles{k});
    xlabel('Re W(j\omega)');
    ylabel('Im W(j\omega)');
end
xlim(nyquist_zoom);
ylim(nyquist_zoom);

omega_end = [1.5 0.6];
figure('Name', 'Михайлов', 'Position', [50 50 800 420]);
layout_mikhailov = tiledlayout(1, 2, 'TileSpacing', 'compact');
title(layout_mikhailov, 'Кривая Михайлова замкнутой системы');
for k = 1:2
    omega = linspace(0, omega_end(k), 2000);
    D_jw = polyval(D, 1i*omega);
    nexttile;
    plot(real(D_jw), imag(D_jw), 'LineWidth', 1.5);
    hold on;
    plot(real(D_jw(1)), imag(D_jw(1)), 'o', 'LineWidth', 1.5);
    yline(0, 'k');
    axis padded;
    if k == 1
        xline(0, 'k');
        xlim([0 max(xlim)]);
    end
    grid on;
    title(sprintf('\\omega от 0 до %g рад/с', omega_end(k)));
    xlabel('X(\omega)');
    ylabel('Y(\omega)');
end

time_end = 25;
t = linspace(0, time_end, 2500);
step_titles = {'Разомкнутая система', 'Замкнутая система'};
figure('Name', 'Переходные процессы', 'Position', [50 50 800 420]);
layout_step = tiledlayout(1, 2, 'TileSpacing', 'compact');
title(layout_step, 'Переходные процессы при ступенчатом воздействии');
for k = 1:2
    h = step(systems{k}, t);
    nexttile;
    plot(t, h, 'LineWidth', 1.5);
    grid on;
    title(step_titles{k});
    xlabel('t, с');
    ylabel('h(t)');
end
