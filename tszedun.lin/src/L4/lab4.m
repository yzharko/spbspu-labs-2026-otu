function lab4

% Вариант 5
Dopen = [8 4 2 1];
Dclosed = Dopen + [0 0 0 1];   % Единичная отрицательная обратная связь.
W = {tf(1, Dopen), tf(1, Dclosed)};

% Годограф Найквиста разомкнутой системы
wp = 0.5;
radius = 0.025;
w = [linspace(0, wp-radius, 250), logspace(log10(wp+radius), log10(30), 600)];
[re, im] = nyquist(W{1}, w);
z = reshape(complex(re, im), 1, []);
low = z(w < wp);
high = z(w > wp);

% Дуга справа от полюса +0,5j. Для отрицательных частот берём сопряжение.
theta = linspace(-pi/2, pi/2, 181);
arc = 1 ./ polyval(Dopen, 1i*wp + radius*exp(1i*theta));
branch = [low NaN high 0];
positive = [low arc high 0];
contour = [conj(fliplr(positive)) positive];
phase = unwrap(angle(1 + contour));
turns = (phase(1) - phase(end)) / (2*pi);

newFigure('Найквист');
plot(branch, 'b-');
hold on;
plot(conj(branch), 'Color', [0.1 0.53 0.33]);
plot([arc NaN conj(arc)], 'k--');
plot(-1, 0, 'r+', 'MarkerSize', 12);
mid = round(numel(arc)/2);
tail = arc(mid-6);
tip = arc(mid+6);
base = [tail conj(tip)];
direction = [tip conj(tail)] - base;
quiver(real(base), imag(base), real(direction), imag(direction), 0, 'k', ...
       'MaxHeadSize', 0.18, 'HandleVisibility', 'off');
axis equal; axis([-8 8 -8 8]);
xlabel('Re W, отн. ед.'); ylabel('Im W, отн. ед.');
title('Годограф Найквиста разомкнутой системы');
legend('\omega \geq 0', '\omega \leq 0', 'Дуги обхода', 'Точка (−1, 0)', ...
       'Location', 'northeastoutside');

% Кривая Михайлова замкнутой системы
curve = polyval(Dclosed, 1i*(0:0.002:0.9));
points = polyval(Dclosed, 1i*[0 0.5 1/sqrt(2)]);
newFigure('Михайлов');
plot(curve, 'b-');
hold on;
plot(points, 'ko', 'MarkerFaceColor', 'w');
text(real(points), imag(points), {'  \omega = 0', '  \omega = 0,5', '  \omega = 1/√2'}, ...
     'VerticalAlignment', 'top');
axis([-1.5 2.7 -4.3 0.8]);
line(xlim, [0 0], 'Color', 'k', 'LineStyle', ':');
line([0 0], ylim, 'Color', 'k', 'LineStyle', ':');
xlabel('Re D(j\omega), отн. ед.'); ylabel('Im D(j\omega), отн. ед.');
title('Кривая Михайлова замкнутой системы');

% Критерий Гурвица, корни и переходные процессы
den = [Dopen; Dclosed];
names = {'Разомкнутая система', 'Замкнутая система'};
styles = {'b-', 'r-'};
t = 0:0.05:50;
newFigure('Переходные процессы');
for k = 1:2
    d = den(k,:);
    H = [d(2) d(4) 0; d(1) d(3) 0; 0 d(2) d(4)];
    fprintf('\n%s\nМатрица Гурвица:\n', names{k});
    disp(H);
    delta2 = d(2)*d(3) - d(1)*d(4);
    fprintf('Главные миноры: %g, %g, %g\n', d(2), delta2, d(4)*delta2);
    fprintf('Корни характеристического уравнения:\n');
    disp(roots(d));

    subplot(2, 1, k);
    plot(t, step(W{k}, t), styles{k});
    ylabel('y(t), отн. ед.');
    title(names{k});
end
xlabel('t, с');
fprintf('Число обходов точки (−1, 0) по часовой стрелке: %.0f\n', turns);
end

function newFigure(name)
f = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 850 610]);
if isprop(f, 'Theme')
    set(f, 'Theme', 'light');
end
set(f, 'DefaultAxesXGrid', 'on', 'DefaultAxesYGrid', 'on', 'DefaultLineLineWidth', 1.6);
end
