clear; clc; close all;

syms s t

Wopen = 1 / (8*s^3 + 4*s^2 + 2*s + 1);
Wclose = simplifyFraction(Wopen/(1 + Wopen), "Expand", true);
disp(Wclose); % 1/(8*s^3 + 4*s^2 + 2*s + 2)

[Nopen, Dopen] = numden(Wopen);
[Nclose, Dclose] = numden(Wclose);

polyNopen = sym2poly(Nopen);
polyDopen = sym2poly(Dopen);
polyDclose = sym2poly(Dclose);

% Критерий Гурвица
disp(gurv(Dopen));  % 4   0   0 - граница колебательной устойчивости
disp(gurv(Dclose)); % 4  -8 -16 - неустойчива :(

% Кривая Найквиста
figure("Name", "Найквист");
nyquist(tf(polyNopen, polyDopen));

% Кривая Михайлова.
mikhailov(Dclose, 0:0.01:2, [-2 2], [-2 2]);

% Прямая проверка
t_range = 0:0.1:100;
stepcheck("Открытая", Wopen, t, s, t_range); 
stepcheck("Замкнутая", Wclose, t, s, t_range);

% Проверка по корням
disp(roots(polyDopen));  %    -0.5;          +- 0.5i - автоколебания
disp(roots(polyDclose)); % -0.6766; 0.883 +- 0.6014i - расходящийся

function g = gurv(d)
    a = sym2poly(d);
    gurvMatrix = [
        a(2), a(4),    0
        a(1), a(3),    0
           0, a(2), a(4)
    ];
    first = a(2);
    second = det(gurvMatrix(1:2, 1:2));
    third = second * a(4);
    g = [first, second, third];
end

function mikhailov(D, w_range, xLim, yLim)
    Df = matlabFunction(D);
    
    figure("Name", "Михайлов");
    Djw = Df(w_range*1i);
    plot(real(Djw), imag(Djw)); hold on;
    plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', ...
         'MarkerSize', 8);
    clampfilter = @(x, l) x >= l(1) & x <= l(2);
    filter = clampfilter(real(Djw), xLim) & clampfilter(imag(Djw), yLim);
    filtered = Djw(filter);
    endx = real(filtered(end));
    endy = imag(filtered(end));
    plot(endx, endy, 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
     
    text(real(Djw(1)), imag(Djw(1)),   '  \omega = 0',   'Color', 'g');
    text(endx, endy, '  \omega = max', 'Color', 'r');
     
    plot(xLim, [0 0], 'k', 'LineWidth', 0.8);
    plot([0 0], yLim, 'k', 'LineWidth', 0.8);
     
    grid on; box on;
    xlabel('Re D(j\omega)');
    ylabel('Im D(j\omega)');
    xlim(xLim);
    ylim(yLim);
    title('Кривая Михайлова D(j\omega)');
    axis square;
    hold off;
end

function stepcheck(name, w, t, s, t_range)
    WStep = laplace(heaviside(t), t, s);
    WRes = w * WStep;
    f = ilaplace(WRes, s, t);
    figure("Name", name);
    plot(t_range, double(subs(f, t, t_range)));
    title(name);
    xlabel("Время");
    ylabel("Выход")
end

