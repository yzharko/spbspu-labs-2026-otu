%% ЛР №1. Временные характеристики динамических звеньев
K = 3.5;
T = 1.2;
T1 = 1.4;
T2 = 1.8;
T1a = 4.0;
syms s t
W = [
    K / s
    K * s
    K
    K / (T * s + 1)
    K / (T2^2 * s^2 + T1a * s + 1)
    K / (T2^2 * s^2 + T1 * s + 1)
    K / (T * s^2 + 1)
    K / (s * (T * s + 1))
    K * s / (T * s + 1)
    K * (T * s + 1)
    K * (T * s + 1) / s
];
names = ["Интегрирующее", "Идеальное дифференцирующее", ...
    "Усилительное", "Апериодическое 1-го порядка", "Апериодическое 2-го порядка", ...
    "Колебательное", "Консервативное", "Инерционное интегрирующее", ...
    "Реальное дифференцирующее", "Форсирующее", "Изодромное"];
time = (0:0.005:30)';

%% Переходные функции и функции веса
stepFig = lab_figure('Переходные характеристики', [850, 950]);
stepTiles = tiledlayout(stepFig, 4, 3, 'TileSpacing', 'compact');
impulseFig = lab_figure('Импульсные характеристики', [850, 950]);
impulseTiles = tiledlayout(impulseFig, 4, 3, 'TileSpacing', 'compact');
analyticOnly = [2 10];
stepDelta = zeros(size(W));
stepDelta([2, 10]) = [K, K * T];
impulseDelta = zeros(size(W));
impulseDelta([3, 9, 10, 11]) = [K, K / T, K, K * T];
errors = nan(numel(W), 2);
for i = 1:numel(W)
    h = ilaplace(W(i) / s, s, t);
    w = ilaplace(W(i), s, t);
    fprintf('№%d. %s: h = %s; w = %s\n', i, names(i), char(h), char(w));

    hValues = regular_values(h, t, time);
    wValues = regular_values(w, t, time);
    hModel = nan(size(time));
    wModel = nan(size(time));
    if ~ismember(i, analyticOnly)
        [num, den] = numden(W(i));
        model = tf(sym2poly(num), sym2poly(den));
        hModel = step(model, time);
        wModel = impulse(model, time);
        errors(i, :) = [max(abs(hValues - hModel)), max(abs(wValues - wModel))];
    end

    nexttile(stepTiles);
    plot_response(time, hValues, hModel, stepDelta(i), 0, i);
    nexttile(impulseTiles);
    plot_response(time, wValues, wModel, impulseDelta(i), stepDelta(i), i);
end
for tiles = [stepTiles, impulseTiles]
    key = legend(nexttile(tiles, 1), 'Аналитика', 'Модель MATLAB', ...
        'Orientation', 'horizontal');
    key.Layout.Tile = 'south';
end
disp(array2table(errors, 'VariableNames', {'Error_h', 'Error_w'}));
checked = setdiff(1:numel(W), analyticOnly);
assert(all(errors(checked, :) < 1e-8, 'all'), 'Расхождение больше допуска 1e-8.');

%% Влияние K на примере звена №6
lab_figure('Влияние K', [850, 420]);
hold on;
for gain = K * [0.5, 1, 2]
    model = tf(gain, [T2^2, T1, 1]);
    plot(time, step(model, time), 'DisplayName', sprintf('K = %g', gain));
end
title('№6. Влияние K');
xlabel('t, с');
ylabel('h(t)');
grid on;
legend;

%% Влияние T на звено №4
lab_figure('Влияние T', [850, 420]);
hold on;
fprintf('\nT, с\tt95 по stepinfo, с\tt95 по формуле, с\n');
for tau = T * [1/5, 1/2, 1, 2, 5]
    model = tf(K, [tau, 1]);
    y = step(model, time);
    info = stepinfo(y, time, K, 'SettlingTimeThreshold', 0.05);
    fprintf('%g\t%.6f\t\t%.6f\n', tau, info.SettlingTime, -tau * log(0.05));
    plot(time, y, 'DisplayName', sprintf('T = %g', tau));
end
yline(0.95 * K, ':', '0.95K', 'HandleVisibility', 'off');
title('№4. Влияние T');
xlabel('t, с');
ylabel('h(t)');
grid on;
legend;

%% Влияние T1 на звено №6
lab_figure('Влияние T1', [850, 420]);
hold on;
T1values = T1 * [1, 1/2, 1/5, 1/10];
zeta = T1values / (2 * T2);
alpha = T1values / (2 * T2^2);
envelope = log(1 ./ (0.05 * sqrt(1 - zeta.^2))) ./ alpha;
longTime = (0:0.005:ceil(max(envelope) + 2 * T2))';
fprintf('\nT1, с\tzeta\t\tПеререг., %%\tУстановление, с\tОгибающая, с\n');
for j = 1:numel(T1values)
    model = tf(K, [T2^2, T1values(j), 1]);
    y = step(model, longTime);
    info = stepinfo(y, longTime, K, 'SettlingTimeThreshold', 0.05);
    fprintf('%g\t%.6f\t%.3f\t\t%.3f\t\t\t%.3f\n', ...
        T1values(j), zeta(j), info.Overshoot, info.SettlingTime, envelope(j));
    plot(longTime, y, 'DisplayName', sprintf('T_1 = %g', T1values(j)));
end
yline(K * [0.95, 1.05], ':', 'HandleVisibility', 'off');
title('№6. Влияние T_1');
xlabel('t, с');
ylabel('h(t)');
grid on;
legend;

%% Предел T1 = 0 и консервативное звено №7
lab_figure('Консервативный предел', [850, 420]);
limitModel = tf(K, [T2^2, 0, 1]);
conservativeModel = tf(K, [T, 0, 1]);
plot(time, step(limitModel, time), time, step(conservativeModel, time), '--');
xlabel('t, с');
ylabel('h(t)');
grid on;
legend('№6 при T_1 = 0', '№7');
fprintf('\nЧастоты: предел №6 = %.6f; №7 = %.6f рад/с\n', 1 / T2, 1 / sqrt(T));

function f = lab_figure(name, dimensions)
    visible = 'on';
    if batchStartupOptionUsed
        visible = 'off';
    end
    f = figure('Name', name, 'NumberTitle', 'off', 'Theme', 'light', ...
        'Position', [60, 60, dimensions], 'Visible', visible);
end

function y = regular_values(expr, t, time)
    expr = subs(expr, [dirac(1, t), dirac(t), heaviside(t)], [0, 0, 1]);
    evaluate = matlabFunction(expr, 'Vars', t);
    values = evaluate(time);
    imagError = max(abs(imag(values)), [], 'all');
    assert(imagError < 1e-8, 'Мнимая часть отклика достигает %.3g.', imagError);
    y = real(values) + zeros(size(time));
end

function plot_response(time, analytic, model, delta, deltaPrime, link)
    plot(time, analytic, 'b-', time, model, 'r--', 'LineWidth', 1.2);
    if max(analytic) - min(analytic) < 1e-8
        ylim(analytic(1) + [-0.1, 0.1]);
    end
    if delta ~= 0 || deltaPrime ~= 0
        xline(0, ':r', 'HandleVisibility', 'off');
        terms = [string(delta) + " \delta", string(deltaPrime) + " \delta'"];
        label = join(terms([delta, deltaPrime] ~= 0), " + ");
        text(0.05, 0.17, label, 'Units', 'normalized', 'Color', [0.7, 0, 0]);
    end
    title(sprintf('№%d', link));
    xlabel('t, с');
    grid on;
    xlim([-0.5, time(end)]);
end
