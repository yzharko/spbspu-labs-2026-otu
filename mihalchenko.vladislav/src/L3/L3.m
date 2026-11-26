clear;
clc;
close all;

K = 1.1;
T = 0.2;
T1 = 0.3;
T2 = 0.5;
T1a = 1.5;
s = tf('s');
w = logspace(-2, 2, 1601);
out = fullfile(fileparts(mfilename('fullpath')), 'results_L3');
if ~isfolder(out)
    mkdir(out);
end

names = {
    'Интегрирующее'
    'Дифференцирующее'
    'Усилительное'
    'Апериодическое 1-го порядка'
    'Апериодическое 2-го порядка'
    'Колебательное'
    'Консервативное'
    'Инерционное интегрирующее'
    'Реальное дифференцирующее'
    'Форсирующее 1-го порядка'
    'Изодромное'
};
W = {
    K / s
    K * s
    tf(K)
    K / (T * s + 1)
    K / (T2^2 * s^2 + T1a * s + 1)
    K / (T2^2 * s^2 + T1 * s + 1)
    K / (T * s^2 + 1)
    K / (s * (T * s + 1))
    K * s / (T * s + 1)
    K * (T * s + 1)
    K * (T * s + 1) / s
};

%% АФЧХ, ЛАЧХ и ЛФЧХ всех звеньев
for n = 1:numel(W)
    f = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
        'Position', [100 100 1050 290]);
    tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    resonance = 0;
    if n == 7
        resonance = 1 / sqrt(T);
    end
    wf = frequencyGrid(w, resonance);
    valid = isfinite(wf);
    [re, im] = nyquist(W{n}, wf(valid));
    u = nan(size(wf));
    v = u;
    u(valid) = re(:);
    v(valid) = im(:);

    nexttile;
    if n == 3
        plot(K, 0, 'o', 'MarkerFaceColor', [0 0.447 0.741]);
        axis([0 2*K -K K]);
    else
        plot(u, v, 'LineWidth', 1.35);
        if ismember(n, [4 5 6 9])
            axis equal;
        elseif n == 7
            axis([-5*K 5*K -1 1]);
        else
            xlim([min(u) - K/4, max(u) + K/4]);
            ylim([-4*K 4*K]);
        end
    end
    grid on;
    title('АФЧХ');
    xlabel('Re W(j\omega)');
    ylabel('Im W(j\omega)');
    drawBode(W(n), w, resonance, {});
    sgtitle(sprintf('%d. %s', n, names{n}));
    exportgraphics(f, fullfile(out, sprintf('base_%02d.png', n)), 'Resolution', 185);
    close(f);
end

%% Влияние K при фиксированных постоянных времени
Kvalues = K * [0.5 1 2];
for n = 1:numel(W)
    models = cell(size(Kvalues));
    for j = 1:numel(Kvalues)
        models{j} = (Kvalues(j) / K) * W{n};
    end
    resonances = zeros(size(Kvalues));
    if n == 7
        resonances(:) = 1 / sqrt(T);
    end
    f = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
        'Position', [100 100 1050 290]);
    tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    drawBode(models, w, resonances, compose('K = %g', Kvalues));
    sgtitle(sprintf('%d. %s: влияние K', n, names{n}));
    exportgraphics(f, fullfile(out, sprintf('gain_%02d.png', n)), 'Resolution', 185);
    close(f);
end

%% Влияние постоянных времени при фиксированном K
Tvalues = T * [0.5 1 2];
studies = {
    4, 'T', Tvalues, @(p) K / (p * s + 1);
    5, 'T_{1a}', [1.1 T1a 2], @(p) K / (T2^2 * s^2 + p * s + 1);
    5, 'T_2', [0.25 T2 0.7], @(p) K / (p^2 * s^2 + T1a * s + 1);
    6, 'T_1', T1 * [0.5 1 2], @(p) K / (T2^2 * s^2 + p * s + 1);
    6, 'T_2', T2 * [0.5 1 2], @(p) K / (p^2 * s^2 + T1 * s + 1);
    7, 'T', Tvalues, @(p) K / (p * s^2 + 1);
    8, 'T', Tvalues, @(p) K / (s * (p * s + 1));
    9, 'T', Tvalues, @(p) K * s / (p * s + 1);
    10, 'T', Tvalues, @(p) K * (p * s + 1);
    11, 'T', Tvalues, @(p) K * (p * s + 1) / s
};
for k = 1:size(studies, 1)
    [n, parameter, values, makeModel] = studies{k, :};
    models = cell(size(values));
    for j = 1:numel(values)
        models{j} = makeModel(values(j));
    end
    resonances = zeros(size(values));
    if n == 7
        resonances = 1 ./ sqrt(values);
    end
    f = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
        'Position', [100 100 1050 290]);
    tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    drawBode(models, w, resonances, compose([parameter ' = %g'], values));
    sgtitle(sprintf('%d. %s: влияние %s', n, names{n}, parameter));
    exportgraphics(f, fullfile(out, sprintf('time_%02d.png', k)), 'Resolution', 185);
    close(f);
end

fprintf('Графики сохранены в %s\n', out);

function drawBode(models, w, resonances, labels)
    axL = nexttile;
    axP = nexttile;
    hold(axL, 'on');
    hold(axP, 'on');
    styles = {'--', '-', '-.'};
    if isscalar(models)
        styles = {'-'};
    end
    for k = 1:numel(models)
        wf = frequencyGrid(w, resonances(k));
        valid = isfinite(wf);
        [mag, phase] = bode(models{k}, wf(valid));
        L = nan(size(wf));
        phi = L;
        L(valid) = 20 * log10(mag(:));
        phi(valid) = phase(:);
        if resonances(k) > 0
            phi(valid) = mod(phi(valid) + 180, 360) - 180;
        end
        semilogx(axL, wf, L, styles{k}, 'LineWidth', 1.35);
        semilogx(axP, wf, phi, styles{k}, 'LineWidth', 1.35);
    end
    set([axL axP], 'XScale', 'log', 'XLim', [w(1) w(end)]);
    yline(axL, 0, ':', 'HandleVisibility', 'off');
    grid(axL, 'on');
    grid(axP, 'on');
    title(axL, 'ЛАЧХ');
    title(axP, 'ЛФЧХ');
    xlabel(axL, '\omega, рад/с');
    xlabel(axP, '\omega, рад/с');
    ylabel(axL, 'L, дБ');
    ylabel(axP, '\phi, град');
    if ~isempty(labels)
        legend(axL, labels, 'Location', 'best');
        legend(axP, labels, 'Location', 'best');
    end
end

function wf = frequencyGrid(w, resonance)
    wf = w;
    if resonance > 0
        edges = resonance * [0.999 1.001];
        wf = [w(w < edges(1)), edges(1), NaN, edges(2), w(w > edges(2))];
    end
end
