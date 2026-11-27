clear;
clc;
close all;

K0 = 3.2;
T0 = 0.9;
T1Oscillatory = 1.0;
T2 = 1.4;
T1Aperiodic = 3.0;

names = {
    'Интегрирующее'
    'Идеальное дифференцирующее'
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

omegaBase = logspace(-3, 3, 1601).';
s = tf('s');
models = makeLinks(s, K0, T0, T1Aperiodic, T1Oscillatory, T2);
resultsFolder = fullfile(fileparts(mfilename('fullpath')), 'results_L3');

if ~isfolder(resultsFolder)
    mkdir(resultsFolder);
end

[figNyquist, nyquistTiles] = newFigure('АФЧХ: все 11 звеньев');
[figMagnitude, magnitudeTiles] = newFigure('ЛАЧХ: все 11 звеньев');
[figPhase, phaseTiles] = newFigure('ЛФЧХ: все 11 звеньев');

fprintf('Проверка bode и nyquist прямой подстановкой s = j*omega\n');
fprintf('Звено   Ошибка bode    Ошибка nyquist\n');

for id = 1:numel(models)
    model = models{id};
    isConservative = id == 7;
    [omega, poleOmega] = frequencyGrid(model, omegaBase, isConservative);
    [level, phase, bodeError] = bodeData(model, omega, poleOmega);
    [U, V, nyquistError] = nyquistData(model, omega, poleOmega);

    ax = nexttile(nyquistTiles);
    plot(ax, U, V, '-', 'LineWidth', 1.4);
    hold(ax, 'on');
    plot(ax, U, -V, '--', 'Color', [1 0.65 0.20], 'LineWidth', 1.6);
    grid(ax, 'on');
    xlabel(ax, 'Re W(j\omega)');
    ylabel(ax, 'Im W(j\omega)');
    title(ax, linkTitle(id, names));

    if id == 3
        plot(ax, U(1), V(1), 'o', 'MarkerFaceColor', [0 0.45 0.74]);
    elseif ismember(id, [1 2 7 8 10 11])
        xlim(ax, [-4 5]);
        ylim(ax, [-12 12]);
    else
        axis(ax, 'equal');
    end

    ax = nexttile(magnitudeTiles);
    semilogx(ax, omega, level, 'LineWidth', 1.4);
    labelBode(ax, omegaBase, 'L, дБ', linkTitle(id, names));
    yline(ax, 0, ':');
    markPole(ax, poleOmega);

    ax = nexttile(phaseTiles);
    semilogx(ax, omega, phase, 'LineWidth', 1.4);
    labelBode(ax, omegaBase, '\psi, град', linkTitle(id, names));
    markPole(ax, poleOmega);

    fprintf('%5d   %11.3e   %13.3e\n', id, bodeError, nyquistError);
end

title(nyquistTiles, ...
    'АФЧХ: сплошная линия — \omega>0, штриховая — \omega<0', 'Color', 'w');
formatFigure(figNyquist);
formatFigure(figMagnitude);
formatFigure(figPhase);

gainValues = [0.5 1 K0 5];
[figGain, gainTiles] = newFigure('Влияние K на ЛАЧХ всех звеньев');
gainColors = lines(numel(gainValues));
gainData = cell(numel(models) * numel(gainValues), 4);
row = 0;

for id = 1:numel(models)
    model = models{id};
    ax = nexttile(gainTiles);
    hold(ax, 'on');
    [omega, poleOmega] = frequencyGrid(model, omegaBase, id == 7);
    [levelBase, phaseBase] = bodeData(model, omega, poleOmega);
    active = isfinite(levelBase);

    for j = 1:numel(gainValues)
        K = gainValues(j);
        gainModel = (K / K0) * model;
        [level, phase] = bodeData(gainModel, omega, poleOmega);

        assert(max(abs(level(active) - levelBase(active) - 20 * log10(K / K0))) < 1e-6, ...
            'Неверный сдвиг ЛАЧХ при изменении K для звена %d.', id);
        phaseDiff = angle(exp(1i * deg2rad(phase(active) - phaseBase(active))));
        assert(max(abs(phaseDiff)) < 1e-8, ...
            'Фаза изменилась при изменении K для звена %d.', id);

        semilogx(ax, omega, level, 'Color', gainColors(j, :), 'LineWidth', 1.3, ...
            'DisplayName', sprintf('K = %g', K));

        [num, den] = tfdata(gainModel, 'v');
        levelAtOne = 20 * log10(abs(polyval(num, 1i) / polyval(den, 1i)));
        row = row + 1;
        gainData(row, :) = {id, K, levelAtOne, zeroDbFrequencies(num, den)};
    end

    labelBode(ax, omegaBase, 'L, дБ', linkTitle(id, names));
    yline(ax, 0, ':', 'HandleVisibility', 'off');
    xline(ax, 1, ':', 'HandleVisibility', 'off');
    markPole(ax, poleOmega);

    if id == 1
        legend(ax, 'Location', 'best', 'FontSize', 8);
    end
end

formatFigure(figGain);

gainTable = cell2table(gainData, 'VariableNames', ...
    {'Link', 'K', 'L_at_1_rad_s_dB', 'Omega_0dB_rad_s'});
writetable(gainTable, fullfile(resultsFolder, 'gain_intersections.csv'));

timeStudies = {
     4, 'T',  T0 * [0.5 1 2]
     5, 'T1', [2*T2 T1Aperiodic 2*T1Aperiodic]
     5, 'T2', [T2/2 T2 T1Aperiodic/2]
     6, 'T1', T1Oscillatory * [0.5 1 2]
     6, 'T2', T2 * [0.5 1 2]
     7, 'T',  T0 * [0.5 1 2]
     8, 'T',  T0 * [0.5 1 2]
     9, 'T',  T0 * [0.5 1 2]
    10, 'T',  T0 * [0.5 1 2]
    11, 'T',  T0 * [0.5 1 2]
};

[figTimeL, timeLTiles] = newFigure('Влияние постоянных времени на ЛАЧХ');
[figTimeP, timePTiles] = newFigure('Влияние постоянных времени на ЛФЧХ');
timeColors = lines(3);
timeData = cell(size(timeStudies, 1) * 3, 7);
row = 0;

for i = 1:size(timeStudies, 1)
    id = timeStudies{i, 1};
    parameter = timeStudies{i, 2};
    values = timeStudies{i, 3};
    axL = nexttile(timeLTiles);
    axP = nexttile(timePTiles);
    hold(axL, 'on');
    hold(axP, 'on');

    for j = 1:numel(values)
        value = values(j);
        timeModel = varyTime(s, K0, T0, T1Aperiodic, T1Oscillatory, T2, id, parameter, value);
        [omega, poleOmega] = frequencyGrid(timeModel, omegaBase, id == 7);
        [level, phase] = bodeData(timeModel, omega, poleOmega);
        [asymL, asymP, corners] = asymptotes(timeModel, omega, poleOmega);

        label = sprintf('%s = %.3g', parameter, value);
        semilogx(axL, omega, level, 'Color', timeColors(j, :), 'LineWidth', 1.4, ...
            'DisplayName', label);
        semilogx(axP, omega, phase, 'Color', timeColors(j, :), 'LineWidth', 1.4, ...
            'DisplayName', label);
        semilogx(axL, omega, asymL, '--', 'Color', timeColors(j, :), ...
            'HandleVisibility', 'off');
        semilogx(axP, omega, asymP, '--', 'Color', timeColors(j, :), ...
            'HandleVisibility', 'off');

        row = row + 1;

        if isempty(poleOmega)
            phaseStart = corners / 10;
            phaseEnd = 10 * corners;
        else
            phaseStart = poleOmega;
            phaseEnd = poleOmega;
        end

        timeData(row, :) = {id, parameter, value, numberList(corners), ...
            numberList(diff(log10(corners))), numberList(phaseStart), numberList(phaseEnd)};
    end

    caption = sprintf('%d. Изменяем %s', id, parameter);
    labelBode(axL, omegaBase, 'L, дБ', caption);
    labelBode(axP, omegaBase, '\psi, град', caption);
    legend(axL, 'Location', 'best', 'FontSize', 8);
    legend(axP, 'Location', 'best', 'FontSize', 8);
end

title(timeLTiles, ...
    'Влияние постоянных времени на ЛАЧХ: штрихом показаны асимптоты', 'Color', 'w');
title(timePTiles, ...
    'Влияние постоянных времени на ЛФЧХ: штрихом показаны приближения', 'Color', 'w');
formatFigure(figTimeL);
formatFigure(figTimeP);

timeTable = cell2table(timeData, 'VariableNames', ...
    {'Link', 'Parameter', 'Value', 'Corner_rad_s', 'Gap_between_corners_decades', ...
    'Phase_transition_start_rad_s', 'Phase_transition_end_rad_s'});
writetable(timeTable, fullfile(resultsFolder, 'time_corners.csv'));

fprintf('\nЗвено № 4 при K = %.1f и T = %.1f:\n', K0, T0);
fprintf('Сопрягающая частота: %.4f рад/с\n', 1 / T0);
fprintf('Низкочастотный уровень: %.4f дБ\n', 20 * log10(K0));
fprintf('Пересечение уровня 0 дБ: %.4f рад/с\n', sqrt(K0^2 - 1) / T0);

drawnow;
fprintf('\nТаблицы сохранены в %s\n', resultsFolder);

function models = makeLinks(s, K, T, T1Aperiodic, T1Oscillatory, T2)
    models = {
        K / s
        K * s
        tf(K)
        K / (T * s + 1)
        K / (T2^2 * s^2 + T1Aperiodic * s + 1)
        K / (T2^2 * s^2 + T1Oscillatory * s + 1)
        K / (T * s^2 + 1)
        K / (s * (T * s + 1))
        K * s / (T * s + 1)
        K * (T * s + 1)
        K * (T * s + 1) / s
    };
end

function model = varyTime(s, K, T, T1Aperiodic, T1Oscillatory, T2, id, parameter, value)
    switch parameter
        case 'T'
            T = value;
        case 'T1'
            if id == 5
                T1Aperiodic = value;
            else
                T1Oscillatory = value;
            end
        case 'T2'
            T2 = value;
    end

    models = makeLinks(s, K, T, T1Aperiodic, T1Oscillatory, T2);
    model = models{id};
end

function [omega, poleOmega] = frequencyGrid(model, omegaBase, isConservative)
    omega = omegaBase;
    poleOmega = [];

    if isConservative
        [~, den] = tfdata(model, 'v');
        poleOmega = sqrt(den(end) / den(1));
        offsets = logspace(-4, -1, 65).';
        omega = unique([omegaBase; poleOmega * (1 - offsets); poleOmega; poleOmega * (1 + offsets)]);
    end
end

function [level, phase, err] = bodeData(model, omega, poleOmega)
    validFrequencies = frequencyMask(omega, poleOmega);
    [magnitude, phaseDegrees] = bode(model, omega(validFrequencies));
    magnitude = reshape(magnitude, [], 1);
    phaseDegrees = reshape(phaseDegrees, [], 1);
    err = checkResponse(model, omega(validFrequencies), magnitude .* exp(1i * deg2rad(phaseDegrees)), 'bode');

    level = nan(size(omega));
    phase = nan(size(omega));
    level(validFrequencies) = 20 * log10(magnitude);
    phase(validFrequencies) = phaseDegrees;

    if ~isempty(poleOmega)
        phase(validFrequencies & omega > poleOmega) = -180;
    end
end

function [U, V, err] = nyquistData(model, omega, poleOmega)
    validFrequencies = frequencyMask(omega, poleOmega);
    [realPart, imagPart] = nyquist(model, omega(validFrequencies));
    realPart = reshape(realPart, [], 1);
    imagPart = reshape(imagPart, [], 1);
    err = checkResponse(model, omega(validFrequencies), complex(realPart, imagPart), 'nyquist');
    U = nan(size(omega));
    V = nan(size(omega));
    U(validFrequencies) = realPart;
    V(validFrequencies) = imagPart;
end

function validFrequencies = frequencyMask(omega, poleOmega)
    validFrequencies = true(size(omega));

    if ~isempty(poleOmega)
        validFrequencies = abs(omega / poleOmega - 1) > 1e-12;
    end
end

function err = checkResponse(model, omega, actual, command)
    [num, den] = tfdata(model, 'v');
    expected = polyval(num, 1i * omega) ./ polyval(den, 1i * omega);
    err = max(abs(actual - expected) ./ max(1, abs(expected)));
    assert(isfinite(err) && err < 1e-7, ...
        'Результат %s расходится с прямой подстановкой: ошибка %.3g.', command, err);
end

function result = zeroDbFrequencies(num, den)
    p = magnitudeSquaredPolynomial(num);
    q = magnitudeSquaredPolynomial(den);
    n = max(numel(p), numel(q));
    eq = [zeros(1, n - numel(p)), p] - [zeros(1, n - numel(q)), q];
    first = find(abs(eq) > 1e-12 * max(abs([p q])), 1);

    if isempty(first)
        result = 'все положительные';
        return
    end

    candidates = roots(eq(first:end));
    positive = real(candidates) > 0 & ...
        abs(imag(candidates)) < 1e-9 * max(1, abs(real(candidates)));
    crossing = sort(sqrt(real(candidates(positive))));
    result = numberList(crossing);
end

function p = magnitudeSquaredPolynomial(coefficients)
    powers = (numel(coefficients) - 1):-1:0;
    atJ = coefficients .* (1i).^powers;
    squared = real(conv(atJ, conj(atJ)));
    p = squared(1:2:end);
end

function [level, phase, corners] = asymptotes(model, omega, poleOmega)
    [num, den] = tfdata(model, 'v');
    z = roots(num);
    p = roots(den);
    orderAtZero = sum(z == 0) - sum(p == 0);
    gain = abs(num(find(num ~= 0, 1, 'last')) / den(find(den ~= 0, 1, 'last')));
    level = 20 * log10(gain) + 20 * orderAtZero * log10(omega);
    phase = 90 * orderAtZero * ones(size(omega));
    z = z(z ~= 0);
    p = p(p ~= 0);
    cornerValues = abs([z; p]);
    corners = sort(uniquetol(cornerValues, 1e-7));

    for wc = abs(z).'
        [deltaL, deltaP] = factorAsymptotes(omega, wc);
        level = level + deltaL;
        phase = phase + deltaP;
    end

    for wc = abs(p).'
        [deltaL, deltaP] = factorAsymptotes(omega, wc);
        level = level - deltaL;
        phase = phase - deltaP;
    end

    if ~isempty(poleOmega)
        phase = -180 * double(omega > poleOmega);
        phase(~frequencyMask(omega, poleOmega)) = NaN;
    end
end

function [level, phase] = factorAsymptotes(omega, wc)
    x = log10(omega / wc);
    level = 20 * max(x, 0);
    phase = 90 * max(0, min(1, (x + 1) / 2));
end

function labelBode(ax, omegaBase, verticalLabel, caption)
    set(ax, 'XScale', 'log');
    xlim(ax, [omegaBase(1) omegaBase(end)]);
    grid(ax, 'on');
    xlabel(ax, '\omega, рад/с');
    ylabel(ax, verticalLabel);
    title(ax, caption);
end

function markPole(ax, poleOmega)
    if ~isempty(poleOmega)
        xline(ax, poleOmega, ':', '\omega_0', 'Color', 'w', 'Alpha', 1, ...
            'HandleVisibility', 'off');
    end
end

function [fig, tiles] = newFigure(caption)
    fig = figure('Name', caption, 'Color', [0.10 0.10 0.10], 'Position', [35 35 1460 980]);
    tiles = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tiles, caption, 'Color', 'w');
end

function caption = linkTitle(id, names)
    caption = sprintf('%d. %s', id, names{id});
end

function formatFigure(fig)
    axesHandles = findall(fig, 'Type', 'axes');
    axesBackground = [0.07 0.07 0.07];
    gridColor = [0.70 0.70 0.70];

    for i = 1:numel(axesHandles)
        ax = axesHandles(i);
        set(ax, 'Color', axesBackground, 'XColor', 'w', 'YColor', 'w', ...
            'GridColor', gridColor, 'GridAlpha', 0.28, ...
            'MinorGridColor', gridColor, 'MinorGridAlpha', 0.18);
        set([ax.Title, ax.XLabel, ax.YLabel], 'Color', 'w');
    end

    legendHandles = findall(fig, 'Type', 'legend');

    if ~isempty(legendHandles)
        set(legendHandles, 'Color', axesBackground, 'TextColor', 'w', 'EdgeColor', gridColor);
    end
end

function result = numberList(values)
    if isempty(values)
        result = '—';
    else
        result = strtrim(sprintf('%.6g ', values));
    end
end
