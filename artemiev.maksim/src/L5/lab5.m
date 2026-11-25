clear;
clc;
close all;

%% 1. Исходная система и требования варианта 2
plantCorners = [5 10 20];
plantDen = denominator(plantCorners);
plant = tf(32, plantDen);
initialClosed = feedback(plant, 1);      
requirements = [0.6 1.19 20];   
settlingBand = 0.02;

initialHurwitz = cubic_hurwitz(plantDen + [0 0 0 32]);
[gmInitial, pmInitial, ~, wcInitial] = margin(plant);
fprintf('Исходные полюса замкнутой системы:\n'); disp(pole(initialClosed))
fprintf('Миноры Гурвица: %.6g  %.6g  %.6g\n', initialHurwitz);
fprintf('Исходные запасы: %.3f дБ, %.3f град. Срез: %.3f рад/с.\n', ...
    20*log10(gmInitial), pmInitial, wcInitial);
if ~isstable(initialClosed)
    fprintf('У исходной системы нет установившегося отклика((\n');
end

%% 2. Желаемая петлевая характеристика
desiredCorners = [plantCorners(1) 60 100];
desiredDen = denominator(desiredCorners);
phaseReference = asind(1/requirements(2));
crossoverCandidates = [20 25 30];
fprintf('\nПробные частоты среза (ωс, запас по фазе, Tпп, M, σ):\n');
for candidate = crossoverCandidates
    candidateGain = prod(hypot(1, candidate ./ desiredCorners));
    candidateLoop = tf(candidateGain, desiredDen);
    candidateClosed = feedback(candidateLoop, 1);
    candidateQuality = stepinfo(candidateClosed, 'SettlingTimeThreshold', settlingBand);
    candidatePhase = 180 - sum(atand(candidate ./ desiredCorners));
    fprintf('%5.1f  %6.2f  %.4f  %.4f  %6.2f\n', candidate, ...
        candidatePhase, candidateQuality.SettlingTime, ...
        getPeakGain(candidateClosed), candidateQuality.Overshoot);
end
wc = 25;
desiredGain = prod(hypot(1, wc ./ desiredCorners));
desired = tf(desiredGain, desiredDen);
fprintf('Выбранный срез: %.2f рад/с, ориентир по фазе %.2f град.\n', wc, phaseReference);

%% 3. Геометрическое вычитание ЛАЧХ
w = unique([logspace(-1, 3, 1000), plantCorners, desiredCorners, wc]);
[initialAmplitude, initialPhase] = bode_data(plant, w);
[desiredAmplitude, desiredPhase] = bode_data(desired, w);
initialAsymptote = lowpass_asymptote(32, plantCorners, w);
desiredAsymptote = lowpass_asymptote(desiredGain, desiredCorners, w);
deviceAsymptote = desiredAsymptote - initialAsymptote;

deviceGain = desiredGain / 32;
numeratorCorners = plantCorners(2:3);
denominatorCorners = desiredCorners(2:3);
device = tf(deviceGain * denominator(numeratorCorners), denominator(denominatorCorners));
fprintf('\nЖелаемый коэффициент %.8f, усиление ПКУ %.8f.\n', desiredGain, deviceGain);
fprintf('ПКУ:\n'); disp(device)

%% 4. Проверка скорректированной системы
correctedFull = feedback(device * plant, 1);
correctedOpen = minreal(device * plant, 1e-7);
correctedClosed = minreal(correctedFull, 1e-7);
[gmCorrected, pmCorrected, ~, wcCorrected] = margin(correctedOpen);
[correctedAmplitude, correctedPhase] = bode_data(correctedOpen, w);
[closedAmplitude, closedPhase] = bode_data(correctedClosed, w);
[M, resonance] = getPeakGain(correctedClosed, 1e-8);
yFinal = dcgain(correctedClosed);
t = (0:0.0002:2).';
y = step(correctedClosed, t);
quality = stepinfo(y, t, yFinal, 'SettlingTimeThreshold', settlingBand);

fprintf('Полюса скорректированной системы до сокращения:\n'); disp(pole(correctedFull))
fprintf('Миноры Гурвица минимальной модели: %.8g  %.8g  %.8g\n', ...
    cubic_hurwitz(desiredDen + [0 0 0 desiredGain]));
fprintf('Запасы после коррекции: %.3f дБ, %.3f град. Срез: %.3f рад/с.\n', ...
    20*log10(gmCorrected), pmCorrected, wcCorrected);
fprintf('Tпп = %.5f с; σ = %.4f %%.\n', quality.SettlingTime, quality.Overshoot);
fprintf('M = %.6f при ω = %.4f рад/с; M/y∞ = %.6f.\n', M, resonance, M/yFinal);
fprintf('y∞ = %.6f; e∞ = %.6f.\n', yFinal, 1-yFinal);
fprintf('Устойчивость и соответствие (уст., Tпп, M, σ): %d %d %d %d.\n', ...
    isstable(correctedFull), quality.SettlingTime <= requirements(1), ...
    M <= requirements(2), quality.Overshoot <= requirements(3));

%% 5. Графики
figures = gobjects(5,1);
blue = [0 0.32 0.66];
red = [0.78 0.25 0.12];
green = [0.00 0.47 0.37];

figures(1) = lab_figure('Исходная замкнутая САУ', 390);
ax = axes(figures(1));
tInitial = (0:0.001:1).';
plot(ax, tInitial, step(initialClosed, tInitial), 'Color', blue, 'LineWidth', 1.7);
title(ax, 'Исходная замкнутая САУ');
xlabel(ax, 't, с'); ylabel(ax, 'y(t)'); grid(ax, 'on');

figures(2) = lab_figure('Располагаемая и желаемая ЛЧХ', 550);
tiles = tiledlayout(figures(2), 2, 1, 'TileSpacing','compact','Padding','compact');
ax = nexttile(tiles);
semilogx(ax, w, initialAmplitude, 'Color', blue, 'LineWidth', 1.6);
hold(ax, 'on');
semilogx(ax, w, desiredAmplitude, 'Color', red, 'LineWidth', 1.6);
yline(ax, 0, ':', 'HandleVisibility','off');
xline(ax, wc, ':', '25 рад/с', 'HandleVisibility','off');
ylabel(ax, 'L, дБ'); grid(ax, 'on');
legend(ax, 'Исходная W', 'Желаемая W_{ж}', 'Location','southwest');
ax = nexttile(tiles);
semilogx(ax, w, initialPhase, 'Color', blue, 'LineWidth', 1.6);
hold(ax, 'on');
semilogx(ax, w, desiredPhase, 'Color', red, 'LineWidth', 1.6);
yline(ax, -180, ':');
xline(ax, wc, ':');
xlabel(ax, '\omega, рад/с'); ylabel(ax, '\phi, град'); grid(ax, 'on');

figures(3) = lab_figure('Геометрическое вычитание ЛАЧХ', 450);
ax = axes(figures(3));
semilogx(ax, w, initialAsymptote, 'Color', blue, 'LineWidth', 1.6);
hold(ax, 'on');
semilogx(ax, w, desiredAsymptote, '--', 'Color', red, 'LineWidth', 1.6);
semilogx(ax, w, deviceAsymptote, 'Color', green, 'LineWidth', 1.9);
yline(ax, 0, ':', 'HandleVisibility','off');
xticks(ax, [0.1 1 5 10 20 60 100 1000]);
xlabel(ax, '\omega, рад/с'); ylabel(ax, 'L, дБ'); grid(ax, 'on');
title(ax, 'Асимптотические ЛАЧХ');
legend(ax, 'Исходная', 'Желаемая', 'ПКУ: L_{ж} - L_{исх}', 'Location','southwest');

figures(4) = lab_figure('ЛЧХ скорректированной замкнутой САУ', 550);
tiles = tiledlayout(figures(4), 2, 1, 'TileSpacing','compact','Padding','compact');
ax = nexttile(tiles);
semilogx(ax, w, closedAmplitude, 'Color', green, 'LineWidth', 1.7);
hold(ax, 'on');
plot(ax, resonance, 20*log10(M), 'o', 'Color', red, 'MarkerFaceColor', red);
title(ax, 'Скорректированная замкнутая САУ');
ylabel(ax, 'L, дБ'); grid(ax, 'on');
ax = nexttile(tiles);
semilogx(ax, w, closedPhase, 'Color', green, 'LineWidth', 1.7);
xlabel(ax, '\omega, рад/с'); ylabel(ax, '\phi, град'); grid(ax, 'on');

figures(5) = lab_figure('Скорректированная замкнутая САУ', 430);
ax = axes(figures(5));
plot(ax, t, y, 'Color', green, 'LineWidth', 1.7);
hold(ax, 'on');
yline(ax, yFinal, ':', 'y_\infty');
yline(ax, yFinal*(1+settlingBand), '--');
yline(ax, yFinal*(1-settlingBand), '--');
xline(ax, quality.SettlingTime, ':', sprintf('T_{пп} = %.3f с', quality.SettlingTime));
plot(ax, quality.PeakTime, quality.Peak, 'o', 'Color', red, 'MarkerFaceColor', red);
xlim(ax, [0 0.8]); ylim(ax, [0 1.08]);
title(ax, 'Отклик на единичную ступеньку после коррекции');
xlabel(ax, 't, с'); ylabel(ax, 'y(t)'); grid(ax, 'on');

function den = denominator(corners)
    den = 1;
    for corner = corners
        den = conv(den, [1/corner 1]);
    end
end

function minors = cubic_hurwitz(d)
    second = d(2)*d(3) - d(1)*d(4);
    minors = [d(2), second, d(4)*second];
end

function [amplitude, phase] = bode_data(sys, w)
    [mag, phase] = bode(sys, w);
    amplitude = 20*log10(reshape(mag, 1, []));
    phase = reshape(phase, 1, []);
end

function level = lowpass_asymptote(gain, corners, w)
    level = 20*log10(gain) - 20*sum(log10(max(1, w ./ corners(:))), 1);
end

function fig = lab_figure(name, height)
    fig = figure('Name',name,'NumberTitle','off','Theme','light', ...
        'Color','white','Position',[80 80 850 height]);
end
