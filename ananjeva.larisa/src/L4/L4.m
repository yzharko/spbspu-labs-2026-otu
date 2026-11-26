clear; clc; close all;

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figuresDir = fullfile(projectRoot, 'materials', 'figs');
resultsDir = fullfile(projectRoot, 'materials', 'results');
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

gain = 1;
openCoeffs = [0.05 0.1 1 1 1];
closedCoeffs = openCoeffs; closedCoeffs(end) = closedCoeffs(end) + gain;
order = numel(openCoeffs) - 1;
poleTol = 1e-8;
blue = [0.05 0.25 0.60]; lightBlue = [0.55 0.70 0.85];
red = [0.70 0.10 0.12]; green = [0.10 0.45 0.15]; gray = [0.35 0.35 0.35];

hurwitzOpen = [0.1 1 0 0; 0.05 1 1 0; 0 0.1 1 0; 0 0.05 1 1];
hurwitzClosed = hurwitzOpen; hurwitzClosed(2,3) = 2; hurwitzClosed(4,4) = 2;
minorsOpen = minors(hurwitzOpen);
minorsClosed = minors(hurwitzClosed);
polesOpen = sortpoles(roots(openCoeffs)); polesClosed = sortpoles(roots(closedCoeffs));
rightPoleCountOpen = sum(real(polesOpen) > poleTol); rightPoleCountClosed = sum(real(polesClosed) > poleTol);
stableOpen = double(all(minorsOpen > 0)); stableClosed = double(all(minorsClosed > 0));

coeffs = {openCoeffs, closedCoeffs}; hurwitzMinors = {minorsOpen, minorsClosed}; poles = {polesOpen, polesClosed};
systemTags = {'open', 'closed'}; stability = [stableOpen stableClosed]; rightPoleCounts = [rightPoleCountOpen rightPoleCountClosed];

outputFile = fopen(fullfile(resultsDir, 'hurwitz.csv'), 'w');
fprintf(outputFile, 'system,order,a0,a1,a2,a3,a4,a5,d1,d2,d3,d4,d5,stable,rhp\n');
for idx = 1:2
    coeffRow = [coeffs{idx} NaN(1, 6-numel(coeffs{idx}))];
    minorRow = [hurwitzMinors{idx} NaN(1, 5-numel(hurwitzMinors{idx}))];
    fprintf(outputFile, '%s,%d', systemTags{idx}, order);
    fprintf(outputFile, ',%.6g', coeffRow);
    fprintf(outputFile, ',%.6g', minorRow);
    fprintf(outputFile, ',%d,%d\n', stability(idx), rightPoleCounts(idx));
end
fclose(outputFile);

outputFile = fopen(fullfile(resultsDir, 'roots.csv'), 'w');
fprintf(outputFile, 'system,index,re,im\n');
for idx = 1:2
    for poleIdx = 1:numel(poles{idx})
        re = real(poles{idx}(poleIdx)); im = imag(poles{idx}(poleIdx));
        if abs(re) < 1e-10, re = 0; end
        if abs(im) < 1e-10, im = 0; end
        fprintf(outputFile, '%s,%d,%.6f,%.6f\n', systemTags{idx}, poleIdx, re, im);
    end
end
fclose(outputFile);

freqNyquist = [0, logspace(-2, 3, 4000)];
responseNyquist = polyval(gain, 1i*freqNyquist) ./ polyval(openCoeffs, 1i*freqNyquist);
nyquistDistance = abs(responseNyquist + 1);
[minDistance, minDistanceIdx] = min(nyquistDistance);
closestFreq = freqNyquist(minDistanceIdx);
nyquistRightPoleCount = rightPoleCountOpen;
polesOnImagAxis = sum(abs(real(polesOpen)) <= poleTol);
nyquistAngle = unwrap(angle(responseNyquist + 1));
encirclements = -round((nyquistAngle(end) - nyquistAngle(1))/pi);
closedRightPoleCount = nyquistRightPoleCount + encirclements;

outputFile = fopen(fullfile(resultsDir, 'nyquist.csv'), 'w');
fprintf(outputFile, 'system,P,N,Z,poles_on_axis,omega_star,min_dist\n');
fprintf(outputFile, 'open,%d,%d,%d,%d,%.4f,%.4f\n', nyquistRightPoleCount, encirclements, closedRightPoleCount, polesOnImagAxis, closestFreq, minDistance);
fclose(outputFile);

maxMikhailovFreq = 10*max(abs(polesClosed));
freqMikhailov = linspace(0, maxMikhailovFreq, 20000);
polyMikhailov = polyval(closedCoeffs, 1i*freqMikhailov);
mikhailovAngle = unwrap(angle(polyMikhailov));
angleChangePi = (mikhailovAngle(end) - mikhailovAngle(1))/pi;
curvePassesOrigin = double(any(abs(real(polesClosed)) <= poleTol));

outputFile = fopen(fullfile(resultsDir, 'mikhailov.csv'), 'w');
fprintf(outputFile, 'system,order,arg_pi,rhp,stable,through_origin\n');
fprintf(outputFile, 'closed,%d,%.4f,%d,%d,%d\n', order, angleChangePi, rightPoleCountClosed, stableClosed, curvePassesOrigin);
fclose(outputFile);

crossingIdx = find(real(polyMikhailov(1:end-1)).*real(polyMikhailov(2:end)) < 0, 1);
if isempty(crossingIdx)
    crossingFreq = maxMikhailovFreq; crossingImag = imag(polyMikhailov(end));
else
    crossingFreq = freqMikhailov(crossingIdx); crossingImag = imag(polyMikhailov(crossingIdx));
end

s = tf('s');
openLoopTF = gain/(openCoeffs(1)*s^4 + openCoeffs(2)*s^3 + openCoeffs(3)*s^2 + openCoeffs(4)*s + openCoeffs(5));
systems = {openLoopTF, feedback(openLoopTF, 1)};
sNames = {'Разомкнутая система', 'Замкнутая система'};
time = cell(1,2); response = cell(1,2);
dcGain = zeros(1,2); finalValue = zeros(1,2); maxAbsResponse = zeros(1,2);
diverges = zeros(1,2); oscillates = zeros(1,2); simDuration = zeros(1,2);
for idx = 1:2
    systemPoles = poles{idx};
    diverges(idx) = double(any(real(systemPoles) > poleTol));
    oscillates(idx) = double(~diverges(idx) && any(abs(real(systemPoles)) <= poleTol));
    if diverges(idx)
        simDuration(idx) = max(0.5, 3/max(real(systemPoles)));
    elseif oscillates(idx)
        simDuration(idx) = 40;
    else
        simDuration(idx) = min(40, max(12, ceil(8/min(abs(real(systemPoles))))));
    end
    time{idx} = linspace(0, simDuration(idx), 4000).';
    response{idx} = lsim(systems{idx}, ones(size(time{idx})), time{idx});
    dcGain(idx) = dcgain(systems{idx});
    finalValue(idx) = response{idx}(end);
    maxAbsResponse(idx) = max(abs(response{idx}));
end

outputFile = fopen(fullfile(resultsDir, 'step.csv'), 'w');
fprintf(outputFile, 'system,dc_gain,y_end,max_abs,diverges,oscillates,t_end\n');
for idx = 1:2
    fprintf(outputFile, '%s,%.6g,%.6g,%.6g,%d,%d,%.6g\n', ...
        systemTags{idx}, dcGain(idx), finalValue(idx), maxAbsResponse(idx), diverges(idx), oscillates(idx), simDuration(idx));
end
fclose(outputFile);

figureHandle = canvas(17, 10);
plot(real(responseNyquist), imag(responseNyquist), 'Color', blue, 'LineWidth', 1.4); hold on;
plot(real(responseNyquist), -imag(responseNyquist), '--', 'Color', lightBlue, 'LineWidth', 1);
plot(-1, 0, 'x', 'Color', red, 'MarkerSize', 10, 'LineWidth', 1.8);
idx = round(0.62*numel(freqNyquist));
tangent = [real(responseNyquist(idx+1))-real(responseNyquist(idx)), imag(responseNyquist(idx+1))-imag(responseNyquist(idx))];
tangent = tangent/norm(tangent)*0.22;
quiver(real(responseNyquist(idx)), imag(responseNyquist(idx)), tangent(1), tangent(2), 0, 'Color', gray, ...
    'LineWidth', 1, 'MaxHeadSize', 0.6);
text(real(responseNyquist(idx))+0.05, imag(responseNyquist(idx))+0.06, '\omega', 'Color', gray, 'FontSize', 8);
grid on; box on;
xlim([min([real(responseNyquist), -1])-0.15, max(real(responseNyquist))+0.15]);
ylim([min(imag(responseNyquist))-0.15, max(-imag(responseNyquist))+0.15]);
xlabel('Re W(j\omega)'); ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы', 'FontWeight', 'normal');
legend({'\omega \geq 0', '\omega \leq 0 (зеркально)', 'точка (-1; 0)'}, ...
    'Location', 'northeast'); legend boxoff;
savepng(figureHandle, 'fig01_nyquist_full.png', figuresDir);

nearMinusOne = nyquistDistance <= 1.3*minDistance;
if nnz(nearMinusOne) < 5, nearMinusOne = nyquistDistance <= 2*minDistance; end
figureHandle = canvas(15, 10);
plot(real(responseNyquist), imag(responseNyquist), 'Color', blue, 'LineWidth', 1.2); hold on;
plot(real(responseNyquist(nearMinusOne)), imag(responseNyquist(nearMinusOne)), 'Color', red, 'LineWidth', 1.8);
plot(-1, 0, 'x', 'Color', red, 'MarkerSize', 10, 'LineWidth', 1.8);
plot(real(responseNyquist(minDistanceIdx)), imag(responseNyquist(minDistanceIdx)), 'o', 'Color', green, ...
    'MarkerSize', 5, 'MarkerFaceColor', green);
plot([-1, real(responseNyquist(minDistanceIdx))], [0, imag(responseNyquist(minDistanceIdx))], ':', 'Color', gray);
text(-0.95, 0.14, '(-1; 0)', 'Color', red, 'FontSize', 8);
text(real(responseNyquist(minDistanceIdx))+0.03, imag(responseNyquist(minDistanceIdx))-0.13, ...
    sprintf('|W + 1| = %.3g', minDistance), 'Color', green, 'FontSize', 8);
grid on; box on;
xlim([min([real(responseNyquist(nearMinusOne)), -1])-0.12, max(real(responseNyquist(nearMinusOne)))+0.12]);
ylim([min(imag(responseNyquist(nearMinusOne)))-0.12, max(imag(responseNyquist(nearMinusOne)))+0.12]);
xlabel('Re W(j\omega)'); ylabel('Im W(j\omega)');
title('Окрестность точки (-1; 0)', 'FontWeight', 'normal');
legend({'годограф', 'участок вблизи (-1; 0)', 'точка (-1; 0)', ...
    sprintf('ближайшая точка: \\omega = %.3g рад/с', closestFreq)}, ...
'Location', 'southwest'); legend boxoff;
savepng(figureHandle, 'fig02_nyquist_zoom.png', figuresDir);

figureHandle = canvas(15, 12);
plot(real(polyMikhailov), imag(polyMikhailov), 'Color', blue, 'LineWidth', 1.4); hold on;
plot(real(polyMikhailov(1)), imag(polyMikhailov(1)), 'o', 'Color', green, ...
    'MarkerSize', 5, 'MarkerFaceColor', green);
plot(real(polyMikhailov(end)), imag(polyMikhailov(end)), 's', 'Color', red, ...
    'MarkerSize', 5, 'MarkerFaceColor', red);
plot(0, 0, '+', 'Color', gray, 'MarkerSize', 7, 'LineWidth', 1);
grid on; box on;
realBounds = [min(real(polyMikhailov)), max(real(polyMikhailov))];
imagBounds = [min(imag(polyMikhailov)), max(imag(polyMikhailov))];
xlim([realBounds(1)-0.08*diff(realBounds), realBounds(2)+0.08*diff(realBounds)]);
ylim([imagBounds(1)-0.08*diff(imagBounds), imagBounds(2)+0.08*diff(imagBounds)]);
xlabel('Re D_c(j\omega)'); ylabel('Im D_c(j\omega)');
title('Кривая Михайлова замкнутой системы', 'FontWeight', 'normal');
legend({'кривая Михайлова', '\omega = 0', sprintf('\\omega = %.3g рад/с', maxMikhailovFreq), ...
    'начало координат'}, 'Location', 'northwest'); legend boxoff;
savepng(figureHandle, 'fig03_mikhailov_full.png', figuresDir);

freqZoom = linspace(0, 1.6*crossingFreq, 600);
polyZoom = polyval(closedCoeffs, 1i*freqZoom);
figureHandle = canvas(15, 10);
plot(real(polyZoom), imag(polyZoom), 'Color', blue, 'LineWidth', 1.5); hold on;
plot(0, 0, '+', 'Color', gray, 'MarkerSize', 7, 'LineWidth', 1);
plot(real(polyZoom(1)), imag(polyZoom(1)), 'o', 'Color', green, ...
    'MarkerSize', 5, 'MarkerFaceColor', green);
plot(0, crossingImag, 'o', 'Color', red, 'MarkerSize', 5, 'MarkerFaceColor', red);
plot([0, 0], [0, crossingImag], ':', 'Color', gray);
text(real(polyZoom(1))+0.12, imag(polyZoom(1))+0.06, '\omega = 0', 'Color', green, 'FontSize', 8);
text(0.12, crossingImag-0.05, sprintf('Re D_c = 0, \\omega = %.3g рад/с', crossingFreq), ...
    'Color', red, 'FontSize', 8);
grid on; box on;
xlim([min(real(polyZoom))-0.2, max(real(polyZoom))+0.2]);
ylim([min(imag(polyZoom))-0.12, max(imag(polyZoom))+0.15]);
xlabel('Re D_c(j\omega)'); ylabel('Im D_c(j\omega)');
title('Начальный участок кривой Михайлова', 'FontWeight', 'normal');
legend({'кривая Михайлова', 'начало координат', '\omega = 0', ...
    'первое пересечение оси Im D_c'}, 'Location', 'southwest'); legend boxoff;
savepng(figureHandle, 'fig04_mikhailov_zoom.png', figuresDir);

figureHandle = canvas(19, 7.5);
for idx = 1:2
    subplot(1, 2, idx);
    plot(time{idx}, response{idx}, 'Color', blue, 'LineWidth', 1.4); hold on;
    yline(dcGain(idx), ':', 'Color', gray, 'LineWidth', 0.9);
    grid on; box on; xlim([0, simDuration(idx)]);
    xlabel('t, с'); ylabel('y(t)');
    title(sNames{idx}, 'FontWeight', 'normal');
end
savepng(figureHandle, 'fig05_step.png', figuresDir);

figureHandle = canvas(17, 9);
plot(time{1}, response{1}, 'Color', blue, 'LineWidth', 1.4); hold on;
plot(time{2}, response{2}, 'Color', red, 'LineWidth', 1.4);
yline(dcGain(1), ':', 'Color', blue, 'LineWidth', 0.9);
yline(dcGain(2), ':', 'Color', red, 'LineWidth', 0.9);
grid on; box on;
xlabel('t, с'); ylabel('y(t)');
title('Реакции на ступенчатое воздействие', 'FontWeight', 'normal');
legend(sNames{:}, 'Location', 'east'); legend boxoff;
savepng(figureHandle, 'fig06_step_compare.png', figuresDir);

fprintf('Миноры: %s | %s\n', sprintf('%.4g ', minorsOpen), sprintf('%.4g ', minorsClosed));
fprintf('P=%d N=%d Z=%d, darg/pi=%.4f\n', nyquistRightPoleCount, encirclements, closedRightPoleCount, angleChangePi);
fprintf('Готово: %s\n', figuresDir);

function principalMinors = minors(hurwitzMatrix)
    idx = size(hurwitzMatrix, 1);
    principalMinors = zeros(1, idx);
    for minorIdx = 1:idx
        principalMinors(minorIdx) = det(hurwitzMatrix(1:minorIdx, 1:minorIdx));
    end
end

function poles = sortpoles(poles)
    [~, idx] = sortrows([-real(poles), -abs(imag(poles))]);
    poles = poles(idx);
end

function figureHandle = canvas(widthCm, heightCm)
    figureHandle = figure('Visible', 'off', 'Color', 'white', 'Units', 'centimeters', ...
        'Position', [1 1 widthCm heightCm]);
    set(figureHandle, 'DefaultAxesFontName', 'Times New Roman', ...
        'DefaultTextFontName', 'Times New Roman', ...
        'DefaultAxesFontSize', 8, 'DefaultTextFontSize', 8);
end

function savepng(figureHandle, fileName, folderPath)
    exportgraphics(figureHandle, fullfile(folderPath, fileName), 'Resolution', 300);
    close(figureHandle);
end
