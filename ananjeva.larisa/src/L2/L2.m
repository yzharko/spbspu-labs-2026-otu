clear; clc; close all;

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figuresDir = fullfile(projectRoot, 'materials', 'figs');
resultsDir = fullfile(projectRoot, 'materials', 'results');
if ~exist(figuresDir, 'dir'), mkdir(figuresDir); end
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

syms s
W1 = 1/(s + 1);
W2 = sym(2);
W3 = 5/(s^2 + 3*s + 2);
W4 = sym(4);

W = simplify(W3*(1 + W2 + W1*W2) / ...
             (1 + W3*W4*(2 + 2*W2 + W1*W2)));
[numerator, denominator] = numden(W);
numerator = expand(numerator);
denominator = expand(denominator);

fprintf('W(s) = %s\n', char(W));
fprintf('num = %s\n', char(numerator));
fprintf('den = %s\n', char(denominator));

fileId = fopen(fullfile(resultsDir, 'analytic.txt'), 'w');
fprintf(fileId, 'W|%s\n', char(W));
fprintf(fileId, 'num|%s\n', strjoin(string(sym2poly(numerator)), ' '));
fprintf(fileId, 'den|%s\n', strjoin(string(sym2poly(denominator)), ' '));
fclose(fileId);

laplaceVariable = tf('s');
block1 = 1/(laplaceVariable + 1);
block2 = tf(2);
block3 = 5/(laplaceVariable^2 + 3*laplaceVariable + 2);
block4 = tf(4);
block1.InputName = 'e1'; block1.OutputName = 'y1';
block2.InputName = 'e2'; block2.OutputName = 'y2';
block3.InputName = 'e3'; block3.OutputName = 'y';
block4.InputName = 'y';  block4.OutputName = 'y4';
sumBlock1 = sumblk('e1 = u - y4');
sumBlock2 = sumblk('e2 = y1 + e1 - y4');
sumBlock3 = sumblk('e3 = y2 + e1 - y4');
originalSystem = connect(block1, block2, block3, block4, sumBlock1, sumBlock2, sumBlock3, 'u', 'y');
equivalentSystem = tf([15 25], [1 4 125 162]);

t = linspace(0, 10, 2001).';
uStep = ones(size(t));
uImp = zeros(size(t));
uImp(t < 0.01) = 100;
uSine = sin(t);
testInputs = {uStep, uImp, uSine};
testNames = {'step', 'impulse', 'sine'};
outputCurves = cell(numel(testNames), 2);

fileId = fopen(fullfile(resultsDir, 'metrics_ctrl.csv'), 'w');
fprintf(fileId, 'test,max_abs_error,rms_error\n');
for testIndex = 1:numel(testNames)
    outputCurves{testIndex,1} = lsim(originalSystem, testInputs{testIndex}, t);
    outputCurves{testIndex,2} = lsim(equivalentSystem, testInputs{testIndex}, t);
    errorSignal = outputCurves{testIndex,1} - outputCurves{testIndex,2};
    fprintf(fileId, '%s,%.6e,%.6e\n', testNames{testIndex}, max(abs(errorSignal)), sqrt(mean(errorSignal.^2)));
end
fclose(fileId);

for testIndex = 1:numel(testNames)
    figureHandle = figure('Visible', 'off', 'Color', 'white', ...
               'Units', 'centimeters', 'Position', [1 1 16 9]);
    set(figureHandle, 'DefaultAxesFontName', 'Times New Roman', ...
           'DefaultTextFontName', 'Times New Roman', ...
           'DefaultAxesFontSize', 10);
    plot(t, outputCurves{testIndex,1}, 'LineWidth', 1.6, 'Color', [0.05 0.25 0.60]); hold on;
    grid on; box on;
    plot(t, outputCurves{testIndex,2}, '--', 'LineWidth', 1.2, 'Color', [0.75 0.20 0.10]);
    xlabel('t, с'); ylabel('y(t)');
    title(['Проверка: ' testNames{testIndex}]);
    legend({'исходная схема', 'эквивалентная W'}, 'Location', 'southeast');
    legend boxoff;
    xlim([0 10]);
    exportgraphics(figureHandle, fullfile(figuresDir, ['ctrl_' testNames{testIndex} '.png']), 'Resolution', 300);
    close(figureHandle);
end

fprintf('\nГотово\n');
