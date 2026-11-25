clear; clc; close all;

K = 1.5;  T = 0.4;  T1 = 0.6;  T2 = 0.8;  T1A = 2.2*T2;

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'materials', 'figs');
resdir = fullfile(root, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end

s = tf('s');
names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
         'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
         'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
         'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};
W = {K/s, K*s, tf(K,1), K/(T*s + 1), K/(T2^2*s^2 + T1A*s + 1), ...
     K/(T2^2*s^2 + T1*s + 1), K/(T*s^2 + 1), K/(s*(T*s + 1)), ...
     K*s/(T*s + 1), K*(T*s + 1), K*(T*s + 1)/s};

resonanceFrequency  = 1/sqrt(T);
wlist   = unique([logspace(-2, 2.5, 1200), linspace(0.97*resonanceFrequency, 1.03*resonanceFrequency, 60)]);
wlist(abs(wlist - resonanceFrequency) < 1e-6) = [];
resonanceGapMask = abs(wlist - resonanceFrequency) < 0.01*resonanceFrequency;

frequencyResponses  = zeros(numel(names), numel(wlist));
phaseResponses = zeros(numel(names), numel(wlist));
for i = 1:numel(names)
    frequencyResponses(i,:)  = squeeze(freqresp(W{i}, wlist)).';
    [~, phaseOutput]  = bode(W{i}, wlist);
    phaseResponses(i,:) = squeeze(phaseOutput).';
end
magnitudeDb = 20*log10(abs(frequencyResponses));
magnitudeDb(7, resonanceGapMask) = NaN;

for i = 1:numel(names)
    fig = figure('Visible','off','Color','white','Units','centimeters', ...
                 'Position',[1 1 19 6.6]);

    realPart = real(frequencyResponses(i,:));  imaginaryPart = imag(frequencyResponses(i,:));
    subplot(1,3,1);
    if i == 3
        plot(realPart, imaginaryPart, 'o', 'Color', [0.05 0.25 0.60], 'MarkerSize', 4, ...
             'MarkerFaceColor', [0.05 0.25 0.60]);
        text(K*1.03, 0.25, sprintf('(%g; 0)', K), 'FontSize', 7);
        legendEntries = {};
    elseif i == 7
        realPartWithoutResonance = realPart;  realPartWithoutResonance(resonanceGapMask) = NaN;
        plot(realPartWithoutResonance, zeros(size(realPartWithoutResonance)), 'Color', [0.05 0.25 0.60], 'LineWidth', 1.4);
        legendEntries = {'\omega > 0 и \omega < 0'};
    else
        plot(realPart, imaginaryPart, 'Color', [0.05 0.25 0.60], 'LineWidth', 1.4);  hold on;
        plot(realPart, -imaginaryPart, '--', 'Color', [0.55 0.70 0.85], 'LineWidth', 0.9);
        legendEntries = {'\omega > 0', '\omega < 0 (зеркально)'};
    end
    axisLimit = min(max(1.25*max([abs(realPart) abs(imaginaryPart)]), 1.2*K), 4*K);
    if i == 3 || i == 7, axisLimit = 1.6*K; end
    axis equal;  xlim([-axisLimit axisLimit]);  ylim([-axisLimit axisLimit]);
    grid on;  box on;
    xlabel('Re W(j\omega)');  ylabel('Im W(j\omega)');  title('АФЧХ');
    if ~isempty(legendEntries), legend(legendEntries, 'Location', 'northwest');  legend boxoff; end

    subplot(1,3,2);
    semilogx(wlist, magnitudeDb(i,:), 'Color', [0.05 0.25 0.60], 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    finiteMagnitudeDb = magnitudeDb(i, isfinite(magnitudeDb(i,:)));
    if max(finiteMagnitudeDb) - min(finiteMagnitudeDb) < 5, ylim(mean(finiteMagnitudeDb) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('L, дБ');  title('ЛАЧХ');

    subplot(1,3,3);
    semilogx(wlist, phaseResponses(i,:), 'Color', [0.05 0.25 0.60], 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    finitePhaseDeg = phaseResponses(i, isfinite(phaseResponses(i,:)));
    if max(finitePhaseDeg) - min(finitePhaseDeg) < 5, ylim(mean(finitePhaseDeg) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('\psi, град');  title('ЛФЧХ');

    sgtitle(sprintf('№%d. %s', i, names{i}), 'FontSize', 9);
    exportgraphics(fig, fullfile(figdir, sprintf('fig%02d_link.png', i)), ...
                   'Resolution', 300);
    close(fig);
end

Kset = unique([0.5 1 K 2 5]);
for plotTypeIndex = 1:2
    fig = figure('Visible','off','Color','white','Units','centimeters', ...
                 'Position',[1 1 20 14]);
    for i = 1:numel(names)
        subplot(3,4,i);  hold on;
        for curveIndex = 1:numel(Kset)
            transferFunction = W{i}*(Kset(curveIndex)/K);
            if plotTypeIndex == 1
                responseValues = 20*log10(abs(squeeze(freqresp(transferFunction, wlist)).'));
            else
                [~, phaseOutput] = bode(transferFunction, wlist);  responseValues = squeeze(phaseOutput).';
            end
            if i == 7, responseValues(resonanceGapMask) = NaN; end
            semilogx(wlist, responseValues, 'LineWidth', 1);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        xlabel('\omega, рад/с');
        if plotTypeIndex == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d', i), 'FontSize', 8);
    end
    subplot(3,4,12);  axis off;  hold on;
    for curveIndex = 1:numel(Kset), plot(NaN, NaN, 'LineWidth', 1); end
    legend(arrayfun(@(currentParameterValues) sprintf('K = %.3g', currentParameterValues), Kset, 'UniformOutput', false), ...
           'Location', 'west');
    if plotTypeIndex == 1
        sgtitle('ЛАЧХ при разных K');  outputFileName = 'fig12_k_lach.png';
    else
        sgtitle('ЛФЧХ при разных K');  outputFileName = 'fig13_k_lfch.png';
    end
    exportgraphics(fig, fullfile(figdir, outputFileName), 'Resolution', 300);
    close(fig);
end

Tset = [0.1 0.2 0.4 0.8 2];
selectedFunctionIndices  = [4 7 8 9 10 11];
for plotTypeIndex = 1:2
    fig = figure('Visible','off','Color','white','Units','centimeters', ...
                 'Position',[1 1 20 13]);
    for curveIndex = 1:numel(selectedFunctionIndices)
        subplot(2,4,curveIndex);  hold on;
        for parameterIndex = 1:numel(Tset)
            transferFunction = swap(selectedFunctionIndices(curveIndex), K, Tset(parameterIndex), T1, T2, T1A, s);
            if plotTypeIndex == 1
                responseValues = 20*log10(abs(squeeze(freqresp(transferFunction, wlist)).'));
            else
                [~, phaseOutput] = bode(transferFunction, wlist);  responseValues = squeeze(phaseOutput).';
            end
            if selectedFunctionIndices(curveIndex) == 7, responseValues(resonanceGapMask) = NaN; end
            semilogx(wlist, responseValues, 'LineWidth', 1);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        if selectedFunctionIndices(curveIndex) == 7 && plotTypeIndex == 1, ylim([-60 60]); end
        xlabel('\omega, рад/с');
        if plotTypeIndex == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d', selectedFunctionIndices(curveIndex)), 'FontSize', 8);
    end
    subplot(2,4,7);  axis off;  hold on;
    for parameterIndex = 1:numel(Tset), plot(NaN, NaN, 'LineWidth', 1); end
    legend(arrayfun(@(currentParameterValues) sprintf('T = %.3g с', currentParameterValues), Tset, ...
           'UniformOutput', false), 'Location', 'west');
    if plotTypeIndex == 1
        sgtitle('ЛАЧХ при разных T');  outputFileName = 'fig14_t_lach.png';
    else
        sgtitle('ЛФЧХ при разных T');  outputFileName = 'fig15_t_lfch.png';
    end
    exportgraphics(fig, fullfile(figdir, outputFileName), 'Resolution', 300);
    close(fig);
end

functionNumbers = [5 5 6 6];
parameterValueSets = {[1.7 1.76 2.2 3 4], [0.2 0.4 0.6 0.8], [0.12 0.3 0.6 0.9 1.2], ...
       [0.4 0.8 1.6 2.4 3.2]};
parameterNames = {'T1', 'T2', 'T1', 'T2'};
for plotTypeIndex = 1:2
    fig = figure('Visible','off','Color','white','Units','centimeters', ...
                 'Position',[1 1 17 12]);
    for curveIndex = 1:4
        subplot(2,2,curveIndex);  hold on;
        currentParameterValues = parameterValueSets{curveIndex};
        for parameterIndex = 1:numel(currentParameterValues)
            if strcmp(parameterNames{curveIndex}, 'T1')
                transferFunction = swap(functionNumbers(curveIndex), K, T, currentParameterValues(parameterIndex), T2, T1A, s);
            else
                transferFunction = swap(functionNumbers(curveIndex), K, T, T1, currentParameterValues(parameterIndex), T1A, s);
            end
            if plotTypeIndex == 1
                responseValues = 20*log10(abs(squeeze(freqresp(transferFunction, wlist)).'));
            else
                [~, phaseOutput] = bode(transferFunction, wlist);  responseValues = squeeze(phaseOutput).';
            end
            semilogx(wlist, responseValues, 'LineWidth', 1.1);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        xlabel('\omega, рад/с');
        if plotTypeIndex == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d: %s', functionNumbers(curveIndex), parameterNames{curveIndex}));
        legend(arrayfun(@(x) sprintf('%s = %.3g с', parameterNames{curveIndex}, x), currentParameterValues, ...
               'UniformOutput', false), 'Location', 'southwest', 'FontSize', 5);
        legend boxoff;
    end
    if plotTypeIndex == 1
        sgtitle('ЛАЧХ при разных T1 и T2');  outputFileName = 'fig16_s_lach.png';
    else
        sgtitle('ЛФЧХ при разных T1 и T2');  outputFileName = 'fig17_s_lfch.png';
    end
    exportgraphics(fig, fullfile(figdir, outputFileName), 'Resolution', 300);
    close(fig);
end

lowFrequencyGainDb   = 20*log10(K);
discriminantRoot  = sqrt(T1A^2 - 4*T2^2);
dampingRatio = T1/(2*T2);
peakAmplitude = K/(2*dampingRatio*sqrt(1 - dampingRatio^2));
peakFrequency   = sqrt(1 - 2*dampingRatio^2)/T2;

cutoffFrequencies  = [K, 1/K, NaN, 1/T, 1/((T1A + discriminantRoot)/2), 1/T2, 1/sqrt(T), 1/T, 1/T, 1/T, 1/T];
secondCutoffFrequencies = [NaN(1,4), 1/((T1A - discriminantRoot)/2), NaN(1,6)];
lowFrequencyMagnitudeSlopes  = [-20 20 0 0 0 0 0 -20 20 0 -20];
highFrequencyMagnitudeSlopes  = [-20 20 0 -20 -40 -40 -40 -40 0 20 0];
lowFrequencyPhaseAngles  = [-90 90 0 0 0 0 0 -90 90 0 -90];
highFrequencyPhaseAngles = [-90 90 0 -90 -180 -180 -180 -180 0 90 0];
peakAmplitudes = [NaN(1,5), peakAmplitude, NaN(1,5)];
peakFrequencies = [NaN(1,5), peakFrequency, NaN(1,5)];
dampingRatios = [NaN(1,5), dampingRatio, NaN(1,5)];
additionalGains = [NaN(1,8), K/T, NaN, K*T];

fileId = fopen(fullfile(resdir, 'links.csv'), 'w');
fprintf(fileId, 'num,name,omega_c,omega_c2,L0,slope_lo,slope_hi,phi_lo,phi_hi,Amax,omega_p,zeta,K_extra\n');
for i = 1:numel(names)
    fprintf(fileId, '%d,%s,%.6g,%.6g,%.6g,%d,%d,%d,%d,%.6g,%.6g,%.6g,%.6g\n', ...
            i, names{i}, cutoffFrequencies(i), secondCutoffFrequencies(i), lowFrequencyGainDb, ...
            lowFrequencyMagnitudeSlopes(i), highFrequencyMagnitudeSlopes(i), ...
            lowFrequencyPhaseAngles(i), highFrequencyPhaseAngles(i), ...
            peakAmplitudes(i), peakFrequencies(i), dampingRatios(i), additionalGains(i));
end
fclose(fileId);

fileId = fopen(fullfile(resdir, 'points.csv'), 'w');
fprintf(fileId, 'num,K,omega_0db\n');
for curveIndex = 1:numel(Kset)
    for i = 1:numel(names)
        frequencyResponseAtGain = squeeze(freqresp(W{i}*(Kset(curveIndex)/K), wlist)).';
        magnitudeDbAtGain = 20*log10(abs(frequencyResponseAtGain));
        signChanges = sign(magnitudeDbAtGain);
        crossingIndex  = find(signChanges(1:end-1).*signChanges(2:end) < 0, 1);
        zeroDbFrequency = NaN;
        if ~isempty(crossingIndex)
            zeroDbFrequency = 10^(log10(wlist(crossingIndex)) + (-magnitudeDbAtGain(crossingIndex))*(log10(wlist(crossingIndex+1)) - log10(wlist(crossingIndex)))/(magnitudeDbAtGain(crossingIndex+1) - magnitudeDbAtGain(crossingIndex)));
        end
        fprintf(fileId, '%d,%.6g,%.6g\n', i, Kset(curveIndex), zeroDbFrequency);
    end
end
fclose(fileId);

function transferFunction = swap(functionNumbers, K, T, T1, T2, T1A, s)
    switch functionNumbers
        case 4,  transferFunction = K/(T*s + 1);
        case 5,  transferFunction = K/(T2^2*s^2 + T1A*s + 1);
        case 6,  transferFunction = K/(T2^2*s^2 + T1*s + 1);
        case 7,  transferFunction = K/(T*s^2 + 1);
        case 8,  transferFunction = K/(s*(T*s + 1));
        case 9,  transferFunction = K*s/(T*s + 1);
        case 10, transferFunction = K*(T*s + 1);
        case 11, transferFunction = K*(T*s + 1)/s;
    end
end
