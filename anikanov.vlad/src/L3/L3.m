K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;
T1Aperiodic = 2.5;

syms s;
syms k tau tau1 tau2 tauAperiodic positive;

omega = logspace(-2, 2, 801).';

W = {
    k / s
    k * s
    k
    k / (tau * s + 1)
    k / (tau2 ^ 2 * s ^ 2 + tauAperiodic * s + 1)
    k / (tau2 ^ 2 * s ^ 2 + tau1 * s + 1)
    k / (tau * s ^ 2 + 1)
    k / (s * (tau * s + 1))
    k * s / (tau * s + 1)
    k * (tau * s + 1)
    k * (tau * s + 1) / s
};
names = {
    'Ideal integrator'
    'Ideal differentiator'
    'Gain'
    'First-order lag'
    'Second-order aperiodic link'
    'Oscillatory link'
    'Conservative link'
    'Inertial integrator'
    'Real differentiator'
    'First-order lead'
    'Isodromic link'
};
models = cell(size(W));

absoluteTolerance = 1e-10;
relativeTolerance = 1e-9;
maximumError = 0;
comparisonCount = 0;
assert(T1Aperiodic >= 2 * T2 && T1 < 2 * T2);

for index = 1:numel(W)
    expression = subs(W{index}, k, K);
    expression = subs(expression, tau, T);
    expression = subs(expression, tau1, T1);
    expression = subs(expression, tau2, T2);
    expression = subs(expression, tauAperiodic, T1Aperiodic);
    responseFunction = matlabFunction(subs(expression, s, 1i * s), 'Vars', s);
    analytical = responseFunction(omega) + zeros(size(omega));

    [numerator, denominator] = numden(expression);
    models{index} = tf(sym2poly(numerator), sym2poly(denominator));
    numerical = squeeze(freqresp(models{index}, omega));
    [magnitude, phase] = bode(models{index}, omega);
    [realPart, imaginaryPart] = nyquist(models{index}, omega);
    magnitude = magnitude(:);
    phase = phase(:);
    realPart = realPart(:);
    imaginaryPart = imaginaryPart(:);

    error = abs(analytical - numerical);
    assert(all(isfinite(analytical)) && all(isfinite(numerical)));
    assert(all(error <= absoluteTolerance + relativeTolerance * abs(analytical)));
    assert(all(abs(magnitude .* exp(1i * deg2rad(phase)) - numerical) ...
        <= absoluteTolerance + relativeTolerance * abs(numerical)));
    assert(all(abs(realPart + 1i * imaginaryPart - numerical) ...
        <= absoluteTolerance + relativeTolerance * abs(numerical)));
    maximumError = max(maximumError, max(error ./ (1 + abs(analytical))));
    comparisonCount = comparisonCount + numel(omega);

    plotOmega = omega;
    if index == 7
        split = find(omega < 1 / sqrt(T), 1, 'last');
        plotOmega = [omega(1:split); NaN; omega(split + 1:end)];
        magnitude = [magnitude(1:split); NaN; magnitude(split + 1:end)];
        phase = [phase(1:split); NaN; phase(split + 1:end)];
        realPart = [realPart(1:split); NaN; realPart(split + 1:end)];
        imaginaryPart = [imaginaryPart(1:split); NaN; imaginaryPart(split + 1:end)];
    end

    figure('Name', sprintf('Link %02d', index), 'NumberTitle', 'off', ...
        'Theme', 'light', 'Position', [100, 100, 1000, 520]);
    layout = tiledlayout(2, 2);
    nexttile(layout, 1);
    semilogx(plotOmega, 20 * log10(magnitude));
    grid on;
    title(names{index});
    xlabel('Frequency, rad/s');
    ylabel('Magnitude, dB');

    nexttile(layout, 3);
    semilogx(plotOmega, phase);
    grid on;
    xlabel('Frequency, rad/s');
    ylabel('Phase, deg');

    nexttile(layout, 2, [2, 1]);
    plot(realPart, imaginaryPart, '-', realPart, -imaginaryPart, '--');
    if index == 3
        hold on;
        plot(K, 0, 'o');
        ylim([-1, 1]);
    end
    if max(realPart) - min(realPart) < 1e-8
        xlim(realPart(1) + K * [-1, 1]);
    end
    if index == 7
        xlim(K * [-3, 3]);
        ylim([-1, 1]);
    end
    grid on;
    title('Nyquist response');
    xlabel('Real part');
    ylabel('Imaginary part');
    legend('Positive frequency', 'Negative frequency', 'Location', 'best');
end

KSeries = K * [0.5, 1, 2];
KLabels = cell(numel(KSeries), 1);
gainShift = 20 * log10(KSeries / K);

for index = 1:numel(models)
    reference = squeeze(freqresp(models{index}, omega));
    for trial = 1:numel(KSeries)
        varied = squeeze(freqresp(models{index} * (KSeries(trial) / K), omega));
        assert(max(abs(20 * log10(abs(varied ./ reference)) - gainShift(trial))) < relativeTolerance);
        assert(max(abs(angle(varied ./ reference))) < relativeTolerance);
    end
end

figure('Name', 'Gain variation', 'NumberTitle', 'off', 'Theme', 'light', 'Position', [100, 100, 1000, 450]);
layout = tiledlayout(1, 2);
for trial = 1:numel(KSeries)
    [magnitude, phase] = bode(models{4} * (KSeries(trial) / K), omega);
    KLabels{trial} = sprintf('K = %.2f', KSeries(trial));

    nexttile(layout, 1);
    semilogx(omega, 20 * log10(magnitude(:)));
    hold on;
    nexttile(layout, 2);
    semilogx(omega, phase(:));
    hold on;
end
nexttile(layout, 1);
yline(0, ':', 'HandleVisibility', 'off');
grid on;
title('First-order lag');
xlabel('Frequency, rad/s');
ylabel('Magnitude, dB');
legend(KLabels, 'Location', 'best');
nexttile(layout, 2);
grid on;
title('First-order lag');
xlabel('Frequency, rad/s');
ylabel('Phase, deg');
legend(KLabels, 'Location', 'best');

variationLinks = [4, 7, 8, 9, 10, 11, 5, 5, 6, 6];
variationSymbols = {
    tau
    tau
    tau
    tau
    tau
    tau
    tauAperiodic
    tau2
    tau1
    tau2
};
variationNames = {
    'T'
    'T'
    'T'
    'T'
    'T'
    'T'
    'T1'
    'T2'
    'T1'
    'T2'
};
variationSeries = {
    T * [0.5, 1, 2]
    T * [0.5, 1, 2]
    T * [0.5, 1, 2]
    T * [0.5, 1, 2]
    T * [0.5, 1, 2]
    T * [0.5, 1, 2]
    [2.3, 2.5, 3.0]
    [0.8, 1.1, 1.2]
    [0.5, 1.0, 1.5]
    [0.8, 1.1, 1.4]
};

for first = 1:2:numel(variationLinks)
    figure('Name', sprintf('Time constants %02d', (first + 1) / 2), ...
        'NumberTitle', 'off', 'Theme', 'light', 'Position', [100, 100, 1000, 650]);
    layout = tiledlayout(2, 2);
    for experiment = first:min(first + 1, numel(variationLinks))
        index = variationLinks(experiment);
        series = variationSeries{experiment};
        labels = cell(numel(series), 1);
        row = experiment - first;

        for trial = 1:numel(series)
            expression = subs(W{index}, variationSymbols{experiment}, series(trial));
            expression = subs(expression, k, K);
            expression = subs(expression, tau, T);
            expression = subs(expression, tau1, T1);
            expression = subs(expression, tau2, T2);
            expression = subs(expression, tauAperiodic, T1Aperiodic);

            [numerator, denominator] = numden(expression);
            model = tf(sym2poly(numerator), sym2poly(denominator));
            [magnitude, phase] = bode(model, omega);
            magnitude = magnitude(:);
            phase = phase(:);

            plotOmega = omega;
            if index == 7
                split = find(omega < 1 / sqrt(series(trial)), 1, 'last');
                plotOmega = [omega(1:split); NaN; omega(split + 1:end)];
                magnitude = [magnitude(1:split); NaN; magnitude(split + 1:end)];
                phase = [phase(1:split); NaN; phase(split + 1:end)];
            end
            labels{trial} = sprintf('%s = %.2f', variationNames{experiment}, series(trial));

            nexttile(layout, 2 * row + 1);
            semilogx(plotOmega, 20 * log10(magnitude));
            hold on;
            nexttile(layout, 2 * row + 2);
            semilogx(plotOmega, phase);
            hold on;
        end

        for column = 1:2
            nexttile(layout, 2 * row + column);
            grid on;
            title(names{index});
            xlabel('Frequency, rad/s');
            if column == 1
                ylabel('Magnitude, dB');
            else
                ylabel('Phase, deg');
            end
            legend(labels, 'Location', 'best');
        end
    end
end
