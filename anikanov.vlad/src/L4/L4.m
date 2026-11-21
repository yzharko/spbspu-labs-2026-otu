clear;
clc;

K = 3;

W_open = tf(K, [0.1, 0.01, 0.1, 1]);
W = {
    W_open
    feedback(W_open, 1)
};
names = {
    'Open-loop system'
    'Closed-loop system'
};

omega = unique([0, logspace(-3, 3, 6001), 1, 20]);
time = (0:0.005:6).';

tolerance = 1e-8;
denominators = cell(size(W));
rightPoles = zeros(numel(W), 1);
stepError = zeros(numel(W), 1);
stepResponses = zeros(numel(time), numel(W));

for index = 1:numel(W)
    [numerator, denominator] = tfdata(W{index}, 'v');
    denominators{index} = denominator;
    order = numel(denominator) - 1;
    assert(order == 3);
    hurwitz = [
        denominator(2), denominator(4), 0
        denominator(1), denominator(3), 0
        0, denominator(2), denominator(4)
    ];
    minors = zeros(order, 1);
    for minor = 1:order
        minors(minor) = det(hurwitz(1:minor, 1:minor));
    end

    poles = pole(W{index});
    rightPoles(index) = sum(real(poles) > tolerance);
    assert(all(abs(real(poles)) > tolerance));
    assert(all(minors > 0) == isstable(W{index}));
    assert(all(minors > 0) == (rightPoles(index) == 0));

    stepResponses(:, index) = step(W{index}, time);
    [residues, responsePoles, direct] = residue(numerator, conv(denominator, [1, 0]));
    assert(isempty(direct));
    responseCheck = exp(time * responsePoles.') * residues;
    assert(max(abs(imag(responseCheck))) < tolerance);
    responseCheck = real(responseCheck);
    stepError(index) = max(abs(responseCheck - stepResponses(:, index)) ./ (1 + abs(responseCheck)));
    assert(stepError(index) < tolerance);

    disp(names{index});
    disp('Hurwitz matrix:');
    disp(hurwitz);
    disp('Leading Hurwitz minors:');
    disp(minors);
    disp('Poles:');
    disp(poles);
end

expectedDenominator = denominators{1};
expectedDenominator(end) = expectedDenominator(end) + K;
assert(max(abs(denominators{2} - expectedDenominator)) < tolerance);

response = reshape(freqresp(W{1}, omega), 1, []);
fullResponse = [conj(fliplr(response)), response(2:end)];
nyquistPhase = unwrap(angle(1 + fullResponse));
clockwiseTurns = -(nyquistPhase(end) - nyquistPhase(1)) / (2 * pi);
encirclements = round(clockwiseTurns);
assert(abs(clockwiseTurns - encirclements) < 1e-6);
assert(rightPoles(2) == rightPoles(1) + encirclements);

characteristic = polyval(denominators{2}, 1i * omega);
mikhailovPhase = unwrap(angle(characteristic));
phaseChange = mikhailovPhase(end) - mikhailovPhase(1);
assert(abs(phaseChange - (order - 2 * rightPoles(2)) * pi / 2) < 1e-3);

disp('Clockwise Nyquist encirclements:');
disp(encirclements);
disp('Mikhailov phase change divided by pi:');
disp(phaseChange / pi);
disp('Normalized step verification errors:');
disp(stepError);

nyquistFigure = figure('Name', 'Nyquist criterion', 'NumberTitle', 'off', ...
    'Theme', 'light', 'Position', [80, 100, 1150, 500]);
layout = tiledlayout(nyquistFigure, 1, 2);

nexttile(layout);
nyquist(W{1}, omega);
chart = gca;
chart.Title.String = 'Open-loop Nyquist plot';
chart.XLabel.String = 'Re W(j\omega)';
chart.YLabel.String = 'Im W(j\omega)';
chart.XLimits = [-1.3, 3.5];
grid off;

nexttile(layout);
nyquist(W{1}, omega);
chart = gca;
chart.Title.String = 'Critical point region';
chart.XLabel.String = 'Re W(j\omega)';
chart.YLabel.String = 'Im W(j\omega)';
chart.XLimits = [-1.3, 0.2];
chart.YLimits = [-0.5, 0.5];
grid off;

mikhailovFigure = figure('Name', 'Mikhailov criterion', 'NumberTitle', 'off', ...
    'Theme', 'light', 'Position', [80, 100, 1150, 500]);
layout = tiledlayout(mikhailovFigure, 1, 2);
frequencyLimits = [30, 2];
viewNames = {
    'Mikhailov curve'
    'Initial segment'
};

for view = 1:numel(frequencyLimits)
    frequencies = linspace(0, frequencyLimits(view), 3001);
    curve = polyval(denominators{2}, 1i * frequencies);

    ax = nexttile(layout);
    plot(ax, real(curve), imag(curve));
    hold(ax, 'on');
    plot(ax, real(curve(1)), imag(curve(1)), 'o');
    xline(ax, 0, ':');
    yline(ax, 0, ':');
    grid(ax, 'on');
    title(ax, viewNames{view});
    xlabel(ax, 'Re D_c(j\omega)');
    ylabel(ax, 'Im D_c(j\omega)');
    if view == 2
        realLimits = [min(real(curve)), max(real(curve))];
        padding = 0.1 * (realLimits(2) - realLimits(1));
        xlim(ax, realLimits + [-padding, padding]);
    end
end

stepFigure = figure('Name', 'Step responses', 'NumberTitle', 'off', ...
    'Theme', 'light', 'Position', [80, 100, 1150, 500]);
layout = tiledlayout(stepFigure, 1, 2);

for index = 1:numel(W)
    ax = nexttile(layout);
    plot(ax, time, stepResponses(:, index));

    grid(ax, 'on');
    title(ax, names{index});
    xlabel(ax, 't, s');
    ylabel(ax, 'y(t)');
end
