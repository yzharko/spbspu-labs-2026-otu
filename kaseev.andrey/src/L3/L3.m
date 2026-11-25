clearvars;
close all;
clc;

K = 2.5;
T = 0.8;
T1 = 0.8;
T2 = 1.0;
s = tf('s');

W = {
    K/s
    K*s
    tf(K)
    K/(T*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s
};

assert(T1 < 2*T2, 'The variant parameters no longer match the report.');
assert(isequal(W{5}, W{6}), 'Links 5 and 6 should be identical for this variant.');

outDir = getenv('OTU_L3_OUTPUT_DIR');
if isempty(outDir)
    outDir = fullfile(tempdir, 'OTU_L3_figures');
end
if ~isfolder(outDir)
    mkdir(outDir);
end

omega = logspace(-2, 2, 800);
omega0 = 1/sqrt(T);
crossings = cell(numel(W), 1);

for k = 1:numel(W)
    [mag, ~, wBode] = bode(W{k}, omega);
    [realPart, imagPart, wNyquist] = nyquist(W{k}, omega);
    mag = squeeze(mag);
    response = squeeze(realPart) + 1i*squeeze(imagPart);
    assert(max(abs(wBode(:) - omega(:))) < 1e-9);
    assert(max(abs(wNyquist(:) - omega(:))) < 1e-9);
    valid = isfinite(response) & isfinite(mag) & mag > 0;
    mismatch = max(abs(abs(response(valid)) - mag(valid)) ./ max(1, mag(valid)));
    assert(mismatch < 1e-7, 'Bode and Nyquist magnitude disagree.');

    db = 20*log10(mag);
    signChanges = find(db(1:end-1).*db(2:end) < 0);
    roots = zeros(1, numel(signChanges));
    for j = 1:numel(signChanges)
        q = signChanges(j);
        roots(j) = fzero(@(w) 20*log10(abs(evalfr(W{k}, 1i*w))), ...
            [omega(q), omega(q+1)]);
    end
    crossings{k} = roots;

    f = figure('Visible', 'on', 'Color', 'w', 'Position', [100 100 900 530]);
    h = bodeplot(W{k}, omega);
    stylePlot(f, h);
    drawnow;
    frame = getframe(f);
    imwrite(frame.cdata, fullfile(outDir, sprintf('bode_%02d.png', k)));
    close(f);

    f = figure('Visible', 'on', 'Color', 'w', 'Position', [100 100 900 530]);
    h = nyquistplot(W{k}, omega);
    if k == 3
        h.XLimits = [-1 4];
        h.YLimits = [-1 1];
    elseif k == 7
        h.YLimits = [-0.5 0.5];
    end
    stylePlot(f, h);
    if k == 3
        ax = findall(f, 'Type', 'axes');
        hold(ax(1), 'on');
        plot(ax(1), K, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
        text(ax(1), K+0.12, 0.10, '(2.5, 0)', 'Color', 'k', 'FontSize', 12);
    end
    drawnow;
    frame = getframe(f);
    imwrite(frame.cdata, fullfile(outDir, sprintf('nyquist_%02d.png', k)));
    close(f);
end

Kvalues = [1.25 2.5 5.0];
Tvalues = [0.4 0.8 1.6];
for kind = 1:2
    f = figure('Visible', 'on', 'Color', 'w', 'Position', [100 100 900 530]);
    if kind == 1
        sys1 = Kvalues(1)/(T*s+1); sys1.Name = 'K = 1.25';
        sys2 = Kvalues(2)/(T*s+1); sys2.Name = 'K = 2.5';
        sys3 = Kvalues(3)/(T*s+1); sys3.Name = 'K = 5';
        h = bodeplot(sys1, sys2, sys3, omega);
        fileName = 'bode_K_effect.png';
        h.Title.String = 'Effect of K on the first-order link';
    else
        sys1 = K/(Tvalues(1)*s+1); sys1.Name = 'T = 0.4';
        sys2 = K/(Tvalues(2)*s+1); sys2.Name = 'T = 0.8';
        sys3 = K/(Tvalues(3)*s+1); sys3.Name = 'T = 1.6';
        h = bodeplot(sys1, sys2, sys3, omega);
        fileName = 'bode_T_effect.png';
        h.Title.String = 'Effect of T on the first-order link';
    end
    stylePlot(f, h);
    h.LegendVisible = 'on';
    drawnow;
    legends = findall(f, 'Type', 'legend');
    if ~isempty(legends)
        set(legends, 'Color', 'w', 'TextColor', 'k');
    end
    frame = getframe(f);
    imwrite(frame.cdata, fullfile(outDir, fileName));
    close(f);
end

fid = fopen(fullfile(outDir, 'results.txt'), 'w');
fprintf(fid, 'K = %.2f, T = %.2f, T1 = %.2f, T2 = %.2f\n', K, T, T1, T2);
fprintf(fid, 'T1 >= 2*T2: %d; damping ratio: %.4f\n', T1 >= 2*T2, T1/(2*T2));
fprintf(fid, 'First-order corner frequency: %.6f rad/s\n', 1/T);
fprintf(fid, 'Conservative resonance: %.6f rad/s\n', omega0);
for k = 1:numel(W)
    fprintf(fid, '%02d zero-dB frequencies:', k);
    fprintf(fid, ' %.6f', crossings{k});
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('MATLAB screenshots and results saved in %s\n', outDir);

function stylePlot(f, h)
    h.Title.Color = [0 0 0];
    h.XLabel.Color = [0 0 0];
    set(h.YLabel, 'Color', [0 0 0]);
    drawnow;
    ax = findall(f, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.7 0.7 0.7], 'FontSize', 12);
    set(findall(f, 'Type', 'text'), 'Color', 'k');
end
