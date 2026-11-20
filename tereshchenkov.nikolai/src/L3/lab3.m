clc;
clear;
close all;

K = 1.1;
T = 0.2;
T1 = 0.3;
T2 = 0.5;
T1ap = 2 * T2;

s = tf('s');
w = logspace(-2, 2, 2000);

names = {
    'Интегрирующее (идеальное)'
    'Дифференцирующее (идеальное)'
    'Усилительное (безынерционное)'
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
    K/s
    K*s
    tf(K)
    K/(T*s + 1)
    K/(T2^2*s^2 + T1ap*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s
};

rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir)
    rootDir = pwd;
end
figDir = fullfile(rootDir, 'fig');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

Klist = [0.55 1.1 2.2];
Tlist = round(T * [0.2 0.5 1 2 5], 10);

fid = 1;
fprintf(fid, 'K = %g\nT = %g\nT1 = %g\nT2 = %g\nT1ap = %g\n', K, T, T1, T2, T1ap);
fprintf(fid, '20lgK = %.4f\n', 20*log10(K));
fprintf(fid, 'wc4 = %.6f\nwc5 = %.6f\nwn6 = %.6f\nwc7 = %.6f\n', 1/T, 1/T2, 1/T2, 1/sqrt(T));

for i = 1:numel(W)
    saveBode(fullfile(figDir, sprintf('bode_%02d.png', i)), W{i}, w, names{i});
    saveNyquist(fullfile(figDir, sprintf('nyq_%02d.png', i)), W{i}, w, names{i}, i, T, K);

    [mag, phase] = bode(W{i}, w);
    L = 20*log10(squeeze(mag));
    phase = squeeze(phase);
    A1 = abs(evalfr(W{i}, 1i));
    fprintf(fid, '\n=== %d. %s ===\n', i, names{i});
    fprintf(fid, 'A(1) = %.6f  L(1) = %.4f  phi(1) = %.4f\n', A1, 20*log10(A1), interpLog(w, phase, 1));
    fprintf(fid, 'L0dB: %s\n', fmtList(zeroCross(w, L)));

    WK = cell(1, numel(Klist));
    labK = cell(1, numel(Klist));
    for j = 1:numel(Klist)
        WK{j} = (Klist(j) / K) * W{i};
        labK{j} = sprintf('K = %g', Klist(j));
        [mK, ~] = bode(WK{j}, w);
        fprintf(fid, 'K=%g L0dB: %s\n', Klist(j), fmtList(zeroCross(w, 20*log10(squeeze(mK)))));
    end
    saveBode(fullfile(figDir, sprintf('k_%02d.png', i)), WK, w, [names{i} ', K'], labK);
end

% постоянные времени: K фиксирован
saveTime(fid, figDir, 4, 'T', Tlist, @(Tv) K/(Tv*s + 1), w, names{4});
saveTime(fid, figDir, 7, 'T', Tlist, @(Tv) K/(Tv*s^2 + 1), w, names{7});
saveTime(fid, figDir, 8, 'T', Tlist, @(Tv) K/(s*(Tv*s + 1)), w, names{8});
saveTime(fid, figDir, 9, 'T', Tlist, @(Tv) K*s/(Tv*s + 1), w, names{9});
saveTime(fid, figDir, 10, 'T', Tlist, @(Tv) K*(Tv*s + 1), w, names{10});
saveTime(fid, figDir, 11, 'T', Tlist, @(Tv) K*(Tv*s + 1)/s, w, names{11});

T1list5 = [1 1.5 2 4];
saveTime(fid, figDir, 5, 'T1', T1list5, @(Tv) K/(T2^2*s^2 + Tv*s + 1), w, [names{5} ', T2 = 0.5']);
T2list5 = [0.25 0.5 1];
saveTime(fid, figDir, 51, 'T2', T2list5, @(Tv) K/(Tv^2*s^2 + (2*Tv)*s + 1), w, [names{5} ', T1 = 2 T2']);

T1list6 = [0.9 0.6 0.3 0.15 0.05];
saveTime(fid, figDir, 6, 'T1', T1list6, @(Tv) K/(T2^2*s^2 + Tv*s + 1), w, [names{6} ', T2 = 0.5']);
T2list6 = [0.2 0.35 0.5 0.8 1];
saveTime(fid, figDir, 61, 'T2', T2list6, @(Tv) K/(Tv^2*s^2 + T1*s + 1), w, [names{6} ', T1 = 0.3']);

fprintf(fid, '\ncorners\n');
for Tv = Tlist
    fprintf(fid, 'T=%g wc_Ts+1=%g wc_Ts2+1=%g\n', Tv, 1/Tv, 1/sqrt(Tv));
end
for Tv = T1list5
    disc = Tv^2 - 4*T2^2;
    tau1 = (Tv + sqrt(disc)) / 2;
    tau2 = (Tv - sqrt(disc)) / 2;
    fprintf(fid, 'T1=%g tau1=%g tau2=%g w1=%g w2=%g decades=%g\n', ...
        Tv, tau1, tau2, 1/tau1, 1/tau2, log10((1/tau2)/(1/tau1)));
end
zeta = T1 / (2*T2);
wn = 1/T2;
wr = wn * sqrt(1 - 2*zeta^2);
Mr = 1 / (2*zeta*sqrt(1 - zeta^2));
fprintf(fid, 'osc zeta=%g wn=%g wr=%g Mr=%g Lpeak=%g L0=%g\n', ...
    zeta, wn, wr, Mr, 20*log10(K*Mr), 20*log10(K));
for Tv = T1list6
    z = Tv / (2*T2);
    fprintf(fid, 'T1=%g zeta=%g wr=%g MrdB=%g\n', Tv, z, wn*sqrt(max(0, 1 - 2*z^2)), 20*log10(1/(2*z*sqrt(1 - z^2))));
end

fprintf(fid, 'DONE\n');

function saveBode(path, sys, w, name, labels)
    fig = figure('Color', 'w', 'Position', [80 80 720 480]);
    useLight(fig);
    if iscell(sys)
        bode(sys{:}, w);
        grid on;
        legend(labels, 'Location', 'best');
    else
        bode(sys, w);
        grid on;
    end
    title(name, 'Interpreter', 'none');
    useLight(fig);
    exportgraphics(fig, path, 'Resolution', 130);
    close(fig);
end

function saveNyquist(path, sys, w, name, id, T, K)
    fig = figure('Color', 'w', 'Position', [80 80 640 480]);
    useLight(fig);
    if id == 7
        wr = 1 / sqrt(T);
        ww = w(abs(log(w / wr)) > 0.05);
        [re, im] = nyquist(sys, ww);
        re = squeeze(re);
        im = squeeze(im);
        cut = find(abs(diff(re)) > 5);
        re(cut) = NaN;
        im(cut) = NaN;
        plot(re, im, 'b', re, -im, 'b');
        grid on;
        xlabel('Real Axis');
        ylabel('Imaginary Axis');
    else
        nyquist(sys, w);
        grid on;
        if id == 3
            hold on;
            plot(K, 0, 'bo', 'MarkerFaceColor', 'b');
        end
    end
    title(name, 'Interpreter', 'none');
    useLight(fig);
    if id == 7
        xlim([-12 12]);
        ylim([-1.5 1.5]);
    elseif id == 1 || id == 2
        xlim([-0.35 0.35]);
    elseif id == 3
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
    elseif id == 8
        xlim([-0.5 0.15]);
        ylim([-8 8]);
    elseif id == 10
        xlim([0.4 1.8]);
    elseif id == 11
        xlim([-0.2 0.7]);
    end
    exportgraphics(fig, path, 'Resolution', 130);
    close(fig);
end

function useLight(fig)
    try
        theme(fig, 'light');
    catch
    end
    set(fig, 'Color', 'w');
    ax = findall(fig, 'Type', 'axes');
    set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    set(findall(fig, 'Type', 'text'), 'Color', 'k');
    lg = findall(fig, 'Type', 'legend');
    set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.3 0.3 0.3]);
end

function saveTime(fid, figDir, id, label, vals, makeW, w, name)
    WT = cell(1, numel(vals));
    labs = cell(1, numel(vals));
    fprintf(fid, '\n=== time %d %s ===\n', id, label);
    for j = 1:numel(vals)
        WT{j} = makeW(vals(j));
        labs{j} = sprintf('%s = %g', label, vals(j));
        [mag, ~] = bode(WT{j}, w);
        fprintf(fid, '%s=%g L0dB: %s\n', label, vals(j), fmtList(zeroCross(w, 20*log10(squeeze(mag)))));
    end
    saveBode(fullfile(figDir, sprintf('t_%02d.png', id)), WT, w, name, labs);
end

function y = interpLog(x, v, x0)
    y = interp1(log(x), v, log(x0), 'linear');
end

function z = zeroCross(w, L)
    z = nan(1, numel(L) - 1);
    n = 0;
    for k = 1:numel(L) - 1
        if ~isfinite(L(k)) || ~isfinite(L(k + 1))
            continue;
        end
        if L(k) == 0
            n = n + 1;
            z(n) = w(k);
        elseif L(k) * L(k + 1) < 0
            t = L(k) / (L(k) - L(k + 1));
            n = n + 1;
            z(n) = exp(log(w(k)) + t * (log(w(k + 1)) - log(w(k))));
        end
    end
    z = z(1:n);
end

function s = fmtList(z)
    if isempty(z)
        s = 'none';
        return;
    end
    s = strjoin(compose('%.4f', z), ', ');
end
