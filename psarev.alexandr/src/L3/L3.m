clc; clear; close all;
s = tf('s');

K = 2.6;  T = 0.3;  T1 = 0.6;  T2 = 0.5;
T1a = 1.5;

names = {'W1 интегрирующее','W2 дифференцирующее','W3 усилительное', ...
         'W4 апериодическое 1-го пор.','W5 апериодическое 2-го пор.', ...
         'W6 колебательное','W7 консервативное','W8 инерц. интегрирующее', ...
         'W9 реальное дифф.','W10 форсирующее','W11 изодромное'};

W = {K/s, K*s, tf(K), K/(T*s+1), K/(T2^2*s^2+T1a*s+1), K/(T2^2*s^2+T1*s+1), ...
     K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), K*(T*s+1), K*(T*s+1)/s};

omega_dec_min = -2;
omega_dec_max = 2;
num_freq_points = 2000;
w = logspace(omega_dec_min, omega_dec_max, num_freq_points);

links_num = numel(W);

fam0 = struct('tfs', {}, 'labels', {}, 'title', {});
for n = 1:links_num
    fam0(n).tfs = W(n);  fam0(n).labels = {};  fam0(n).title = names{n};
end
drawAll(fam0, w, 'nyq', 'АФЧХ звеньев');
drawAll(fam0, w, 'L', 'ЛАЧХ звеньев');
drawAll(fam0, w, 'phi', 'ЛФЧХ звеньев');

K_values = [0.5 1 2 2.6 5];
famK = struct('tfs', {}, 'labels', {}, 'title', {});
for n = 1:links_num
    famK(n).tfs = cellfun(@(k) W{n}*k/K, num2cell(K_values), 'UniformOutput', false);
    famK(n).labels = compose('K = %.2g', K_values);
    famK(n).title = names{n};
end
drawAll(famK, w, 'L', 'Влияние K на ЛАЧХ');
drawAll(famK, w, 'phi', 'Влияние K на ЛФЧХ');

T_values = [0.1 0.3 1 2 5];
T1_5 = [1 1.5 2 3 5];
T2_5 = [0.1 0.25 0.5 0.6 0.75];
T1_6 = [0.1 0.3 0.6 0.8 0.95];
T2_6 = [0.4 0.5 0.7 1 2];

panels = { ...
    'W4 (T)', T_values, 'T',  @(p) K/(p*s+1); ...
    'W5 (T1)', T1_5, 'T1', @(p) K/(T2^2*s^2+p*s+1); ...
    'W5 (T2)', T2_5, 'T2', @(p) K/(p^2*s^2+T1a*s+1); ...
    'W6 (T1)', T1_6, 'T1', @(p) K/(T2^2*s^2+p*s+1); ...
    'W6 (T2)', T2_6, 'T2', @(p) K/(p^2*s^2+T1*s+1); ...
    'W7 (T)', T_values, 'T',  @(p) K/(p*s^2+1); ...
    'W8 (T)', T_values, 'T',  @(p) K/(s*(p*s+1)); ...
    'W9 (T)', T_values, 'T',  @(p) K*s/(p*s+1); ...
    'W10 (T)', T_values, 'T',  @(p) K*(p*s+1); ...
    'W11 (T)', T_values, 'T',  @(p) K*(p*s+1)/s};

famT = struct('tfs', {}, 'labels', {}, 'title', {});
for n = 1:size(panels,1)
    famT(n).tfs = cellfun(panels{n,4}, num2cell(panels{n,2}), 'UniformOutput', false);
    famT(n).labels = compose([panels{n,3} ' = %.2g'], panels{n,2});
    famT(n).title = panels{n,1};
end
drawAll(famT, w, 'L', 'Влияние постоянных времени на ЛАЧХ');
drawAll(famT, w, 'phi', 'Влияние постоянных времени на ЛФЧХ');

function drawAll(fam, w, mode, name)
    figure('Name', name, 'NumberTitle', 'off');
    n_rows = 3;
    n_cols = 4;
    gain_link_idx = 3;

    for n = 1:numel(fam)
        subplot(n_rows, n_cols, n); hold on;
        for i = 1:numel(fam(n).tfs)
            if strcmp(mode, 'nyq')
                [re, im] = nyquist(fam(n).tfs{i}, w);
                re = squeeze(re); im = squeeze(im);
                if n == gain_link_idx
                    plot(re(1), im(1), 'o', 'MarkerFaceColor', 'b');
                else
                    h = plot(re, im);
                    plot(re, -im, '--', 'Color', h.Color);
                end
            else
                [mag, ph] = bode(fam(n).tfs{i}, w);
                if strcmp(mode, 'L'), semilogx(w, 20*log10(squeeze(mag)));
                else,                 semilogx(w, squeeze(ph)); end
            end
        end
        if ~strcmp(mode, 'nyq'), set(gca, 'XScale', 'log'); end
        if strcmp(mode, 'nyq')
            switch n
                case {1, 2}, xlim([-1 1]);
                case 3, axis([0 5 -1 1]);
                case 7, axis([-10 10 -1 1]);
                case 10, xlim([0 5]);
                case 11, xlim([0 2]);
            end
        end
        grid on; title(fam(n).title);
        if strcmp(mode, 'nyq'), xlabel('Re'); ylabel('Im');
        else, xlabel('\omega, рад/с'); ylabel(ternary(strcmp(mode,'L'), 'L, дБ', '\psi, °')); end
        if ~isempty(fam(n).labels), legend(fam(n).labels, 'Location', 'best', 'FontSize', 7); end
    end
    sgtitle(name);
end

function r = ternary(c, a, b)
    if c, r = a; else, r = b; end
end
