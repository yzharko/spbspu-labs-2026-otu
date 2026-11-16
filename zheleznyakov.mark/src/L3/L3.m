clear; clc; close all;

cfg.K = 1.8;  cfg.T = 0.6;  cfg.T1 = 0.9;  cfg.T2 = 1.2;
cfg.T1a = 2.2*cfg.T2;
cfg.font = 'Liberation Serif';
cfg.col  = [0.00 0.45 0.55; 0.85 0.33 0.10; 0.93 0.69 0.13; ...
            0.49 0.18 0.56; 0.47 0.67 0.19];
assert(cfg.T1a >= 2*cfg.T2 && cfg.T1 < 2*cfg.T2, 'параметры звеньев 5 и 6');

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.figdir = fullfile(root, 'materials', 'figs');
cfg.resdir = fullfile(root, 'materials', 'results');
if ~exist(cfg.figdir, 'dir'), mkdir(cfg.figdir); end
if ~exist(cfg.resdir, 'dir'), mkdir(cfg.resdir); end

K = cfg.K;  T = cfg.T;  T1 = cfg.T1;  T2 = cfg.T2;  T1a = cfg.T1a;
s = tf('s');

names = {'Интегрирующее'; 'Дифференцирующее'; 'Усилительное'; ...
         'Апериодическое 1-го порядка'; 'Апериодическое 2-го порядка'; ...
         'Колебательное'; 'Консервативное'; 'Инерционное интегрирующее'; ...
         'Реальное дифференцирующее'; 'Форсирующее 1-го порядка'; 'Изодромное'};
W = {K/s, K*s, tf(K,1), K/(T*s + 1), K/(T2^2*s^2 + T1a*s + 1), ...
     K/(T2^2*s^2 + T1*s + 1), K/(T*s^2 + 1), K/(s*(T*s + 1)), ...
     K*s/(T*s + 1), K*(T*s + 1), K*(T*s + 1)/s};

w0  = 1/sqrt(T);
w   = unique([logspace(-2, 2.5, 1200), linspace(0.97*w0, 1.03*w0, 60)]);
w(abs(w - w0) < 1e-6) = [];
gap = abs(w - w0) < 0.01*w0;

for i = 1:11
    [m, p]   = bode(W{i}, w);
    MAG(i,:) = squeeze(m).';
    PH(i,:)  = squeeze(p).';
    H(i,:)   = squeeze(freqresp(W{i}, w)).';
end
L = 20*log10(MAG);
L(7, gap) = NaN;

for i = 1:11
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 19 6.6]);
    set(f,'DefaultAxesFontName',cfg.font,'DefaultTextFontName',cfg.font, ...
          'DefaultAxesFontSize',8);

    re = real(H(i,:));  im = imag(H(i,:));
    subplot(1,3,1);
    switch i
        case 3
            plot(re, im, 'o', 'Color', cfg.col(1,:), 'MarkerSize', 4, ...
                 'MarkerFaceColor', cfg.col(1,:));
            text(K*1.03, 0.25, sprintf('(%g; 0)', K), 'FontSize', 7);
            leg = {};
        case 7
            rr = re;  rr(gap) = NaN;
            plot(rr, zeros(size(rr)), 'Color', cfg.col(1,:), 'LineWidth', 1.4);
            leg = {'\omega > 0 и \omega < 0'};
        otherwise
            plot(re, im, 'Color', cfg.col(1,:), 'LineWidth', 1.4);  hold on;
            plot(re, -im, '--', 'Color', [0.62 0.82 0.82], 'LineWidth', 0.9);
            leg = {'\omega > 0', '\omega < 0 (зеркально)'};
    end
    r = min(max(1.25*max([abs(re) abs(im)]), 1.2*K), 4*K);
    if i == 3 || i == 7, r = 1.6*K; end
    axis equal;  xlim([-r r]);  ylim([-r r]);
    grid on;  box on;
    xlabel('Re W(j\omega)');  ylabel('Im W(j\omega)');  title('АФЧХ');
    if ~isempty(leg), legend(leg, 'Location', 'northwest');  legend boxoff; end

    subplot(1,3,2);
    semilogx(w, L(i,:), 'Color', cfg.col(1,:), 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    lf = L(i, isfinite(L(i,:)));
    if max(lf) - min(lf) < 5, ylim(mean(lf) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('L, дБ');  title('ЛАЧХ');

    subplot(1,3,3);
    semilogx(w, PH(i,:), 'Color', cfg.col(1,:), 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
    set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    pf = PH(i, isfinite(PH(i,:)));
    if max(pf) - min(pf) < 5, ylim(mean(pf) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('\psi, град');  title('ЛФЧХ');

    sgtitle(sprintf('№%d. %s', i, names{i}), 'FontName', cfg.font, 'FontSize', 9);
    exportgraphics(f, fullfile(cfg.figdir, sprintf('fig%02d_link.png', i)), ...
                   'Resolution', 300);
    close(f);
end

Ks  = unique([0.5 1 K 2 5]);
Ts  = [0.15 0.3 0.6 1.2 3];
idx = [4 7 8 9 10 11];

fam = struct('title', {}, 'sys', {}, 'lbl', {});
for i = 1:11
    fam(i).title = sprintf('№%d', i);
    fam(i).sys   = arrayfun(@(k) W{i}*(k/K), Ks, 'UniformOutput', false);
    fam(i).lbl   = arrayfun(@(k) sprintf('K = %.3g', k), Ks, 'UniformOutput', false);
end
plot_family(fam, w, 'L',   'ЛАЧХ при разных K', 'fig12_k_lach.png', cfg, [3 4], 12);
plot_family(fam, w, 'phi', 'ЛФЧХ при разных K', 'fig13_k_lfch.png', cfg, [3 4], 12);

fam = struct('title', {}, 'sys', {}, 'lbl', {});
for j = 1:numel(idx)
    fam(j).title = sprintf('№%d', idx(j));
    fam(j).sys   = arrayfun(@(t) rebuild(idx(j),K,t,T1,T2,T1a,s), ...
                            Ts, 'UniformOutput', false);
    fam(j).lbl   = arrayfun(@(t) sprintf('T = %.3g с', t), Ts, 'UniformOutput', false);
end
plot_family(fam, w, 'L',   'ЛАЧХ при разных T', 'fig14_t_lach.png', cfg, [2 4], 7);
plot_family(fam, w, 'phi', 'ЛФЧХ при разных T', 'fig15_t_lfch.png', cfg, [2 4], 7);

spec = {5,'T1',[2.4 2.64 3.2 4 5.5];  5,'T2',[0.3 0.6 0.9 1.2];
        6,'T1',[0.18 0.45 0.9 1.35 1.8]; 6,'T2',[0.6 1.2 2.4 3.6 4.8]};
fam = struct('title', {}, 'sys', {}, 'lbl', {});
for j = 1:size(spec,1)
    fam(j).title = sprintf('№%d: %s', spec{j,1}, spec{j,2});
    for m = 1:numel(spec{j,3})
        T1m = T1;  T2m = T2;
        if strcmp(spec{j,2}, 'T1'), T1m = spec{j,3}(m); else, T2m = spec{j,3}(m); end
        fam(j).sys{m} = rebuild(spec{j,1}, K, T, T1m, T2m, T1a, s);
        fam(j).lbl{m} = sprintf('%s = %.3g с', spec{j,2}, spec{j,3}(m));
    end
end
plot_family(fam, w, 'L',   'ЛАЧХ при разных T1 и T2', 'fig16_s_lach.png', cfg, [2 2], 0);
plot_family(fam, w, 'phi', 'ЛФЧХ при разных T1 и T2', 'fig17_s_lfch.png', cfg, [2 2], 0);

L0   = 20*log10(K);
d    = sqrt(T1a^2 - 4*T2^2);
zeta = T1/(2*T2);
Amax = K/(2*zeta*sqrt(1 - zeta^2));

wc  = [K, 1/K, NaN, 1/T, 1/((T1a + d)/2), 1/T2, 1/sqrt(T), 1/T, 1/T, 1/T, 1/T];
wc2 = [NaN(1,4), 1/((T1a - d)/2), NaN(1,6)];
sl  = [-20 20 0 0 0 0 0 -20 20 0 -20];
sh  = [-20 20 0 -20 -40 -40 -40 -40 0 20 0];
pl  = [-90 90 0 0 0 0 0 -90 90 0 -90];
phi = [-90 90 0 -90 -180 -180 -180 -180 0 90 0];
amp = [NaN(1,5), Amax, NaN(1,5)];
wpp = [NaN(1,5), sqrt(1 - 2*zeta^2)/T2, NaN(1,5)];
zta = [NaN(1,5), zeta, NaN(1,5)];
gex = [NaN(1,8), K/T, NaN, K*T];

Tlink = table((1:11).', names, wc.', wc2.', repmat(L0,11,1), sl.', sh.', ...
              pl.', phi.', amp.', wpp.', zta.', gex.', 'VariableNames', ...
              {'num','name','omega_c','omega_c2','L0','slope_lo','slope_hi', ...
               'phi_lo','phi_hi','Amax','omega_p','zeta','gain_extra'});
Tlink{:,3:13} = round(Tlink{:,3:13}, 6);
writetable(Tlink, fullfile(cfg.resdir, 'links.csv'));

rows = zeros(11*numel(Ks), 3);
r = 0;
for j = 1:numel(Ks)
    for i = 1:11
        [m, ~] = bode(W{i}*(Ks(j)/K), w);
        Lk = 20*log10(squeeze(m).');
        sc = sign(Lk);
        k  = find(sc(1:end-1).*sc(2:end) < 0, 1);
        om = NaN;
        if ~isempty(k)
            om = 10^(log10(w(k)) + (-Lk(k))*(log10(w(k+1)) - log10(w(k)))/(Lk(k+1) - Lk(k)));
        end
        r = r + 1;
        rows(r,:) = [i, Ks(j), om];
    end
end
Tpt = table(rows(:,1), rows(:,2), rows(:,3), 'VariableNames', {'num','K','omega_0db'});
writetable(Tpt, fullfile(cfg.resdir, 'points.csv'));


function sys = rebuild(num, K, T, T1, T2, T1a, s)
    switch num
        case 4,  sys = K/(T*s + 1);
        case 5,  sys = K/(T2^2*s^2 + T1a*s + 1);
        case 6,  sys = K/(T2^2*s^2 + T1*s + 1);
        case 7,  sys = K/(T*s^2 + 1);
        case 8,  sys = K/(s*(T*s + 1));
        case 9,  sys = K*s/(T*s + 1);
        case 10, sys = K*(T*s + 1);
        case 11, sys = K*(T*s + 1)/s;
    end
end

function plot_family(fam, w, kind, ttl, fname, cfg, layout, legPos)
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 20 13 + 2*(layout(1) > 2)]);
    set(f,'DefaultAxesFontName',cfg.font,'DefaultTextFontName',cfg.font, ...
          'DefaultAxesFontSize',7);
    for j = 1:numel(fam)
        subplot(layout(1), layout(2), j);  hold on;
        for m = 1:numel(fam(j).sys)
            [mg, ph] = bode(fam(j).sys{m}, w);
            if strcmp(kind, 'L'), y = 20*log10(squeeze(mg).'); else, y = squeeze(ph).'; end
            if strcmp(fam(j).title, '№7')
                y(abs(w - 1/sqrt(cfg.T)) < 0.01/sqrt(cfg.T)) = NaN;
            end
            semilogx(w, y, 'Color', cfg.col(m,:), 'LineWidth', 1.05);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        if strcmp(kind, 'L') && strcmp(fam(j).title, '№7'), ylim([-60 60]); end
        xlabel('\omega, рад/с');
        if strcmp(kind, 'L'), ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(fam(j).title, 'FontSize', 8);
        if legPos == 0
            legend(fam(j).lbl, 'Location', 'southwest', 'FontSize', 5);
            legend boxoff;
        end
    end
    if legPos > 0
        subplot(layout(1), layout(2), legPos);  axis off;  hold on;
        for m = 1:numel(fam(1).lbl)
            plot(NaN, NaN, 'Color', cfg.col(m,:), 'LineWidth', 1.05);
        end
        legend(fam(1).lbl, 'Location', 'west', 'FontSize', 7);
    end
    sgtitle(ttl, 'FontName', cfg.font, 'FontSize', 10);
    exportgraphics(f, fullfile(cfg.figdir, fname), 'Resolution', 300);
    close(f);
end
