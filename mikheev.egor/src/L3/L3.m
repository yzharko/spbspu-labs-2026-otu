clear; clc; close all;

K = 2.8;
T = 0.6;
T1 = 0.6;
T2 = 0.5;
T1a = 2.2*T2;

base   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(base, 'materials', 'figs');
resdir = fullfile(base, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end

s = tf('s');

link(1).name = 'Интегрирующее';
link(2).name = 'Дифференцирующее';
link(3).name = 'Усилительное';
link(4).name = 'Апериодическое 1-го порядка';
link(5).name = 'Апериодическое 2-го порядка';
link(6).name = 'Колебательное';
link(7).name = 'Консервативное';
link(8).name = 'Инерционное интегрирующее';
link(9).name = 'Реальное дифференцирующее';
link(10).name = 'Форсирующее 1-го порядка';
link(11).name = 'Изодромное';

link(1).W  = K/s;
link(2).W  = K*s;
link(3).W  = tf(K, 1);
link(4).W  = K/(T*s + 1);
link(5).W  = K/(T2^2*s^2 + T1a*s + 1);
link(6).W  = K/(T2^2*s^2 + T1*s + 1);
link(7).W  = K/(T*s^2 + 1);
link(8).W  = K/(s*(T*s + 1));
link(9).W  = K*s/(T*s + 1);
link(10).W = K*(T*s + 1);
link(11).W = K*(T*s + 1)/s;

w0  = 1/sqrt(T);
w   = unique([logspace(-2, 2.5, 1200), linspace(0.97*w0, 1.03*w0, 60)]);
w(abs(w - w0) < 1e-6) = [];
gap = abs(w - w0) < 0.01*w0;

for k = 1:11
    [m, p] = bode(link(k).W, w);
    link(k).mag   = squeeze(m).';
    link(k).phase = squeeze(p).';
    link(k).L     = 20*log10(link(k).mag);
    h             = squeeze(freqresp(link(k).W, w)).';
    link(k).re    = real(h);
    link(k).im    = imag(h);
end
link(7).L(gap) = NaN;

for k = 1:11
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 18 6.4]);
    set(f,'DefaultAxesFontName','Liberation Serif','DefaultTextFontName', ...
          'Liberation Serif','DefaultAxesFontSize',8);
    t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    switch k
        case 3
            plot(link(k).re, link(k).im, 'o', 'Color', [0.35 0.10 0.55], ...
                 'MarkerFaceColor', [0.35 0.10 0.55], 'MarkerSize', 4);
            text(K*1.03, 0.25, sprintf('(%g; 0)', K), 'FontSize', 7);
            leg = {};
        case 7
            rr = link(k).re;  rr(gap) = NaN;
            plot(rr, zeros(size(rr)), 'Color', [0.35 0.10 0.55], 'LineWidth', 1.4);
            leg = {'\omega > 0 и \omega < 0'};
        otherwise
            plot(link(k).re, link(k).im, 'Color', [0.35 0.10 0.55], ...
                 'LineWidth', 1.4);  hold on;
            plot(link(k).re, -link(k).im, '--', 'Color', [0.72 0.62 0.82], ...
                 'LineWidth', 0.9);
            leg = {'\omega > 0', '\omega < 0 (зеркально)'};
    end
    r = min(max(1.25*max([abs(link(k).re) abs(link(k).im)]), 1.2*K), 4*K);
    if k == 3 || k == 7, r = 1.6*K; end
    axis equal;  xlim([-r r]);  ylim([-r r]);
    grid on;  box on;
    xlabel('Re W(j\omega)');  ylabel('Im W(j\omega)');  title('АФЧХ');
    if ~isempty(leg), legend(leg, 'Location', 'northwest');  legend boxoff; end

    nexttile;
    semilogx(w, link(k).L, 'Color', [0.35 0.10 0.55], 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    lf = link(k).L(isfinite(link(k).L));
    if max(lf) - min(lf) < 5, ylim(mean(lf) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('L, дБ');  title('ЛАЧХ');

    nexttile;
    semilogx(w, link(k).phase, 'Color', [0.35 0.10 0.55], 'LineWidth', 1.4);
    set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    pf = link(k).phase(isfinite(link(k).phase));
    if max(pf) - min(pf) < 5, ylim(mean(pf) + [-10 10]); end
    xlabel('\omega, рад/с');  ylabel('\psi, град');  title('ЛФЧХ');

    title(t, sprintf('№%d. %s', k, link(k).name), 'FontSize', 10, ...
          'FontName', 'Liberation Serif');
    exportgraphics(f, fullfile(figdir, sprintf('fig%02d_link.png', k)), ...
                   'Resolution', 300);
    close(f);
end

Kset = unique([0.5 1 K 2 5]);
colK = lines(numel(Kset));
for mode = 1:2
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 20 14]);
    set(f,'DefaultAxesFontName','Liberation Serif','DefaultTextFontName', ...
          'Liberation Serif','DefaultAxesFontSize',7);
    t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:11
        nexttile;  hold on;
        for j = 1:numel(Kset)
            [m, p] = bode(link(k).W*(Kset(j)/K), w);
            if mode == 1, y = 20*log10(squeeze(m).'); else, y = squeeze(p).'; end
            if k == 7, y(gap) = NaN; end
            semilogx(w, y, 'Color', colK(j,:), 'LineWidth', 1);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        xlabel('\omega, рад/с');
        if mode == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d', k), 'FontSize', 8);
    end
    nexttile(12);  axis off;  hold on;
    for j = 1:numel(Kset), plot(NaN, NaN, 'Color', colK(j,:), 'LineWidth', 1); end
    legend(arrayfun(@(v) sprintf('K = %.3g', v), Kset, 'UniformOutput', false), ...
           'Location', 'west');
    if mode == 1
        title(t, 'ЛАЧХ при разных K', 'FontName', 'Liberation Serif');
        fname = 'fig12_k_lach.png';
    else
        title(t, 'ЛФЧХ при разных K', 'FontName', 'Liberation Serif');
        fname = 'fig13_k_lfch.png';
    end
    exportgraphics(f, fullfile(figdir, fname), 'Resolution', 300);
    close(f);
end

Tset = [0.15 0.3 0.6 1.2 3];
idxT = [4 7 8 9 10 11];
for mode = 1:2
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 20 13]);
    set(f,'DefaultAxesFontName','Liberation Serif','DefaultTextFontName', ...
          'Liberation Serif','DefaultAxesFontSize',7);
    t = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:numel(idxT)
        k = idxT(p);
        nexttile;  hold on;
        for j = 1:numel(Tset)
            sys = make_link(k, Tset(j), T1, T2, T1a, K, s);
            [m, ph] = bode(sys, w);
            if mode == 1, y = 20*log10(squeeze(m).'); else, y = squeeze(ph).'; end
            if k == 7, y(gap) = NaN; end
            semilogx(w, y, 'Color', colK(j,:), 'LineWidth', 1.05);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        if mode == 1 && k == 7, ylim([-60 60]); end
        xlabel('\omega, рад/с');
        if mode == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d', k), 'FontSize', 8);
    end
    nexttile(7);  axis off;  hold on;
    for j = 1:numel(Tset), plot(NaN, NaN, 'Color', colK(j,:), 'LineWidth', 1.05); end
    legend(arrayfun(@(v) sprintf('T = %.3g с', v), Tset, 'UniformOutput', false), ...
           'Location', 'west', 'FontSize', 7);
    if mode == 1
        title(t, 'ЛАЧХ при разных T', 'FontName', 'Liberation Serif');
        fname = 'fig14_t_lach.png';
    else
        title(t, 'ЛФЧХ при разных T', 'FontName', 'Liberation Serif');
        fname = 'fig15_t_lfch.png';
    end
    exportgraphics(f, fullfile(figdir, fname), 'Resolution', 300);
    close(f);
end

panel(1).num = 5;  panel(1).par = 'T1'; panel(1).val = [1.0 1.1 1.4 1.8 2.4];
panel(2).num = 5;  panel(2).par = 'T2'; panel(2).val = [0.15 0.25 0.35 0.5];
panel(3).num = 6;  panel(3).par = 'T1'; panel(3).val = [0.1 0.3 0.6 0.8 0.95];
panel(4).num = 6;  panel(4).par = 'T2'; panel(4).val = [0.35 0.5 0.8 1.2 1.8];
for mode = 1:2
    f = figure('Visible','off','Color','white','Units','centimeters', ...
               'Position',[1 1 17 12]);
    set(f,'DefaultAxesFontName','Liberation Serif','DefaultTextFontName', ...
          'Liberation Serif','DefaultAxesFontSize',8);
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for p = 1:4
        nexttile;  hold on;
        for j = 1:numel(panel(p).val)
            if strcmp(panel(p).par, 'T1')
                sys = make_link(panel(p).num, T, panel(p).val(j), T2, T1a, K, s);
            else
                sys = make_link(panel(p).num, T, T1, panel(p).val(j), T1a, K, s);
            end
            [m, ph] = bode(sys, w);
            if mode == 1, y = 20*log10(squeeze(m).'); else, y = squeeze(ph).'; end
            semilogx(w, y, 'Color', colK(j,:), 'LineWidth', 1.1);
        end
        set(gca, 'XScale', 'log');  grid on;  box on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');
        xlabel('\omega, рад/с');
        if mode == 1, ylabel('L, дБ'); else, ylabel('\psi, град'); end
        title(sprintf('№%d: %s', panel(p).num, panel(p).par));
        legend(arrayfun(@(v) sprintf('%s = %.3g с', panel(p).par, v), ...
               panel(p).val, 'UniformOutput', false), ...
               'Location', 'southwest', 'FontSize', 5);
        legend boxoff;
    end
    if mode == 1
        title(t, 'ЛАЧХ при разных T1 и T2', 'FontName', 'Liberation Serif');
        fname = 'fig16_s_lach.png';
    else
        title(t, 'ЛФЧХ при разных T1 и T2', 'FontName', 'Liberation Serif');
        fname = 'fig17_s_lfch.png';
    end
    exportgraphics(f, fullfile(figdir, fname), 'Resolution', 300);
    close(f);
end

L0   = 20*log10(K);
dsc  = sqrt(T1a^2 - 4*T2^2);
zeta = T1/(2*T2);
Amax = K/(2*zeta*sqrt(1 - zeta^2));

wc  = [K, 1/K, NaN, 1/T, 1/((T1a + dsc)/2), 1/T2, 1/sqrt(T), 1/T, 1/T, 1/T, 1/T];
wc2 = [NaN(1,4), 1/((T1a - dsc)/2), NaN(1,6)];
sl  = [-20 20 0 0 0 0 0 -20 20 0 -20];
sh  = [-20 20 0 -20 -40 -40 -40 -40 0 20 0];
pl  = [-90 90 0 0 0 0 0 -90 90 0 -90];
ph_ = [-90 90 0 -90 -180 -180 -180 -180 0 90 0];
amp = [NaN(1,5), Amax, NaN(1,5)];
wpp = [NaN(1,5), sqrt(1 - 2*zeta^2)/T2, NaN(1,5)];
zta = [NaN(1,5), zeta, NaN(1,5)];
gex = [NaN(1,8), K/T, NaN, K*T];

fid = fopen(fullfile(resdir, 'links.csv'), 'w');
fprintf(fid, 'num,name,omega_c,omega_c2,L0,slope_lo,slope_hi,phi_lo,phi_hi,Amax,omega_p,zeta,gain_extra\n');
M = [wc; wc2; repmat(L0,1,11); sl; sh; pl; ph_; amp; wpp; zta; gex].';
for k = 1:11
    fprintf(fid, '%d,%s,%s\n', k, link(k).name, strjoin(compose('%.6g', M(k,:)), ','));
end
fclose(fid);

fid = fopen(fullfile(resdir, 'points.csv'), 'w');
fprintf(fid, 'num,K,omega_0db\n');
for j = 1:numel(Kset)
    for k = 1:11
        [m, ~] = bode(link(k).W*(Kset(j)/K), w);
        Lk = 20*log10(squeeze(m).');
        sc = sign(Lk);
        q  = find(sc(1:end-1).*sc(2:end) < 0, 1);
        om = NaN;
        if ~isempty(q)
            om = 10^(log10(w(q)) + (-Lk(q))*(log10(w(q+1)) - log10(w(q)))/(Lk(q+1) - Lk(q)));
        end
        fprintf(fid, '%d,%.6g,%.6g\n', k, Kset(j), om);
    end
end
fclose(fid);

function sys = make_link(num, T, T1, T2, T1a, K, s)
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
