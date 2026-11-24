clear; clc; close all;

K0 = 1.9;
T  = 0.7;
T1 = 1.1;
T2 = 1.3;

T1_5 = T1;
if T1_5 < 2*T2
    T1_5 = 3.0;
    fprintf('Звено 5: T1 = %.1f < 2*T2 = %.1f, принято T1 = %.1f с\n', T1, 2*T2, T1_5);
end
T1_6 = T1;
if T1_6 >= 2*T2
    T1_6 = T2;
    fprintf('Звено 6: T1 = %.1f >= 2*T2 = %.1f, принято T1 = %.1f с\n', T1, 2*T2, T1_6);
end

Ta = (T1_5 + sqrt(T1_5^2 - 4*T2^2))/2;
Tb = (T1_5 - sqrt(T1_5^2 - 4*T2^2))/2;

w = logspace(-2, 2, 2000);

syms s
syms K om positive

links = {
    'Интегрирующее',                K/s
    'Идеальное дифференцирующее',   K*s
    'Усилительное',                 K
    'Апериодическое 1-го порядка',  K/(T*s + 1)
    'Апериодическое 2-го порядка',  K/(T2^2*s^2 + T1_5*s + 1)
    'Колебательное',                K/(T2^2*s^2 + T1_6*s + 1)
    'Консервативное',               K/(T*s^2 + 1)
    'Инерционное интегрирующее',    K/(s*(T*s + 1))
    'Реальное дифференцирующее',    K*s/(T*s + 1)
    'Форсирующее 1-го порядка',     K*(T*s + 1)
    'Изодромное',                   K*(T*s + 1)/s
};
names = links(:, 1);
n = numel(names);
titles = compose('%d. %s', (1:n)', string(names));

A_an   = sym(zeros(n, 1));
psi_an = sym(zeros(n, 1));
U_an   = sym(zeros(n, 1));
V_an   = sym(zeros(n, 1));

for i = 1:n
    Wjw       = subs(links{i, 2}, s, 1i*om);
    U_an(i)   = simplify(real(Wjw));
    V_an(i)   = simplify(imag(Wjw));
    A_an(i)   = simplify(abs(Wjw));
    psi_an(i) = simplify(atan2(V_an(i), U_an(i)));

    fprintf('\n%2d. %s\n    W(s)   = %s\n    U(w)   = %s\n    V(w)   = %s\n    A(w)   = %s\n    psi(w) = %s\n', ...
        i, names{i}, char(links{i, 2}), char(U_an(i)), char(V_an(i)), char(A_an(i)), char(psi_an(i)));
end

p = tf('s');
W_tf = @(k) {k/p, k*p, tf(k), k/(T*p + 1), ...
             k/(T2^2*p^2 + T1_5*p + 1), k/(T2^2*p^2 + T1_6*p + 1), ...
             k/(T*p^2 + 1), k/(p*(T*p + 1)), k*p/(T*p + 1), ...
             k*(T*p + 1), k*(T*p + 1)/p};
W0 = W_tf(K0);

fprintf('\nПроверка: аналитика против bode/nyquist (K = %.1f)\n', K0);
for i = 1:n
    [mag, ph] = bode(W0{i}, w);
    [re, im]  = nyquist(W0{i}, w);
    L_num   = 20*log10(squeeze(mag))';
    psi_num = squeeze(ph)';
    UV_num  = squeeze(re)' + 1i*squeeze(im)';

    L_an   = 20*log10(evalOnGrid(subs(A_an(i), K, K0), om, w));
    psi_d  = rad2deg(evalOnGrid(subs(psi_an(i), K, K0), om, w));
    UV_an  = evalOnGrid(subs(U_an(i), K, K0), om, w) + 1i*evalOnGrid(subs(V_an(i), K, K0), om, w);

    errL   = max(abs(L_num - L_an));
    errPsi = max(abs(mod(psi_num - psi_d + 180, 360) - 180));
    errUV  = max(abs(UV_num - UV_an));
    fprintf('%2d. %-30s max|dL| = %.1e дБ, max|dpsi| = %.1e град, max|dW| = %.1e\n', ...
        i, names{i}, errL, errPsi, errUV);
end

asym = {
    {-1, [], [], []}
    {1, [], [], []}
    {0, [], [], []}
    {0, [], T, []}
    {0, [], [Ta Tb], []}
    {0, [], [], T2}
    {0, [], [], sqrt(T)}
    {-1, [], T, []}
    {1, [], T, []}
    {0, T, [], []}
    {-1, T, [], []}
};

tlL   = newFigure('ЛАЧХ типовых звеньев', 3, 4);
tlPsi = newFigure('ЛФЧХ типовых звеньев', 3, 4);
tlAF  = newFigure('АФЧХ типовых звеньев', 3, 4);

for i = 1:n
    [mag, ph] = bode(W0{i}, w);
    L   = 20*log10(squeeze(mag));
    psi = squeeze(ph);

    ax = nexttile(tlL); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    semilogx(ax, w, L, 'b', 'LineWidth', 1.3, 'DisplayName', 'точная');
    semilogx(ax, w, asymptote(w, K0, asym{i}), 'r--', 'LineWidth', 1, 'DisplayName', 'асимптотическая');
    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');
    ylim(ax, limitsL(L));
    title(ax, titles{i}, 'FontSize', 9);
    xlabel(ax, '\omega, рад/с'); ylabel(ax, 'L, дБ');

    ax = nexttile(tlPsi); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    semilogx(ax, w, psi, 'b', 'LineWidth', 1.3);
    set(ax, 'XScale', 'log');
    ylim(ax, [min(psi) - 20, max(psi) + 20]);
    title(ax, titles{i}, 'FontSize', 9);
    xlabel(ax, '\omega, рад/с'); ylabel(ax, '\psi, град');

    ax = nexttile(tlAF); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    plotAFC(ax, W0{i}, w, K0);
    title(ax, titles{i}, 'FontSize', 9);
end
lg = legend(nexttile(tlL, 1)); lg.Layout.Tile = 12;

K_set  = K0*[0.5 1 2];
cols   = lines(numel(K_set));
styles = {'-', '--', ':'};
labelsK = compose('K = %.2f', K_set);

tlKL   = newFigure('Влияние K на ЛАЧХ', 3, 4);
tlKPsi = newFigure('Влияние K на ЛФЧХ', 3, 4);

wc_all = cell(n, numel(K_set));
L1     = zeros(n, numel(K_set));
for i = 1:n
    axL   = nexttile(tlKL);   hold(axL, 'on');   grid(axL, 'on');   box(axL, 'on');
    axPsi = nexttile(tlKPsi); hold(axPsi, 'on'); grid(axPsi, 'on'); box(axPsi, 'on');

    for k = 1:numel(K_set)
        Wk = W_tf(K_set(k));
        [mag, ph] = bode(Wk{i}, w);
        L = 20*log10(squeeze(mag));

        semilogx(axL, w, L, 'Color', cols(k, :), 'LineWidth', 1.3);
        wc_all{i, k} = crossings(w, L);
        semilogx(axL, wc_all{i, k}, zeros(size(wc_all{i, k})), 'o', 'Color', cols(k, :), ...
            'MarkerFaceColor', cols(k, :), 'HandleVisibility', 'off');
        semilogx(axPsi, w, squeeze(ph), styles{k}, 'Color', cols(k, :), 'LineWidth', 2.5 - 0.5*k);

        L1(i, k) = 20*log10(abs(evalfr(Wk{i}, 1i)));
    end

    set([axL axPsi], 'XScale', 'log');
    yline(axL, 0, 'k:', 'HandleVisibility', 'off');
    if i == 3
        ylim(axL, [-5 15]);
    elseif i == 7
        ylim(axL, [-60 60]);
    end
    title(axL, titles{i}, 'FontSize', 9); title(axPsi, titles{i}, 'FontSize', 9);
    xlabel(axL, '\omega, рад/с'); ylabel(axL, 'L, дБ');
    xlabel(axPsi, '\omega, рад/с'); ylabel(axPsi, '\psi, град');
end
lg = legend(nexttile(tlKL, 1), labelsK);   lg.Layout.Tile = 12;
lg = legend(nexttile(tlKPsi, 1), labelsK); lg.Layout.Tile = 12;

fprintf('\nЧастоты среза (пересечение ЛАЧХ с осью 0 дБ), рад/с\n');
fprintf('%-34s', 'Звено'); fprintf('%-18s', labelsK{:}); fprintf('\n');
for i = 1:n
    fprintf('%-34s', titles{i});
    for k = 1:numel(K_set)
        fprintf('%-18s', listStr(wc_all{i, k}));
    end
    fprintf('\n');
end

fprintf('\nПересечение ЛАЧХ с осью ординат (w = 1 рад/с), дБ\n');
fprintf('%-34s', 'Звено'); fprintf('%-18s', labelsK{:}); fprintf('\n');
for i = 1:n
    fprintf('%-34s', titles{i}); fprintf('%-18.2f', L1(i, :)); fprintf('\n');
end
fprintf('Сдвиг ЛАЧХ при K -> 2K: 20lg2 = %.2f дБ\n', 20*log10(2));

m_set = [0.5 1 2];
W_m = @(m) {[], [], [], K0/(m*T*p + 1), ...
            K0/((m*T2)^2*p^2 + m*T1_5*p + 1), K0/((m*T2)^2*p^2 + m*T1_6*p + 1), ...
            K0/(m*T*p^2 + 1), K0/(p*(m*T*p + 1)), K0*p/(m*T*p + 1), ...
            K0*(m*T*p + 1), K0*(m*T*p + 1)/p};
w_corner = @(m) {[], [], [], 1/(m*T), [1/(m*Ta), 1/(m*Tb)], 1/(m*T2), ...
                 1/sqrt(m*T), 1/(m*T), 1/(m*T), 1/(m*T), 1/(m*T)};
labelsM = compose('%.1f*T', m_set);

tlTL   = newFigure('Влияние T на ЛАЧХ', 2, 4);
tlTPsi = newFigure('Влияние T на ЛФЧХ', 2, 4);

fprintf('\nСопрягающие частоты, рад/с (звенья 1-3 не содержат T)\n');
fprintf('%-34s', 'Звено'); fprintf('%-18s', labelsM{:}); fprintf('\n');
for i = 4:n
    axL   = nexttile(tlTL);   hold(axL, 'on');   grid(axL, 'on');   box(axL, 'on');
    axPsi = nexttile(tlTPsi); hold(axPsi, 'on'); grid(axPsi, 'on'); box(axPsi, 'on');
    fprintf('%-34s', titles{i});

    for k = 1:numel(m_set)
        Wm = W_m(m_set(k));
        wk = w_corner(m_set(k));
        [mag, ph] = bode(Wm{i}, w);

        semilogx(axL, w, 20*log10(squeeze(mag)), 'Color', cols(k, :), 'LineWidth', 1.3);
        semilogx(axPsi, w, squeeze(ph), 'Color', cols(k, :), 'LineWidth', 1.3);
        xline(axL, wk{i}, '--', 'Color', cols(k, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        xline(axPsi, wk{i}, '--', 'Color', cols(k, :), 'LineWidth', 1, 'HandleVisibility', 'off');

        fprintf('%-18s', listStr(wk{i}));
    end
    fprintf('\n');

    set([axL axPsi], 'XScale', 'log');
    if i == 7
        ylim(axL, [-60 60]);
    end
    title(axL, titles{i}, 'FontSize', 9); title(axPsi, titles{i}, 'FontSize', 9);
    xlabel(axL, '\omega, рад/с'); ylabel(axL, 'L, дБ');
    xlabel(axPsi, '\omega, рад/с'); ylabel(axPsi, '\psi, град');
end
labelsT = compose('T = %.1f\\cdotT_{вар}', m_set);
lg = legend(nexttile(tlTL, 1), labelsT, 'Orientation', 'horizontal');   lg.Layout.Tile = 'south';
lg = legend(nexttile(tlTPsi, 1), labelsT, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
fprintf('Звено 5: ширина участка -20 дБ/дек lg(Ta/Tb) = %.3f дек при любом масштабе T\n', log10(Ta/Tb));

T1_list = T1_6*[0.5 1 2];
tl6   = newFigure(sprintf('Колебательное звено: T_2 = %.1f с, изменение T_1', T2), 2, 1);
axL   = nexttile(tl6); hold(axL, 'on');   grid(axL, 'on');   box(axL, 'on');
axPsi = nexttile(tl6); hold(axPsi, 'on'); grid(axPsi, 'on'); box(axPsi, 'on');

fprintf('\nЗвено 6 (T2 = %.1f с): резонансный пик ЛАЧХ, 20lgK = %.2f дБ\n', T2, 20*log10(K0));
fprintf('   T1, с    zeta    Lmax, дБ   w, рад/с\n');
for k = 1:numel(T1_list)
    W6 = K0/(T2^2*p^2 + T1_list(k)*p + 1);
    [mag, ph] = bode(W6, w);
    L = 20*log10(squeeze(mag));
    [Lmax, idx] = max(L);
    semilogx(axL, w, L, 'Color', cols(k, :), 'LineWidth', 1.3);
    semilogx(axPsi, w, squeeze(ph), 'Color', cols(k, :), 'LineWidth', 1.3);
    fprintf('  %6.2f   %5.3f   %7.2f    %6.3f\n', T1_list(k), T1_list(k)/(2*T2), Lmax, w(idx));
end
set([axL axPsi], 'XScale', 'log');
xline(axL, 1/T2, 'k--', 'HandleVisibility', 'off');
xline(axPsi, 1/T2, 'k--', 'HandleVisibility', 'off');
xlabel(axPsi, '\omega, рад/с'); ylabel(axL, 'L, дБ'); ylabel(axPsi, '\psi, град');
legend(axL, compose('T_1 = %.2f с', T1_list), 'Location', 'southwest');

function tl = newFigure(figName, rows, cols)
figure('Name', figName, 'Color', 'w', 'Position', [40 40 1400 900]);
tl = tiledlayout(rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, figName, 'FontWeight', 'bold');
end

function y = evalOnGrid(expr, x, xs)
f = matlabFunction(expr, 'Vars', x);
y = f(xs);
if isscalar(y)
    y = repmat(y, size(xs));
end
end

function La = asymptote(w, K, a)
[q, Tn, Td, Td2] = a{:};
La = 20*log10(K) + 20*q*log10(w);
for Ti = Tn,  La = La + max(0, 20*log10(w*Ti)); end
for Ti = Td,  La = La - max(0, 20*log10(w*Ti)); end
for Ti = Td2, La = La - max(0, 40*log10(w*Ti)); end
end

function wc = crossings(w, L)
lw  = log10(w(:));
L   = L(:);
idx = find(L(1:end-1).*L(2:end) < 0);
wc  = 10.^(lw(idx) - L(idx).*(lw(idx+1) - lw(idx))./(L(idx+1) - L(idx)));
wc  = wc';
end

function str = listStr(x)
if isempty(x)
    str = '-';
else
    str = strjoin(compose('%.3f', x(:)'), ', ');
end
end

function plotAFC(ax, Wi, w, K)
[re, im] = nyquist(Wi, w);
re = squeeze(re); im = squeeze(im);
R = 3*K;
far = abs(re) > 1.1*R | abs(im) > 1.1*R;
re(far) = NaN; im(far) = NaN;

plot(ax, re, im, 'b', 'LineWidth', 1.3);
if ~far(1)
    plot(ax, re(1), im(1), 'bo', 'MarkerFaceColor', 'b');
end
xline(ax, 0, 'k:'); yline(ax, 0, 'k:');
if any(far)
    axis(ax, [-R R -R R]);
else
    xr = [min([re; 0]), max([re; 0])];
    yr = [min([im; 0]), max([im; 0])];
    half = 0.6*max(diff(xr), diff(yr));
    axis(ax, [mean(xr) + [-half half], mean(yr) + [-half half]]);
end
axis(ax, 'square');
xlabel(ax, 'Re W(j\omega)'); ylabel(ax, 'Im W(j\omega)');
end

function yl = limitsL(L)
lo = max(min(L), -80);
hi = min(max(L), 60);
if hi - lo < 10
    lo = lo - 10; hi = hi + 10;
end
yl = [lo - 5, hi + 5];
end
