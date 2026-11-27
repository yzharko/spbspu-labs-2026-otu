clear; clc; close all;

K0 = 2.4;
T0 = 0.4;
T1 = 0.8;
T2 = 0.9;

K_values = [0.5 1 2 3 5];
T_values = [0.1 0.5 1 2 5];
w = logspace(-2, 3, 2000);

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
         'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
         'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
         'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};
kinds = {'АФЧХ', 'ЛАЧХ', 'ЛФЧХ'};

for k = 1:3
    figure('Name', [kinds{k} ': все звенья'], 'NumberTitle', 'off', ...
        'Position', [50 50 1400 850]);
    for n = 1:11
        subplot(3, 4, n);
        W = {make_W(n, K0, T0, T1, T2)};
        draw(kinds{k}, n, W, {''}, T0, w);
        title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
    end
    sgtitle([kinds{k} ' всех типовых звеньев']);
end

for k = 2:3
    figure('Name', [kinds{k} ': влияние K'], 'NumberTitle', 'off', ...
        'Position', [50 50 1400 850]);
    for n = 1:11
        subplot(3, 4, n);
        W = cell(1, numel(K_values));
        labs = cell(1, numel(K_values));
        for j = 1:numel(K_values)
            W{j} = make_W(n, K_values(j), T0, T1, T2);
            labs{j} = sprintf('K = %.2g', K_values(j));
        end
        draw(kinds{k}, n, W, labs, T0*ones(1, numel(W)), w);
        title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
    end
    sgtitle([kinds{k} ': влияние коэффициента K']);
end

tile_n   = [4   5     5     6     6     7  8  9  10 11];
tile_par = {'T','T1','T2','T1','T2','T','T','T','T','T'};
tile_val = {T_values, [2 2.5 3], [0.1 0.5 0.9], ...
            [0.2 0.5 0.8 1.2 1.6], [0.5 0.9 1.5 2], ...
            T_values, T_values, T_values, T_values, T_values};

for k = 2:3
    figure('Name', [kinds{k} ': влияние T'], 'NumberTitle', 'off', ...
        'Position', [50 50 1500 700]);
    for q = 1:numel(tile_n)
        subplot(2, 5, q);
        n = tile_n(q);
        par = tile_par{q};
        vals = tile_val{q};
        W = cell(1, numel(vals));
        labs = cell(1, numel(vals));
        Tres = T0*ones(1, numel(vals));
        for j = 1:numel(vals)
            T = T0;  a = T1;  b = T2;
            if strcmp(par, 'T'),  T = vals(j); Tres(j) = T; end
            if strcmp(par, 'T1'), a = vals(j); end
            if strcmp(par, 'T2'), b = vals(j); end
            W{j} = make_W(n, K0, T, a, b);
            labs{j} = sprintf('%s = %.2g', par, vals(j));
        end
        draw(kinds{k}, n, W, labs, Tres, w);
        title(sprintf('%d. %s: %s', n, names{n}, par), 'FontSize', 8);
    end
    sgtitle([kinds{k} ': влияние постоянных времени']);
end

function W = make_W(n, K, T, T1, T2)
    s = tf('s');
    switch n
        case 1,  W = K/s;
        case 2,  W = K*s;
        case 3,  W = tf(K);
        case 4,  W = K/(T*s + 1);
        case 5,  W = K/(T2^2*s^2 + T1*s + 1);
        case 6,  W = K/(T2^2*s^2 + T1*s + 1);
        case 7,  W = K/(T*s^2 + 1);
        case 8,  W = K/(s*(T*s + 1));
        case 9,  W = K*s/(T*s + 1);
        case 10, W = K*(T*s + 1);
        case 11, W = K*(T*s + 1)/s;
    end
end

function draw(kind, n, W, labs, Tres, w)
    colors = lines(numel(W));
    h = zeros(1, numel(W));
    rmax = 0;
    hold on;
    grid on;

    for j = 1:numel(W)
        [mag, ph] = bode(W{j}, w);
        [re, im] = nyquist(W{j}, w);
        mag = squeeze(mag);
        ph  = squeeze(ph);
        re  = squeeze(re);
        im  = squeeze(im);

        if n == 7
            im(:) = 0;
            [~, idx] = min(abs(w - 1/sqrt(Tres(j))));
            mag(idx) = NaN;
            ph(idx)  = NaN;
            re(idx)  = NaN;
        end
        rmax = max(rmax, max(abs(re)));

        if strcmp(kind, 'ЛАЧХ')
            h(j) = plot(w, 20*log10(mag), 'Color', colors(j, :));
        elseif strcmp(kind, 'ЛФЧХ')
            h(j) = plot(w, ph, 'Color', colors(j, :));
        else
            if n == 3
                h(j) = plot(re(1), im(1), 'o', 'Color', colors(j, :), ...
                    'MarkerFaceColor', colors(j, :));
            else
                h(j) = plot(re, im, 'Color', colors(j, :));
                plot(re, -im, ':', 'Color', colors(j, :));
            end
        end
    end

    if strcmp(kind, 'ЛАЧХ')
        set(gca, 'XScale', 'log');
        xlim([w(1) w(end)]);
        xlabel('\omega, рад/с');
        ylabel('L(\omega), дБ');
    elseif strcmp(kind, 'ЛФЧХ')
        set(gca, 'XScale', 'log');
        xlim([w(1) w(end)]);
        xlabel('\omega, рад/с');
        ylabel('\psi(\omega), град');
    else
        xlabel('Re W(j\omega)');
        ylabel('Im W(j\omega)');
        if n == 1 || n == 2
            xlim([-1 1]);
        elseif n == 7
            xlim([-10 10]);
            ylim([-1 1]);
        elseif n == 10 || n == 11
            xlim([0 2*rmax]);
        end
    end

    if numel(W) > 1
        legend(h, labs, 'Location', 'best', 'FontSize', 6);
    end
end
