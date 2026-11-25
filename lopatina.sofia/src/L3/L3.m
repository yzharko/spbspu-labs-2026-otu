clc; clear; close all;

s = tf('s');

K0 = 0.8;
T0 = 0.5;
T1 = 0.7;
T2 = 0.5;
T1a = 2*T2;
K_values = [0.8 1 2 5];
T_values = [0.5 1 2 5];

w = logspace(-2, 2, 1000);

names = {
    'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', 'Реальное дифференцирующее', ...
    'Форсирующее', 'Изодромное'
    };

W = {
    @(K,T) K/s
    @(K,T) K*s
    @(K,T) tf(K)
    @(K,T) K/(T*s + 1)
    @(K,T) K/(T2^2*s^2 + T1a*s + 1)
    @(K,T) K/(T2^2*s^2 + T1*s + 1)
    @(K,T) K/(T*s^2 + 1)
    @(K,T) K/(s*(T*s + 1))
    @(K,T) K*s/(T*s + 1)
    @(K,T) K*(T*s + 1)
    @(K,T) K*(T*s + 1)/s
};

% Группировка звеньев для графиков
groups = {1:4, 5:8, 9:11};
groupsT = {[4 7 8 9], [10 11]};

%% АФЧХ
for g = 1:length(groups)
    figure('Color', 'w', 'Position', [100 80 1100 800]);
    for k = 1:length(groups{g})
        n = groups{g}(k);
        subplot(2, 2, k);
        nyquist(W{n}(K0, T0), w);
        grid on;
        title(names{n});
    end
end

%% Вспомогательная функция для построения графиков ЛАЧХ и ЛФЧХ
function bode_window(idxs, W, names, K, T, labels, w)
figure('Color', 'w', 'Position', [100 80 1100 800]);
outer = tiledlayout(2, 2, 'TileSpacing', 'loose');
for k = 1:length(idxs)
    n = idxs(k);
    inner = tiledlayout(outer, 2, 1, 'TileSpacing', 'compact');
    inner.Layout.Tile = k;
    ax1 = nexttile(inner);  hold(ax1, 'on');  grid(ax1, 'on');
    ax2 = nexttile(inner);  hold(ax2, 'on');  grid(ax2, 'on');
    for j = 1:length(K)
        [mag, ph] = bode(W{n}(K(j), T(j)), w);
        semilogx(ax1, w, 20*log10(squeeze(mag)));
        semilogx(ax2, w, squeeze(ph));
    end
    set([ax1 ax2], 'XScale', 'log');
    xlim(ax1, [w(1) w(end)]);
    xlim(ax2, [w(1) w(end)]);
    title(ax1, names{n});
    ylabel(ax1, 'L, дБ');
    ylabel(ax2, '\psi, град');
    xlabel(ax2, '\omega, рад/с');
    if k == 1
        legend(ax1, labels, 'Location', 'best');
    end
end
end

%% ЛАЧХ и ЛФЧХ при изменении K
labelsK = arrayfun(@(k) sprintf('K = %.2g', k), K_values, 'UniformOutput', false);
for g = 1:length(groups)
    bode_window(groups{g}, W, names, K_values, T0*ones(size(K_values)), labelsK, w);
end

%% ЛАЧХ и ЛФЧХ при изменении T
labelsT = arrayfun(@(t) sprintf('T = %.2g', t), T_values, 'UniformOutput', false);
for g = 1:length(groupsT)
    bode_window(groupsT{g}, W, names, K0*ones(size(T_values)), T_values, labelsT, w);
end

%% Сопрягающие частоты
fprintf('Сопрягающие частоты w_c = 1/T:\n');
for i = 1:length(T_values)
    fprintf('T = %.2f c, omega_c = %.3f рад/с\n', T_values(i), 1/T_values(i));
end
