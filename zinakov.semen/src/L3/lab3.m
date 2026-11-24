clc; clear; close all;

s = tf('s');

K   = 2.2;
T   = 0.2;
T1  = 0.4;
T2  = 0.3;
% для апериодического звена 2-го порядка, тк там надо T1 >= 2*T2 = 0.6 с
T1a = 0.7;

names = {'интегрирующего звена', ...
         'дифференцирующего звена', ...
         'усилительного звена', ...
         'апериодического звена 1-го порядка', ...
         'апериодического звена 2-го порядка', ...
         'колебательного звена', ...
         'консервативного звена', ...
         'инерционного интегрирующего звена', ...
         'реального дифференцирующего звена', ...
         'форсирующего звена 1-го порядка', ...
         'изодромного звена'};

W = {K/s, ...
     K*s, ...
     tf(K), ...
     K/(T*s + 1), ...
     K/(T2^2*s^2 + T1a*s + 1), ...
     K/(T2^2*s^2 + T1*s + 1), ...
     K/(T*s^2 + 1), ...
     K/(s*(T*s + 1)), ...
     K*s/(T*s + 1), ...
     K*(T*s + 1), ...
     K*(T*s + 1)/s};

links_num = numel(W);

for i = 1:links_num
    figure;
    nyquist(W{i});
    grid on;
    title(['АФЧХ ' names{i}]);

    figure;
    bode(W{i});
    grid on;
    title(['ЛАЧХ и ЛФЧХ ' names{i}]);
end

omega_dec_min   = -2;
omega_dec_max   = 2;
num_freq_points = 2000;
w = logspace(omega_dec_min, omega_dec_max, num_freq_points);

K_values = [0.5 1 2.2 3.5 5];
famK = struct('tfs', {}, 'labels', {}, 'title', {});
for n = 1:links_num
    famK(n).tfs    = cellfun(@(k) W{n}*k/K, num2cell(K_values), 'UniformOutput', false);
    famK(n).labels = compose('K = %.2g', K_values);
    famK(n).title  = ['№' num2str(n) ' ' names{n}];
end

draw_canvas(famK, w, 'L',   'Влияние K на ЛАЧХ типовых звеньев');
draw_canvas(famK, w, 'phi', 'Влияние K на ЛФЧХ типовых звеньев');

T_values = [0.1 0.2 0.5 1 2];
T1_5     = [0.6 0.7 1.0 1.5 2.5];
T2_5     = [0.1 0.18 0.25 0.3 0.35];
T1_6     = [0.05 0.15 0.25 0.4 0.55];
T2_6     = [0.22 0.3 0.5 0.8 1.2];

panels = { ...
    '4 (T)',   T_values, 'T',  @(p) K/(p*s + 1); ...
    '5 (T1)',  T1_5,     'T1', @(p) K/(T2^2*s^2 + p*s + 1); ...
    '5 (T2)',  T2_5,     'T2', @(p) K/(p^2*s^2 + T1a*s + 1); ...
    '6 (T1)',  T1_6,     'T1', @(p) K/(T2^2*s^2 + p*s + 1); ...
    '6 (T2)',  T2_6,     'T2', @(p) K/(p^2*s^2 + T1*s + 1); ...
    '7 (T)',   T_values, 'T',  @(p) K/(p*s^2 + 1); ...
    '8 (T)',   T_values, 'T',  @(p) K/(s*(p*s + 1)); ...
    '9 (T)',   T_values, 'T',  @(p) K*s/(p*s + 1); ...
    '10 (T)',  T_values, 'T',  @(p) K*(p*s + 1); ...
    '11 (T)',  T_values, 'T',  @(p) K*(p*s + 1)/s};

panels_num = size(panels, 1);
famT = struct('tfs', {}, 'labels', {}, 'title', {});
for n = 1:panels_num
    famT(n).tfs    = cellfun(panels{n, 4}, num2cell(panels{n, 2}), 'UniformOutput', false);
    famT(n).labels = compose([panels{n, 3} ' = %.2g'], panels{n, 2});
    famT(n).title  = panels{n, 1};
end

draw_canvas(famT, w, 'L',   'Влияние постоянных времени на ЛАЧХ (K = 2.2)');
draw_canvas(famT, w, 'phi', 'Влияние постоянных времени на ЛФЧХ (K = 2.2)');

function draw_canvas(fam, w, mode, name)
    figure('Name', name, 'NumberTitle', 'off');
    n_rows    = 3;
    n_cols    = 4;
    dB_factor = 20;
    font_sz   = 7;

    for n = 1:numel(fam)
        subplot(n_rows, n_cols, n);
        hold on;
        for i = 1:numel(fam(n).tfs)
            [mag, phase] = bode(fam(n).tfs{i}, w);
            if strcmp(mode, 'L')
                semilogx(w, dB_factor*log10(squeeze(mag)));
            else
                semilogx(w, squeeze(phase));
            end
        end
        set(gca, 'XScale', 'log');
        grid on;
        title(fam(n).title);
        xlabel('\omega, рад/с');
        if strcmp(mode, 'L')
            ylabel('L, дБ');
        else
            ylabel('\psi, °');
        end
        legend(fam(n).labels, 'Location', 'best', 'FontSize', font_sz);
    end
    sgtitle(name);
end
