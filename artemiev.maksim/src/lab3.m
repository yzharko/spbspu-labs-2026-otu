clear; close all; clc;

K = 3.2;
T = 0.9;
T1 = 1.0;
T2 = 1.4;
T1_ap = T2 * 2;

s = tf('s');

names = {
    'Интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1го порядка';
    'Апериодическое 2го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1го порядка';
    'Изодромное'};

W = {
    tf(K/s);
    tf(K*s);
    tf(K);
    tf(K/(T*s + 1));
    tf(K/(T2^2*s^2 + T1_ap*s + 1));
    tf(K/(T2^2*s^2 + T1*s + 1));
    tf(K/(T*s^2 + 1));
    tf(K/(s*(T*s + 1)));
    tf(K*s/(T*s + 1));
    tf(K*(T*s + 1));
    tf(K*(T*s + 1)/s)
};

for i = 1:length(W)
    figure('Name', ['АФЧХ: ', names{i}], 'NumberTitle', 'off');
    w_range = logspace(-2, 2, 500); 

    if i == 7
        w_range = logspace(-1, 0.5, 500);
    end
    
    nyquist(W{i}, w_range);
    grid on;
    title(['АФЧХ: ', names{i}]);
    
    ax = gca;
    if ismember(i, [1, 2, 7, 8, 10, 11])
        ax.XLim = [-5 5];
        ax.YLim = [-5 5];
    elseif i == 6
        ax.XLim = [-6 6];
        ax.YLim = [-8 8];
    end

    figure('Name', ['ЛЧХ: ', names{i}], 'NumberTitle', 'off');
    bode(W{i});
    grid on;
    title(['ЛАЧХ и ЛФЧХ: ', names{i}]);
end

figure('Name','Сводный график всех АФЧХ', 'NumberTitle', 'off');
for i = 1:length(W)
    subplot(4,3,i);
    w_range = logspace(-2, 2, 500);
    if i == 7
        w_range = logspace(-1, 0.5, 500);
    end
    
    nyquist(W{i}, w_range);
    grid on;
    title(names{i});

    ax = gca;
    if ismember(i, [1, 2, 7, 8, 10, 11])
        ax.XLim = [-3 3];
        ax.YLim = [-3 3];
    end
end

figure('Name','Сводный график всех ЛЧХ', 'NumberTitle', 'off');
for i = 1:length(W)
    % амплитуда
    subplot(11,2,2*i-1);
    bodemag(W{i});
    grid on;
    title(['ЛАЧХ: ', names{i}]);
    xlabel('');
    
    % фаза
    subplot(11,2,2*i);
    [mag, phase, wout] = bode(W{i});
    semilogx(wout, squeeze(phase), 'LineWidth', 1.2);
    grid on;
    title(['ЛФЧХ: ', names{i}]);
    xlabel('\omega, рад/с');
    ylabel('\phi, град');
    ylim([-180 90]);
end