clc;
clear;
close all;

s = tf('s');

K = 3.2;
T = 0.9;
T1 = 1.0;
T2 = 1.4;
T1_ap = T2 * 2;

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

groups = {1:6, 7:11};         
ncol   = 3;
w_bode = logspace(-2, 2, 1000);

%% АФЧХ
for g = 1:numel(groups)
    idx  = groups{g};
    nrow = ceil(numel(idx) / ncol);

    figure('Name', sprintf('АФЧХ, окно %d', g), ...
           'NumberTitle', 'off', 'WindowState', 'maximized');

    for k = 1:numel(idx)
        i = idx(k);
        subplot(nrow, ncol, k);

        w_range = logspace(-2, 2, 500);
        if i == 7
            w_range = logspace(-1, 0.5, 500);   % консервативное звено: суженный диапазон
        end

        nyquist(W{i}, w_range);
        grid on;
        title(names{i}, 'FontSize', 9);

        ax = gca;
        if ismember(i, [1, 2, 7, 8, 10, 11])
            ax.XLim = [-5 5];
            ax.YLim = [-5 5];
        elseif i == 6
            ax.XLim = [-6 6];
            ax.YLim = [-8 8];
        end
    end
    put_sgtitle(gcf, sprintf('АФЧХ звеньев, окно %d', g));
end

%% ЛАЧХ и ЛФЧХ
for g = 1:numel(groups)
    idx  = groups{g};
    n    = numel(idx);
    nrow = ceil(n / ncol);

    figure('Name', sprintf('ЛЧХ, окно %d', g), ...
           'NumberTitle', 'off', 'WindowState', 'maximized');

    for k = 1:n
        i = idx(k);
        r = floor((k-1) / ncol);
        c = mod(k-1, ncol);

        axM = subplot(2*nrow, ncol, 2*r*ncol + c + 1);   % сверху ЛАЧХ
        axP = subplot(2*nrow, ncol, (2*r+1)*ncol + c + 1);   % снизу ЛФЧХ

        [mag, phase] = bode(W{i}, w_bode);
        semilogx(axM, w_bode, 20*log10(squeeze(mag)), 'LineWidth', 1.2);
        semilogx(axP, w_bode, squeeze(phase),         'LineWidth', 1.2);

        grid(axM, 'on');  grid(axP, 'on');
        axM.XLim = [1e-2 1e2];   axP.XLim = [1e-2 1e2];
        ylim(axM, [-80 80]);
        ylim(axP, [-200 120]);

        title(axM, names{i}, 'FontSize', 9);
        ylabel(axM, 'L, дБ');
        ylabel(axP, '\phi, град');
        xlabel(axP, '\omega, рад/с');
    end
    put_sgtitle(gcf, sprintf('ЛАЧХ и ЛФЧХ звеньев, окно %d', g));
end

%% Определяем влияние коэффициентов K и T
K_set = K * [0.5 1 2];    
T_set = T * [0.5 1 2];     
zeta  = T1 / (2*T2);      

% Параметрические модели звеньев W(K, Tv)
Wp = {
    @(Kv,Tv) tf(Kv/s);                              
    @(Kv,Tv) tf(Kv*s);                      
    @(Kv,Tv) tf(Kv);                             
    @(Kv,Tv) tf(Kv/(Tv*s + 1));                    
    @(Kv,Tv) tf(Kv/(Tv^2*s^2 + 2*Tv*s + 1));        
    @(Kv,Tv) tf(Kv/(Tv^2*s^2 + 2*zeta*Tv*s + 1));   
    @(Kv,Tv) tf(Kv/(Tv*s^2 + 1));                   
    @(Kv,Tv) tf(Kv/(s*(Tv*s + 1)));                  
    @(Kv,Tv) tf(Kv*s/(Tv*s + 1));                   
    @(Kv,Tv) tf(Kv*(Tv*s + 1));                      
    @(Kv,Tv) tf(Kv*(Tv*s + 1)/s)                    
};

w    = logspace(-3, 3, 20000);
cols = lines(3);                                      
lblK = arrayfun(@(k) sprintf('K = %.2f', k), K_set, 'UniformOutput', false);
lblT = arrayfun(@(t) sprintf('T = %.2f', t), T_set, 'UniformOutput', false);

%% Графики: влияние K и T на ЛАЧХ и ЛФЧХ
for mode = 1:2
    if mode == 1
        groupsKT = {1:6, 7:11};        % влияние K
        lbl = lblK;  par = 'K';
    else
        groupsKT = {4:7, 8:11};        % влияние T
        lbl = lblT;  par = 'T';
    end

    for g = 1:numel(groupsKT)
        idx  = groupsKT{g};
        n    = numel(idx);
        ncolKT = 3;  if n <= 4, ncolKT = 2; end
        nrow = ceil(n / ncolKT);

        figure('Name', sprintf('Влияние %s, окно %d', par, g), ...
               'NumberTitle', 'off', 'WindowState', 'maximized');

        for k = 1:n
            i = idx(k);
            r = floor((k-1) / ncolKT);
            c = mod(k-1, ncolKT);
            axM = subplot(2*nrow, ncolKT, 2*r*ncolKT + c + 1);   % сверху ЛАЧХ
            axP = subplot(2*nrow, ncolKT, (2*r+1)*ncolKT + c + 1);   % снизу ЛФЧХ
            hold(axM, 'on');  hold(axP, 'on');

            for j = 1:3
                if mode == 1
                    sys = Wp{i}(K_set(j), T);
                else
                    sys = Wp{i}(K, T_set(j));
                end
                [m, ph] = bode(sys, w);
                plot(axM, w, 20*log10(squeeze(m)), 'Color', cols(j,:), 'LineWidth', 1.2);
                plot(axP, w, squeeze(ph),          'Color', cols(j,:), 'LineWidth', 1.2);
            end

            set([axM axP], 'XScale', 'log', 'XLim', [1e-3 1e3]);
            grid(axM, 'on');  grid(axP, 'on');
            ylim(axM, [-80 80]);
            ylim(axP, [-200 120]);
            title(axM, names{i}, 'FontSize', 9);
            ylabel(axM, 'L, дБ');
            ylabel(axP, '\phi, град');
            xlabel(axP, '\omega, рад/с');
            if k == 1
                legend(axM, lbl, 'Location', 'best', 'FontSize', 7);
            end
        end
        put_sgtitle(gcf, sprintf('Влияние %s на ЛАЧХ и ЛФЧХ', par));
    end
end

function put_sgtitle(fig, txt)
    axs = findobj(fig, 'Type', 'axes');
    for a = axs'
        p = a.Position;
        a.Position = [p(1), p(2)*0.93, p(3), p(4)*0.93];
    end
    annotation(fig, 'textbox', [0 0.95 1 0.05], 'String', txt, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');
end
