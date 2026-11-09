clc;
clear;
close all;
s = tf('s');

K0 = 1.2;
T0 = 0.3;
T1 = 0.5;
T2 = 0.6;
T15 = 1.5; % Для звена 5: T15 > 2*T2.
K_values = [0.5 1 1.2 2 3 5];
T_values = [0.1 0.3 0.5 1 2 5];
w = logspace(-3, 3, 1500); 

names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

W_all = {K0/s, K0*s, tf(K0, 1), K0/(T0*s+1), ...
    K0/(T2^2*s^2+T15*s+1), K0/(T2^2*s^2+T1*s+1), ...
    K0/(T0*s^2+1), K0/(s*(T0*s+1)), K0*s/(T0*s+1), ...
    K0*(T0*s+1), K0*(T0*s+1)/s};

plot_names = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ'};
group_names = {'Исходные параметры', 'Влияние K', ...
    'Влияние T (для звеньев 5 и 6 — T1)', 'Влияние T2 (звенья 5 и 6)'};
for groupIndex = 1:4
    for plotIndex = 1:3
        figure(3*(groupIndex-1)+plotIndex);
        set(gcf, 'Name', [plot_names{plotIndex} ': ' group_names{groupIndex}], ...
            'NumberTitle', 'off', 'Position', [50 50 1400 850]);
    end
end

for n = 1:numel(W_all) 
    W = W_all{n};
    name = sprintf('%d. %s', n, names{n});

    WK = cell(1, numel(K_values));
    labelsK = cell(1, numel(K_values));

    for i = 1:numel(K_values)
        K = K_values(i);
        WK{i} = (K/K0)*W;
        labelsK{i} = sprintf('K = %.2g', K);
    end

    W_groups = {{W}, WK, {W}};
    labels_groups = {{sprintf('K = %g', K0)}, labelsK, {'T отсутствует'}};

    if n == 4 || n >= 7
        WT = cell(1, numel(T_values));
        labelsT = cell(1, numel(T_values));

        for i = 1:numel(T_values)
            T = T_values(i);
            switch n
                case 4
                    WT{i} = K0/(T*s+1);
                case 7
                    WT{i} = K0/(T*s^2+1);
                case 8
                    WT{i} = K0/(s*(T*s+1));
                case 9
                    WT{i} = K0*s/(T*s+1);
                case 10
                    WT{i} = K0*(T*s+1);
                case 11
                    WT{i} = K0*(T*s+1)/s;
            end
            labelsT{i} = sprintf('T = %.2g', T);
        end
        W_groups{3} = WT;
        labels_groups{3} = labelsT;
    end

    if n == 5 || n == 6
        T1_fixed = T1;
        T1_values = T1*[0.2 0.5 1 2 5];
        T2_values = T2*[0.2 0.5 1 2 5];

        if n == 5
            T1_fixed = T15;
            T1_values = [1.3 1.5 2 3 5];
            T2_values = [0.1 0.3 0.5 0.6 0.7];
        end

        WT1 = cell(1, numel(T1_values));
        WT2 = cell(1, numel(T2_values));

        labelsT1 = cell(1, numel(T1_values));
        labelsT2 = cell(1, numel(T2_values));

        for i = 1:numel(T1_values)
            WT1{i} = K0/(T2^2*s^2+T1_values(i)*s+1);
            labelsT1{i} = sprintf('T1 = %.2g, ζ = %.2g', ...
                T1_values(i), T1_values(i)/(2*T2));
        end

        for i = 1:numel(T2_values)
            WT2{i} = K0/(T2_values(i)^2*s^2+T1_fixed*s+1);
            labelsT2{i} = sprintf('T2 = %.2g, ζ = %.2g', ...
                T2_values(i), T1_fixed/(2*T2_values(i)));
        end

        W_groups{3} = WT1;
        W_groups{4} = WT2;
        labels_groups{3} = labelsT1;
        labels_groups{4} = labelsT2;
    end


    for groupIndex = 1:numel(W_groups)
        layout = [3 4 n];
        if groupIndex == 4 
            layout = [1 2 n-4];
        end

        models = W_groups{groupIndex};
        colors = lines(numel(models));
        
        for i = 1:numel(models)
            [mag, phase] = bode(models{i}, w);
            [re, im] = nyquist(models{i}, w);
            mag = mag(:);
            phase = phase(:);
            re = re(:);
            im = im(:);
            if n == 7 
                gap = find(diff(sign(re)) ~= 0) + 1;
                re(gap) = NaN;
                mag(gap) = NaN;
                phase(gap) = NaN;
            end

            figure(3*(groupIndex-1)+1);
            subplot(layout(1), layout(2), layout(3));
            semilogx(w, 20*log10(mag), 'Color', colors(i,:));
            hold on;
            grid on;
            title(name, 'FontSize', 8);
            xlabel('Частота, рад/с');
            ylabel('L, дБ');

            figure(3*(groupIndex-1)+2);
            subplot(layout(1), layout(2), layout(3));
            semilogx(w, phase, 'Color', colors(i,:));
            hold on;
            grid on;
            title(name, 'FontSize', 8);
            xlabel('Частота, рад/с');
            ylabel('Фаза, град');

            figure(3*(groupIndex-1)+3);
            subplot(layout(1), layout(2), layout(3));
            if n == 3 
                plot(re, im, 'o', 'Color', colors(i,:));
            else
                plot([re; NaN; re], [im; NaN; -im], 'Color', colors(i,:));
            end
            hold on;
            grid on;
            title(name, 'FontSize', 8);
            xlabel('Re W');
            ylabel('Im W');
        end
        
        for p = 1:3
            figure(3*(groupIndex-1)+p);
            subplot(layout(1), layout(2), layout(3));
            legend(labels_groups{groupIndex}, 'Location', 'best', 'FontSize', 6);
        end
    end
end
