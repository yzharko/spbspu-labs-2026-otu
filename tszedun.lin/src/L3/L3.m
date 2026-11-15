clc;
clear;
close all;

K = 4;
T = 1.5;
T1 = 1.6;
T2 = 2;
T1_5 = 5; % Для звена 5: T1 >= 2*T2.

s = tf('s');
make_models = @(k,t,t1,t2,t15) {
    k/s;
    k*s;
    tf(k);
    k/(t*s + 1);
    k/(t2^2*s^2 + t15*s + 1);
    k/(t2^2*s^2 + t1*s + 1);
    k/(t*s^2 + 1);
    k/(s*(t*s + 1));
    k*s/(t*s + 1);
    k*(t*s + 1);
    k*(t*s + 1)/s
};

names = {
    'Интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1-го порядка';
    'Апериодическое 2-го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1-го порядка';
    'Изодромное'
};

omega = logspace(-3, 3, 2000);
W = make_models(K,T,T1,T2,T1_5);
outDir = fullfile(fileparts(mfilename('fullpath')), 'results_L3');
if ~isfolder(outDir)
    mkdir(outDir);
end

%% АФЧХ, ЛАЧХ и ЛФЧХ всех 11 звеньев
for n = 1:numel(W)
    draw_plots(W(n), omega, n, sprintf('%02d. %s', n,names{n}), {}, ...
        fullfile(outDir, sprintf('base_%02d.png', n)));
end

%% Влияние K при неизменных постоянных времени
K_values = [2 4 8];
labelsK = {'K = 2', 'K = 4', 'K = 8'};
for n = 1:numel(W)
    WK = cell(1, numel(K_values));
    for j = 1:numel(K_values)
        WK{j} = (K_values(j)/K)*W{n};
    end
    draw_plots(WK, omega, n, sprintf('%02d. %s: влияние K', n,names{n}), labelsK, ...
        fullfile(outDir, sprintf('vary_K_%02d.png', n)));
end

%% Влияние постоянных времени при K = 4
factors = [0.5 1 2];
WT = cell(numel(W), numel(factors));
for j = 1:numel(factors)
    a = factors(j);
    WT(:,j) = make_models(K,a*T,a*T1,a*T2,a*T1_5);
end

% У звеньев 1-3 постоянных времени нет.
for n = 4:numel(W)
    labelsT = cell(1, numel(factors));
    for j = 1:numel(factors)
        a = factors(j);
        if n == 5
            labelsT{j} = sprintf('T1 = %g, T2 = %g', a*T1_5,a*T2);
        elseif n == 6
            labelsT{j} = sprintf('T1 = %g, T2 = %g', a*T1,a*T2);
        else
            labelsT{j} = sprintf('T = %g', a*T);
        end
    end
    draw_plots(WT(n,:), omega, n, ...
        sprintf('%02d. %s: влияние постоянных времени', n,names{n}), labelsT, ...
        fullfile(outDir, sprintf('vary_T_%02d.png', n)));
end

fprintf('Графики сохранены в папке:\n%s\n', outDir);

function draw_plots(models, omega, number, caption, labels, filename)
    base = numel(models) == 1;
    columns = 2 + base;
    fig = figure('Visible','off', 'Color','white', ...
        'Position',[50 50 450*columns 470]);
    colors = [0 0.447 0.741; 0.85 0.325 0.098; 0.15 0.6 0.25];
    styles = {'-', '--', ':'};

    for j = 1:numel(models)
        [mag, phase] = bode(models{j}, omega);
        mag = mag(:);
        phase = phase(:);
        gap = [];
        if number == 7
            [~, den] = tfdata(models{j}, 'v');
            resonance = sqrt(den(end)/den(1));
            gap = find(omega >= resonance, 1);
            % Не соединяем точки по разные стороны разрыва.
            mag(gap) = NaN;
            phase(gap) = NaN;
        end

        subplot(1,columns,1 + base);
        semilogx(omega,20*log10(mag),styles{j}, 'Color',colors(j,:), 'LineWidth',1.5);
        hold on;
        subplot(1,columns,2 + base);
        semilogx(omega,phase,styles{j}, 'Color',colors(j,:), 'LineWidth',1.5);
        hold on;

        if base
            subplot(1,columns,1);
            [re, im] = nyquist(models{j}, omega);
            re = re(:);
            im = im(:);
            re(gap) = NaN;
            im(gap) = NaN;
            if number == 3
                plot(re(1),im(1),'o', 'Color',colors(1,:), 'MarkerFaceColor',colors(1,:));
            else
                plot(re,im, 'Color',colors(1,:), 'LineWidth',1.5);
            end
            hold on;
            if ismember(number,[4 5 6 9])
                axis equal;
            end
            % Видимая область для ветвей, уходящих в бесконечность.
            if ismember(number,[1 8 11])
                ylim([-16 0]);
            elseif ismember(number,[2 10])
                ylim([0 16]);
            elseif number == 7
                xlim([-8 8]);
                ylim([-1 1]);
            end
            if ismember(number,[1 2])
                xlim([-1 1]);
            elseif ismember(number,[3 10])
                xlim([0 8]);
            elseif number == 11
                xlim([0 12]);
            end
            if number ~= 3
                arrow_frequencies = 1;
                if number == 7
                    arrow_frequencies = [0.3 1.2];
                end
                for frequency = arrow_frequencies
                    [~, q] = min(abs(omega - frequency));
                    d = [re(q+1)-re(q), im(q+1)-im(q)];
                    span = [diff(xlim), diff(ylim)];
                    d = 0.08*d/norm(d./span);
                    quiver(re(q),im(q),d(1),d(2),0, 'Color',colors(1,:), ...
                        'LineWidth',1.5, 'MaxHeadSize',1);
                end
            end
            title('АФЧХ, \omega > 0');
            xlabel('Re W(j\omega), отн. ед.');
            ylabel('Im W(j\omega), отн. ед.');
        end
    end

    for p = 1:2
        subplot(1,columns,p + base);
        xlim([omega(1) omega(end)]);
        xlabel('\omega, рад/с');
        if p == 1
            title('ЛАЧХ');
            ylabel('L(\omega), дБ');
        else
            title('ЛФЧХ');
            ylabel('\phi(\omega), град');
        end
        if ~isempty(labels)
            legend(labels,'Location','best', 'Color','white','TextColor','black');
        end
    end
    for p = 1:columns
        ax = subplot(1,columns,p);
        grid(ax,'on');
        box(ax,'on');
        set(ax,'Color','white','XColor','black','YColor','black','FontSize',11);
        set([ax.Title ax.XLabel ax.YLabel],'Color','black');
    end
    sgtitle(caption,'Color','black','FontSize',14);
    exportgraphics(fig,filename,'Resolution',200,'BackgroundColor','white');
    close(fig);
end
