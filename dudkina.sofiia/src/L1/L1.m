K = 1;
T = 0.2;
T1 = 0.4;
T2 = 0.5;
T1_ap2 = 1;

s = tf('s');

links(1).name = 'Интегрирующее';
links(1).W = K/s;

links(2).name = 'Дифференцирующее';
links(2).W = K*s;

links(3).name = 'Усилительное';
links(3).W = tf(K, 1);

links(4).name = 'Апериодическое 1-го порядка';
links(4).W = K/(T*s + 1);

links(5).name = 'Апериодическое 2-го порядка';
links(5).W = K/(T2^2*s^2 + T1_ap2*s + 1);

links(6).name = 'Колебательное';
links(6).W = K/(T2^2*s^2 + T1*s + 1);

links(7).name = 'Консервативное';
links(7).W = K/(T*s^2 + 1);

links(8).name = 'Инерционное интегрирующее';
links(8).W = K/(T*s^2 + s);

links(9).name = 'Реальное дифференцирующее';
links(9).W = K*s/(T*s + 1);

links(10).name = 'Форсирующее 1-го порядка';
links(10).W = K*(T*s + 1);

links(11).name = 'Изодромное';
links(11).W = K*(T*s + 1)/s;

t = 0:0.01:10;
tBefore = (-1:0.1:0).';
yBefore = zeros(size(tBefore));
tPlot = [tBefore; t(:)];

fig = figure('Name', 'Отклики звеньев', 'NumberTitle', 'off');

tabs = uitabgroup(fig);
stepTab = uitab(tabs, 'Title', 'Скачок');
impulseTab = uitab(tabs, 'Title', 'Импульс');

stepLayout = tiledlayout(stepTab, 3, 4);
impulseLayout = tiledlayout(impulseTab, 3, 4);

for i = 1:numel(links)
    link = links(i);

    if i == 2 || i == 10
        h = zeros(size(t));
        w = zeros(size(t));

        if i == 10
            h = K*ones(size(t));
        end
    else
        h = step(link.W, t);
        w = impulse(link.W, t);
    end

    plotTitle = sprintf('%d. %s', i, link.name);

    ax = nexttile(stepLayout, i);
    plot(ax, tPlot, [yBefore; h(:)])
    grid(ax, 'on')
    xlabel(ax, 'Время, с')
    ylabel(ax, 'h(t)')
    title(ax, plotTitle)

    ax = nexttile(impulseLayout, i);
    plot(ax, tPlot, [yBefore; w(:)])
    grid(ax, 'on')
    xlabel(ax, 'время, с')
    ylabel(ax, 'w(t)')
    title(ax, plotTitle)
end

