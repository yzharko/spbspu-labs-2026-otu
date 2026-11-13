clear;
close all;
clc;

K = 0.5;
T = 0.1;
T1 = 0.3;
T2 = 0.4;

T1a = 0.9;

s = tf('s');
conservative_model = K / (T*s^2 + 1);
links = [
    struct('name', 'Интегрирующее', ...
        'model', K / s);
    struct('name', 'Идеальное дифференцирующее', ...
        'model', K * s);
    struct('name', 'Усилительное', ...
        'model', tf(K));
    struct('name', 'Апериодическое первого порядка', ...
        'model', K / (T*s + 1));
    struct('name', 'Апериодическое второго порядка', ...
        'model', K / (T2^2*s^2 + T1a*s + 1));
    struct('name', 'Колебательное', ...
        'model', K / (T2^2*s^2 + T1*s + 1));
    struct('name', 'Консервативное', ...
        'model', conservative_model);
    struct('name', 'Инерционное интегрирующее', ...
        'model', K / (s*(T*s + 1)));
    struct('name', 'Реальное дифференцирующее', ...
        'model', K*s / (T*s + 1));
    struct('name', 'Форсирующее', ...
        'model', K * (T*s + 1));
    struct('name', 'Изодромное', ...
        'model', K * (T*s + 1) / s)
];

frequency_range = {1e-2, 1e3};
K_values = K * [1/2, 1, 2];
plot_options = bodeoptions;
plot_options.XLabel.String = '\omega';
plot_options.YLabel.String = {'L', 'phi'};

bode_window = figure('Name', 'ЛАЧХ и ЛФЧХ');
nyquist_window = figure('Name', 'АФЧХ');
gain_window = figure('Name', 'Влияние K');
time_window = figure('Name', 'Влияние постоянных времени');
set([bode_window nyquist_window gain_window time_window], ...
    'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85], 'Color', 'w');

for i = 1:numel(links)
    model = links(i).model;
    link_title = sprintf('%d. %s', i, links(i).name);

    group_plot(bode_window, i, numel(links));
    bode(model, frequency_range, plot_options);
    grid on;
    title(link_title);

    group_plot(nyquist_window, i, numel(links));
    if isstatic(model)
        nyquist(model, 'o', frequency_range);
    elseif isequal(model, conservative_model)
        nyquist(model, '.', frequency_range);
        axis([-2 2 -0.5 0.5]);
    else
        nyquist(model, frequency_range);
    end
    grid off;
    title(link_title);

    gain_models = cell(1, numel(K_values));
    for j = 1:numel(K_values)
        gain_models{j} = (K_values(j) / K) * model;
    end

    group_plot(gain_window, i, numel(links));
    bode(gain_models{:}, frequency_range, plot_options);
    grid on;
    title(link_title);
end

figure(gain_window);
legend(compose('K = %g', K_values), 'Orientation', 'horizontal', ...
    'Units', 'normalized', 'Position', [0.25 0.02 0.5 0.025]);

T_values = T * [1/2, 1, 2];
experiments = [
    struct('name', '4. Апериодическое I: T', ...
        'values', T_values, ...
        'make_model', @(value) K / (value*s + 1));
    struct('name', '7. Консервативное: T', ...
        'values', T_values, ...
        'make_model', @(value) K / (value*s^2 + 1));
    struct('name', '8. Инерционное интегрирующее: T', ...
        'values', T_values, ...
        'make_model', @(value) K / (s*(value*s + 1)));
    struct('name', '9. Реальное дифференцирующее: T', ...
        'values', T_values, ...
        'make_model', @(value) K*s / (value*s + 1));
    struct('name', '10. Форсирующее: T', ...
        'values', T_values, ...
        'make_model', @(value) K * (value*s + 1));
    struct('name', '11. Изодромное: T', ...
        'values', T_values, ...
        'make_model', @(value) K * (value*s + 1) / s);
    struct('name', '5. Апериодическое II: T1a', ...
        'values', [0.8 T1a 1.2], ...
        'make_model', @(value) K / (T2^2*s^2 + value*s + 1));
    struct('name', '5. Апериодическое II: T2', ...
        'values', [0.2 T2 0.45], ...
        'make_model', @(value) K / (value^2*s^2 + T1a*s + 1));
    struct('name', '6. Колебательное: T1', ...
        'values', [0.15 T1 0.75], ...
        'make_model', @(value) K / (T2^2*s^2 + value*s + 1));
    struct('name', '6. Колебательное: T2', ...
        'values', [0.2 T2 0.6], ...
        'make_model', @(value) K / (value^2*s^2 + T1*s + 1))
];

for i = 1:numel(experiments)
    experiment = experiments(i);
    values = experiment.values;
    make_model = experiment.make_model;
    models = cell(1, numel(values));

    for j = 1:numel(values)
        models{j} = make_model(values(j));
    end

    group_plot(time_window, i, numel(experiments));
    bode(models{:}, frequency_range, plot_options);
    grid on;
    title(experiment.name);
end

figure(time_window);
legend({'Меньшее значение', 'Исходное значение', 'Большее значение'}, ...
    'Orientation', 'horizontal', 'Units', 'normalized', ...
    'Position', [0.15 0.02 0.7 0.025]);

function group_plot(window, index, count)
    figure(window);
    n = ceil(sqrt(count));
    subplot(n, n, index);
end
