clear;
close all;
clc;

K = 0.5;
T = 0.1;
T1 = 0.3;
T2 = 0.4;

syms s real;
s_tf = tf('s');

links = [
    struct('name', 'Интегрирующее', ...
        'W', K / s, ...
        'model', K / s_tf, ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Идеальное дифференцирующее', ...
        'W', K * s, ...
        'model', tf(0), ...
        'step_delta', K, ...
        'impulse_delta', 0);

    struct('name', 'Усилительное', ...
        'W', sym(K), ...
        'model', tf(K), ...
        'step_delta', 0, ...
        'impulse_delta', K);

    struct('name', 'Апериодическое первого порядка', ...
        'W', K / (T * s + 1), ...
        'model', K / (T * s_tf + 1), ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Апериодическое второго порядка', ...
        'W', K / (T2^2 * s^2 + T1 * s + 1), ...
        'model', K / (T2^2 * s_tf^2 + T1 * s_tf + 1), ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Колебательное', ...
        'W', K / (T2^2 * s^2 + T1 * s + 1), ...
        'model', K / (T2^2 * s_tf^2 + T1 * s_tf + 1), ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Консервативное', ...
        'W', K / (T * s^2 + 1), ...
        'model', K / (T * s_tf^2 + 1), ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Инерционное интегрирующее', ...
        'W', K / (s * (T * s + 1)), ...
        'model', K / (s_tf * (T * s_tf + 1)), ...
        'step_delta', 0, ...
        'impulse_delta', 0);

    struct('name', 'Реальное дифференцирующее', ...
        'W', K * s / (T * s + 1), ...
        'model', K * s_tf / (T * s_tf + 1), ...
        'step_delta', 0, ...
        'impulse_delta', K / T);

    struct('name', 'Форсирующее', ...
        'W', K * (T * s + 1), ...
        'model', tf(K), ...
        'step_delta', K * T, ...
        'impulse_delta', K);

    struct('name', 'Изодромное', ...
        'W', K * (T * s + 1) / s, ...
        'model', K * (T * s_tf + 1) / s_tf, ...
        'step_delta', 0, ...
        'impulse_delta', K * T)
];

plot_responses(links);
study_time_constant(K, T);
study_damping(K, T1, T2);

function plot_responses(links)
    disp('Аналитические отклики, сравнение с моделированием и влияние K.');

    syms s t real;
    n = numel(links);
    time_step = 0.002;
    observation_time = 6;
    plot_columns = 3;
    plot_rows = ceil((n + 1) / plot_columns);
    gain_factors = [1 / 2, 1, 2];
    response_titles = {'Регулярные части h(t): изменение K',  'Регулярные части w(t): изменение K'};
    model_labels = {'step при K', 'impulse при K'};

    U_step = laplace(heaviside(t), t, s);
    U_impulse = laplace(dirac(t), t, s);
    time = 0:time_step:observation_time;

    for i = 1:n
        link = links(i);
        W = link.W;
        W_check = link.model;
        a = link.step_delta;
        b = link.impulse_delta;

        Y_step = W * U_step;
        Y_impulse = W * U_impulse;
        h = ilaplace(Y_step, s, t);
        w = ilaplace(Y_impulse, s, t);

        regular = {simplify(h - a * dirac(t)),  simplify(w - b * dirac(t) - a * diff(dirac(t), t))};
        numerical = {step(W_check, time).', impulse(W_check - b, time).'};

        for r = 1:numel(regular)
            expression = subs(regular{r}, heaviside(t), 1);
            values = double(subs(expression, t, time)) + zeros(size(time));

            figure(r);
            subplot(plot_rows, plot_columns, i);
            plot(time, values.' * gain_factors);
            hold on;
            plot(time, numerical{r}, 'k--');
            set_plot_format(link.name);
        end
    end

    for r = 1:numel(response_titles)
        figure(r);
        sgtitle(response_titles{r});
        last_plot = gca;
        legend_cell = subplot(plot_rows, plot_columns, n + 1);
        axis off;
        legend(last_plot, 'K/2', 'K', '2K', model_labels{r},  'Position', legend_cell.Position);
    end
end

function study_time_constant(K, T)
    disp('Влияние T на апериодическое звено первого порядка.');

    time_step = 0.002;
    settling_tolerance = 0.02;
    time_constants_to_observe = 6;
    parameter_figure_number = 3;
    T_values = T * [1 / 5, 1 / 2, 1, 2, 5];
    time_T = 0:time_step:time_constants_to_observe * max(T_values);

    figure(parameter_figure_number);
    subplot(1, 2, 1);
    hold on;

    for T_value = T_values
        plot(time_T, K * (1 - exp(-time_T / T_value)));
    end

    set_plot_format('Первый порядок: изменение T');
    ylabel('h(t)');
    legend(string(T_values));

    settling_T = -log(settling_tolerance) * T_values;
    disp('T, с | Время установления, с');
    disp([T_values.', settling_T.']);
end

function study_damping(K, T1, T2)
    disp('Влияние уменьшения T1 при фиксированном T2.');

    s_tf = tf('s');
    time_step = 0.002;
    observation_time = 30;
    settling_tolerance = 0.02;
    parameter_figure_number = 3;
    T1_values = [T1, T1 / 2, T1 / 5, 0];
    time_T1 = 0:time_step:observation_time;
    overshoot = zeros(size(T1_values));
    settling_T1 = zeros(size(T1_values));

    figure(parameter_figure_number);
    subplot(1, 2, 2);
    hold on;

    for j = 1:numel(T1_values)
        model = K / (T2^2 * s_tf^2 + T1_values(j) * s_tf + 1);
        response = step(model, time_T1);
        plot(time_T1, response);

        if T1_values(j) == 0
            overshoot(j) = 100;
            settling_T1(j) = Inf;
        else
            info = stepinfo(response, time_T1, K,  'SettlingTimeThreshold', settling_tolerance);
            overshoot(j) = info.Overshoot;
            settling_T1(j) = info.SettlingTime;
        end
    end

    set_plot_format('Колебательное звено: уменьшение T1');
    ylabel('h(t)');
    legend(string(T1_values));
    disp('T1, с | Перерегулирование, % | Время установления, с');
    disp([T1_values.', overshoot.', settling_T1.']);
end

function set_plot_format(plot_title)
    grid on;
    title(plot_title);
    xlabel('Время, с');
end
