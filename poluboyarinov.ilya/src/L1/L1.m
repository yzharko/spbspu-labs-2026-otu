clear;
clc;
close all;

K0 = 1.8;
T = 0.6;
T1 = 0.9;
T2 = 1.2;
T11 = 3.0; %чтобы модели 5 и 6 отличались

K_values = [0.9, 1.8, 3.6];
time = 0:0.01:10;
time_check = linspace(0, 10, 1000);
time_long = linspace(0, 40, 1500);
T_values = [T/5, T/2, T, 2*T, 5*T];
T1_values = [T1, T1/2, T1/5, 0];

syms p t K real
s_tf = tf('s');

names = {
    'Интегрирующее'
    'Идеальное дифференцирующее'
    'Усилительное'
    'Апериодическое 1-го порядка'
    'Апериодическое 2-го порядка'
    'Колебательное'
    'Консервативное'
    'Инерционное интегрирующее'
    'Реальное дифференцирующее'
    'Форсирующее'
    'Изодромное'
};

n = numel(names);

W = {
    K/p
    K*p
    K
    K/(T*p + 1)
    K/(T2^2*p^2 + T11*p + 1)
    K/(T2^2*p^2 + T1*p + 1)
    K/(T*p^2 + 1)
    K/(p*(T*p + 1))
    K*p/(T*p + 1)
    K*(T*p + 1)
    K*(T*p + 1)/p
};

Wtf = cell(n, 1);

for i = 1:n
    expression = subs(W{i}, K, K0);
    [numerator, denominator] = numden(expression);

    numerator_coefficients = double(sym2poly(expand(numerator)));
    denominator_coefficients = double(sym2poly(expand(denominator)));

    Wtf{i} = tf(numerator_coefficients, denominator_coefficients);
end

fprintf('K = %.2f, T = %.2f, T1 = %.2f, T2 = %.2f, T11 = %.2f\n', ...
    K0, T, T1, T2, T11);

Y_step = cell(n, 1);
Y_impulse = cell(n, 1);
h = cell(n, 1);
w = cell(n, 1);
R = cell(n, 2);

for i = 1:n
    Y_step{i} = simplify(W{i}/p);
    Y_impulse{i} = simplify(W{i});

    h{i} = simplify(ilaplace(Y_step{i}, p, t));
    w{i} = simplify(ilaplace(Y_impulse{i}, p, t));

    R{i, 1} = simplify(subs(subs(subs( ...
        h{i}, dirac(1, t), 0), dirac(t), 0), heaviside(t), 1));

    R{i, 2} = simplify(subs(subs(subs( ...
        w{i}, dirac(1, t), 0), dirac(t), 0), heaviside(t), 1));

    fprintf('\nЗвено №%d — %s\n', i, names{i});
    fprintf('W(s) = %s\n', char(W{i}));
    fprintf('h(t) = %s\n', char(h{i}));
    fprintf('w(t) = %s\n', char(w{i}));
end
fprintf('\n');

for mode = 1:2
    if mode == 1
        label = 'h(t)';
        figure_name = 'Переходные характеристики';
    else
        label = 'w(t)';
        figure_name = 'Импульсные характеристики';
    end

    figure('Name', figure_name, 'Color', 'w', ...
        'Position', [50, 50, 1450, 900]);

    for i = 1:n
        subplot(4, 3, i);
        hold on;
        grid on;

        singular = (mode == 1 && i == 2) || ...
            (mode == 2 && ismember(i, [1, 2, 3, 10]));

        if ~singular
            for K_plot = K_values
                expression = subs(R{i, mode}, K, K_plot);
                response = double(subs(expression, t, time));

                if isscalar(response)
                    response = repmat(response, size(time));
                end

                plot(time, response, 'LineWidth', 1.3, ...
                    'DisplayName', sprintf('K = %.1f', K_plot));
            end
        end

        title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
        xlabel('Время, с');
        ylabel(label);

        if mode == 1 && i == 2
            text(0.08, 0.85, 'K\delta(t)', ...
                'Units', 'normalized', 'FontSize', 8);
        elseif mode == 2 && i == 2
            text(0.08, 0.85, 'K\delta''(t)', ...
                'Units', 'normalized', 'FontSize', 8);
        elseif mode == 2 && i == 3
            text(0.08, 0.85, 'K\delta(t)', ...
                'Units', 'normalized', 'FontSize', 8);
        elseif mode == 2 && i == 9
            text(0.08, 0.85, '+ (K/T)\delta(t)', ...
                'Units', 'normalized', 'FontSize', 8);
        elseif mode == 2 && i == 11
            text(0.08, 0.85, '+ KT\delta(t)', ...
                'Units', 'normalized', 'FontSize', 8);
        end

        if i == 1
            legend('show', 'Location', 'best');
        end

        hold off;
    end

    sgtitle(figure_name);
end

step_check = [1, 3, 4, 5, 6, 7, 8, 9, 11];

for i = step_check
    expression = subs(R{i, 1}, K, K0);
    analytic_response = double(subs(expression, t, time_check));
    analytic_response = analytic_response(:);

    matlab_response = step(Wtf{i}, time_check);
    matlab_response = squeeze(matlab_response);
    matlab_response = matlab_response(:);

    if numel(matlab_response) ~= numel(analytic_response)
        fprintf('h%d: размер\n', i);
        continue;
    end

    error_value = max(abs(matlab_response - analytic_response));
    fprintf('h%d: %.3e\n', i, error_value);
end

impulse_check = [4, 5, 6, 7, 8, 9, 11];

for i = impulse_check
    expression = subs(R{i, 2}, K, K0);
    analytic_response = double(subs(expression, t, time_check));

    if isscalar(analytic_response)
        analytic_response = repmat(analytic_response, size(time_check));
    end

    if i == 9
        model = -(K0/T)/(T*s_tf + 1);
    elseif i == 11
        model = K0/s_tf;
    else
        model = Wtf{i};
    end

    matlab_response = impulse(model, time_check);
    error_value = max(abs(matlab_response(:) - analytic_response(:)));

    fprintf('w%d: %.3e\n', i, error_value);
end

figure('Name', 'Влияние T', 'Color', 'w');
hold on;
grid on;

for T_plot = T_values
    model = K0/(T_plot*s_tf + 1);

    plot(time, step(model, time), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('T = %.2f', T_plot));
end

title('Влияние T');
xlabel('Время, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

figure('Name', 'Влияние T1', 'Color', 'w');
hold on;
grid on;

for T1_plot = T1_values
    model = K0/(T2^2*s_tf^2 + T1_plot*s_tf + 1);

    plot(time_long, step(model, time_long), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('T1 = %.2f', T1_plot));
end

title('Влияние T1');
xlabel('Время, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;
