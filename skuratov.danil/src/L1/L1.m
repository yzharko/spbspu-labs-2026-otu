syms s t K;
K0 = 1.9;
T = 0.7;
T1 = 1.1;
T2 = 1.3;
T1_5 = 3;
K_values = [0.95 1.9 3.8];
dt = 0.01;
output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
W = {K/s; K*s; K; K/(T*s+1); K/(T2^2*s^2+T1_5*s+1);
    K/(T2^2*s^2+T1*s+1); K/(T*s^2+1); K/(s*(T*s+1));
    K*s/(T*s+1); K*(T*s+1); K*(T*s+1)/s};
s_check = tf('s');
models = {1/s_check; tf(0); tf(1); 1/(T*s_check+1);
    1/(T2^2*s_check^2+T1_5*s_check+1); 1/(T2^2*s_check^2+T1*s_check+1);
    1/(T*s_check^2+1); 1/(s_check*(T*s_check+1));
    s_check/(T*s_check+1); tf(1); (T*s_check+1)/s_check};
U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
end_times = [5 0.3 0.3 5 20 25 25 6 5 0.3 5];
labels = compose('K = %.2f', K_values);

for number = 1:numel(W)
    Y_step = W{number} * U_step;
    Y_impulse = W{number} * U_impulse;
    h = simplify(ilaplace(Y_step, s, t));
    w = simplify(ilaplace(Y_impulse, s, t));
    fprintf('\nЗвено %d\nh(t) = %s\nw(t) = %s\n', number, char(h), char(w));
    h_regular = subs(h, [dirac(1,t), dirac(t), heaviside(t)], [0 0 1]);
    w_regular = subs(w, [dirac(1,t), dirac(t), heaviside(t)], [0 0 1]);
    full_h = '';
    full_w = '';
    switch number
        case 2
            full_h = 'K\delta(t)';
            full_w = 'K\delta^{\prime}(t)';
        case 3
            full_w = 'K\delta(t)';
        case 9
            full_w = '(K/T)\delta(t) - (K/T^2)e^{-t/T}';
        case 10
            full_h = 'K + KT\delta(t)';
            full_w = 'K\delta(t) + KT\delta^{\prime}(t)';
        case 11
            full_w = 'K + KT\delta(t)';
    end
    time = (0:round(end_times(number)/dt)).' * dt;
    hs = zeros(numel(time), numel(K_values));
    ws = zeros(size(hs));
    for j = 1:numel(K_values)
        h_k = subs(h_regular, K, K_values(j));
        w_k = subs(w_regular, K, K_values(j));
        h_function = matlabFunction(h_k, 'Vars', t);
        w_function = matlabFunction(w_k, 'Vars', t);
        hs(:,j) = h_function(time);
        ws(:,j) = w_function(time);
    end
    fig = draw_pair(time, hs, ws, labels);
    set(fig, 'Position', [100 100 1000 400]);
    if ~isempty(full_h)
        subplot(1, 2, 1);
        title({'Обычная часть h(t)', ['Полный отклик: ' full_h]}, 'FontSize', 13);
        if number == 2
            title({'Обычная часть h(t) = 0', ['Полный отклик: ' full_h]}, 'FontSize', 13);
        end
    end
    if ~isempty(full_w)
        subplot(1, 2, 2);
        title({'Обычная часть w(t)', ['Полный отклик: ' full_w]}, 'FontSize', 13);
        if ismember(number, [2 3 10])
            title({'Обычная часть w(t) = 0', ['Полный отклик: ' full_w]}, 'FontSize', 13);
        end
    end
    print(fig, fullfile(output_dir, sprintf('link%02d.png', number)), '-dpng', '-r150');
    close(fig);

    fig = figure('Visible', 'off', 'Position', [100 100 1000 400]);
    check_titles = {'h(t): step', 'w(t): impulse'};
    if ~isempty(full_h)
        check_titles{1} = 'Обычная часть h(t): step';
    end
    if ~isempty(full_w)
        check_titles{2} = 'Обычная часть w(t): impulse';
    end
    for panel = 1:numel(check_titles)
        subplot(1, 2, panel);
        hold on;
        for j = 1:numel(K_values)
            W_check = K_values(j) * models{number};
            if panel == 1
                step(W_check, time);
            else
                impulse(W_check, time);
            end
        end
        hold off;
        grid on;
        title(check_titles{panel});
        xlabel('Время t');
        ylabel('Выходной сигнал');
        legend("\fontsize{16}" + string(labels));
        chart = gca;
        chart.LegendLocation = 'northoutside';
        chart.AxesStyle.FontSize = 16;
        chart.Title.FontSize = 16;
        chart.XLabel.FontSize = 16;
        chart.YLabel.FontSize = 16;
        for r = 1:numel(chart.Responses)
            chart.Responses(r).LineWidth = 1.5;
        end
    end
    print(fig, fullfile(output_dir, sprintf('check%02d.png', number)), '-dpng', '-r150');
    close(fig);
end

T_values = T * [1/5 1/2 1 2 5];
time_T = (0:2000).' * dt;
h_T = zeros(numel(time_T), numel(T_values));
w_T = zeros(size(h_T));
for j = 1:numel(T_values)
    W_T = sym(K0)/(sym(T_values(j))*s + 1);
    h_function = matlabFunction(ilaplace(W_T*U_step, s, t), 'Vars', t);
    w_function = matlabFunction(ilaplace(W_T*U_impulse, s, t), 'Vars', t);
    h_T(:,j) = h_function(time_T);
    w_T(:,j) = w_function(time_T);
    model = tf(K0, [T_values(j) 1]);
    h_check = step(model, time_T);
    first98 = find(h_check >= 0.98*K0, 1);
    fprintf('T = %.2f: t98 по формуле = %.6f с; по модели = %.2f с\n', T_values(j), -T_values(j)*log(0.02), time_T(first98));
end
fig = draw_pair(time_T, h_T, w_T, compose('T = %.2f', T_values));
print(fig, fullfile(output_dir, 'exp04_T.png'), '-dpng', '-r150');
close(fig);

T1_values = T1 * [1 1/2 1/5 0];
time_T1 = (0:10000).' * dt;
h_T1 = zeros(numel(time_T1), numel(T1_values));
w_T1 = zeros(size(h_T1));
for j = 1:numel(T1_values)
    W_T1 = sym(K0)/(sym(T2)^2*s^2 + sym(T1_values(j))*s + 1);
    h_function = matlabFunction(ilaplace(W_T1*U_step, s, t), 'Vars', t);
    w_function = matlabFunction(ilaplace(W_T1*U_impulse, s, t), 'Vars', t);
    h_T1(:,j) = h_function(time_T1);
    w_T1(:,j) = w_function(time_T1);
    model = tf(K0, [T2^2 T1_values(j) 1]);
    h_check = step(model, time_T1);
    overshoot = 100*(max(h_check) - K0)/K0;
    last_outside = find(abs(h_check - K0) > 0.02*K0, 1, 'last');
    settling = Inf;
    if T1_values(j) > 0 && last_outside < numel(time_T1)
        settling = time_T1(last_outside + 1);
    end
    fprintf('T1 = %.2f: перерегулирование = %.4f %%; установление в полосе 2%% = %.2f с\n', T1_values(j), overshoot, settling);
end
fig = draw_pair(time_T1, h_T1, w_T1, compose('T1 = %.2f', T1_values));
print(fig, fullfile(output_dir, 'exp06_T1.png'), '-dpng', '-r150');
close(fig);

function fig = draw_pair(time, h, w, labels)
    fig = figure('Visible', 'off', 'Position', [100 100 1000 700]);
    values = {h, w};
    names = {'h(t)', 'w(t)'};
    for panel = 1:numel(values)
        subplot(1, 2, panel);
        plot(time, values{panel}, 'LineWidth', 1.5);
        grid on;
        title(names{panel});
        xlabel('Время t, с');
        ylabel(names{panel});
        xlim([time(1) time(end)]);
        legend(labels, 'Location', 'northoutside', 'NumColumns', min(2, numel(labels)));
        set(gca, 'FontSize', 16);
    end
end
