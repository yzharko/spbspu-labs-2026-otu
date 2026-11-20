K = 4;
T = 1.5;
T1 = 1.6;
T2 = 2;
T1_5 = 5; % Для звена 5 требуется T1 >= 2*T2 = 4, а в варианте T1 = 1.6.

outDir = fullfile(fileparts(mfilename('fullpath')), 'results_L1');
if ~isfolder(outDir)
    mkdir(outDir);
end
png = @(name) fullfile(outDir, [name '.png']);

syms s t
p = tf('s');
first_order = @(gain,tau) gain / (tau*p + 1);
second_order = @(t1,t2) K / (t2^2*p^2 + t1*p + 1);
conservative = @(tau) K / (tau*p^2 + 1);
inertial_integrator = @(tau) K / (p*(tau*p + 1));
differentiator = @(tau) K*p / (tau*p + 1);
isodromic = @(tau) K*(tau*p + 1) / p;

% Название, W(s) для расчёта по Лапласу, модель tf для проверки.
links = {
    'Интегрирующее', K/s, K/p;
    'Идеальное дифференцирующее', K*s, K*p;
    'Усилительное', K, tf(K);
    'Апериодическое первого порядка', K/(T*s + 1), first_order(K,T);
    'Апериодическое второго порядка', K/(T2^2*s^2 + T1_5*s + 1), second_order(T1_5,T2);
    'Колебательное', K/(T2^2*s^2 + T1*s + 1), second_order(T1,T2);
    'Консервативное', K/(T*s^2 + 1), conservative(T);
    'Инерционное интегрирующее', K/(s*(T*s + 1)), inertial_integrator(T);
    'Реальное дифференцирующее', K*s/(T*s + 1), differentiator(T);
    'Форсирующее', K*(T*s + 1), K*(T*p + 1);
    'Изодромное', K*(T*s + 1)/s, isodromic(T)
};
names = links(:,1);
W = vertcat(links{:,2});
G = links(:,3);

%% Переходные и импульсные переходные функции
U_step = laplace(heaviside(t), t, s)
U_impulse = laplace(dirac(t), t, s)
Y_step = simplify(W * U_step);
Y_impulse = simplify(W * U_impulse);
h = simplify(ilaplace(Y_step, s, t));
w = simplify(ilaplace(Y_impulse, s, t));

fid = fopen(fullfile(outDir, 'formulas.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, 'K=%g; T=%g; T1=%g; T2=%g; T1_5=%g\n', K,T,T1,T2,T1_5);
fprintf(fid, 'U_step=%s; U_impulse=%s\n', char(U_step), char(U_impulse));
for j = 1:numel(W)
    fprintf(fid, '\n%d. %s\nW=%s\nY_step=%s\nh=%s\nY_impulse=%s\nw=%s\n', ...
        j,names{j},char(W(j)),char(Y_step(j)),char(h(j)),char(Y_impulse(j)),char(w(j)));
end
fclose(fid);

%% Графики 11 звеньев и проверка с помощью step и impulse
notes = repmat({''}, numel(W), 2);
notes(2,:) = {'h(t)=4 delta(t)', 'w(t)=4 delta''(t)'};
notes{3,2} = 'w(t)=4 delta(t)';
notes{9,2} = 'При t=0: (8/3) delta(t)';
notes(10,:) = {'При t=0: 6 delta(t)', 'w(t)=4 delta(t)+6 delta''(t)'};
notes{11,2} = 'При t=0: 6 delta(t)';
for j = 1:numel(W)
    time = (0:0.01:10)';
    if ismember(j, [5 6 7])
        time = (0:0.02:40)';
    end
    positive_time = time(2:end);
    h_values = sample_response(h(j), t, positive_time);
    w_values = sample_response(w(j), t, positive_time);
    labels = {'Расчёт по Лапласу'};
    if ~ismember(j, [2 10])
        h_check = step(G{j}, time);
        w_check = impulse(G{j}, time);
        h_values(:,2) = h_check(2:end);
        w_values(:,2) = w_check(2:end);
        labels{2} = 'Проверка step/impulse';
    end
    plot_pair(positive_time, h_values, w_values, labels, ...
        sprintf('%d. %s, t>0', j, names{j}), notes(j,:), png(sprintf('link_%02d', j)));
end

% Для звеньев 2 и 10 проверяем отклики по y2 = K*u' и y10 = K*(u + T*u').
h2_check = K * diff(heaviside(t), t);
w2_check = K * diff(dirac(t), t);
h10_check = K * (heaviside(t) + T*diff(heaviside(t), t));
w10_check = K * (dirac(t) + T*diff(dirac(t), t));
disp('Проверка звеньев 2 и 10: четыре нуля означают совпадение.');
disp(simplify(laplace([h2_check; w2_check; h10_check; w10_check], t, s) ...
    - [Y_step(2); Y_impulse(2); Y_step(10); Y_impulse(10)]));

%% Влияние K, T, T1 и T2
sweep('K', [2 4 8], @(k) first_order(k,T), 0:0.02:12, ...
    'Влияние K, звено 4', png('vary_K'));
sweep('T', T*[1/5 1/2 1 2 5], @(tau) first_order(K,tau), 0:0.02:40, ...
    'Звено 4: изменение T', png('vary_T4'));
sweep('T1', [4 T1_5 8], @(t1) second_order(t1,T2), 0:0.02:40, ...
    sprintf('Звено 5: T2=%g', T2), png('vary_T1'));
sweep('T1', [T1 T1/2 T1/5 0], @(t1) second_order(t1,T2), 0:0.05:100, ...
    sprintf('Звено 6: T2=%g', T2), png('decrease_T1'));
sweep('T2', [1 2 2.5], @(t2) second_order(T1_5,t2), 0:0.02:30, ...
    sprintf('Звено 5: T1=%g', T1_5), png('vary_T2_05'));
sweep('T2', [1 2 4], @(t2) second_order(T1,t2), 0:0.05:100, ...
    sprintf('Звено 6: T1=%g', T1), png('vary_T2_06'));

T_values = T*[1/2 1 2];
time = 0:0.02:20;

% Для звеньев 7–11 выводим графики при t > 0.
sweep('T', T_values, conservative, time, 'Звено 7: изменение T', png('vary_T_07'), 2);
sweep('T', T_values, inertial_integrator, time, 'Звено 8: изменение T', png('vary_T_08'), 2);
sweep('T', T_values, differentiator, time, 'Звено 9: изменение T', png('vary_T_09'), ...
    2, 'При t=0: (K/T) delta(t)');
sweep('T', T_values, isodromic, time, 'Звено 11: изменение T', png('vary_T_11'), ...
    2, 'При t=0: K*T delta(t)');

% У звена 10 меняются коэффициенты при delta в h и при delta' в w.
fprintf('Звено 10: T=%s; K*T=%s.\n', mat2str(T_values), mat2str(K*T_values));

function sweep(name, values, make_model, time, caption, file, first, impulse_note)
    if nargin < 7
        first = 1;
    end
    if nargin < 8
        impulse_note = '';
    end
    time = time(:);
    h = zeros(numel(time), numel(values));
    w = h;
    labels = cell(size(values));
    for j = 1:numel(values)
        model = make_model(values(j));
        h(:,j) = step(model, time);
        w(:,j) = impulse(model, time);
        labels{j} = sprintf('%s=%g', name, values(j));
    end
    rows = first:numel(time);
    plot_pair(time(rows), h(rows,:), w(rows,:), labels, caption, {'', impulse_note}, file);
end

function values = sample_response(formula, t, time)
    formula = subs(formula, [dirac(t), dirac(1,t), heaviside(t)], [0, 0, 1]);
    response = matlabFunction(formula, 'Vars', t);
    values = double(response(time)) + zeros(size(time));
end

function plot_pair(time, h, w, labels, caption, notes, file)
    fig = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 1200 440]);
    values = {h, w};
    ylabels = {'h(t), отн. ед.', 'w(t), отн. ед./с'};
    titles = {'Переходная функция', 'Импульсная переходная функция'};
    for p = 1:2
        ax = subplot(1, 2, p);
        curves = plot(ax, time, values{p}, 'LineWidth', 1.4);
        if numel(labels) == 2 && strcmp(labels{2}, 'Проверка step/impulse')
            curves(2).LineStyle = '--';
        end
        grid(ax, 'on');
        box(ax, 'on');
        set(ax, 'Color', 'white', 'XColor', 'black', 'YColor', 'black', 'FontSize', 14);
        xlabel(ax, 't, с');
        ylabel(ax, ylabels{p});
        title(ax, titles{p}, 'Color', 'black', 'FontSize', 14);
        location = 'best';
        if ~isempty(notes{p})
            location = 'southwest';
            limits = ylim(ax);
            ylim(ax, [limits(1), limits(2) + 0.3*diff(limits)]);
            text(ax, 0.02, 0.97, notes{p}, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
                'Color', 'black', 'FontSize', 13, 'Interpreter', 'none');
        end
        if numel(labels) > 1
            legend(ax, labels, 'Location', location, ...
                'TextColor', 'black', 'Color', 'white', 'FontSize', 13);
        end
    end
    sgtitle(fig, caption, 'Color', 'black', 'FontSize', 16);
    exportgraphics(fig, file, 'Resolution', 200, 'BackgroundColor', 'white');
    close(fig);
end
