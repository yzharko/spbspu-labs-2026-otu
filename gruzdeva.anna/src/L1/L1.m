clear; clc; close all;

K = 0.5; T = 0.1; T1 = 0.3; T2 = 0.4;
T1a = 1.0;

names = {'1. Интегрирующее', '2. Идеальное дифф.', '3. Усилительное', ...
    '4. Апериодическое 1-го', '5. Апериодическое 2-го', '6. Колебательное', ...
    '7. Консервативное', '8. Инерционное интегр.', '9. Реальное дифф.', ...
    '10. Форсирующее', '11. Изодромное'};

syms s t;

W = {
    K / s;
    K * s;
    K;
    K / (T*s + 1);
    K / (T2^2 * s^2 + T1a * s + 1);
    K / (T2^2 * s^2 + T1 * s + 1);
    K / (T*s^2 + 1);
    K / (s * (T*s + 1));
    (K*s) / (T*s + 1);
    K * (T*s + 1);
    K * (T*s + 1) / s
};

h_sym = cell(1, 11);
w_sym = cell(1, 11);

U_step = 1/s;
U_impulse = 1;

for i = 1:11
    fprintf('\n%s\n', names{i});

    Y_step = W{i} * U_step;
    h_sym{i} = simplify(ilaplace(Y_step, s, t));
    fprintf('h(t) = %s\n', char(h_sym{i}));

    Y_impulse = W{i} * U_impulse;
    w_sym{i} = simplify(ilaplace(Y_impulse, s, t));
    fprintf('w(t) = %s\n', char(w_sym{i}));
end

t_anal = (0.01:0.01:5)';

h_singular_note = repmat({''}, 1, 11);
h_singular_note{2} = 'Дополнительно: K*delta(t) в t=0';
h_singular_note{10} = 'Дополнительно: K*T*delta(t) в t=0';

w_singular_note = repmat({''}, 1, 11);
w_singular_note{2} = 'Сингулярная часть: K*delta''(t)';
w_singular_note{3} = 'Сингулярная часть: K*delta(t)';
w_singular_note{9} = 'Сингулярная часть: (K/T)*delta(t)';
w_singular_note{10} = 'K*delta(t) + K*T*delta''(t)';
w_singular_note{11} = 'Сингулярная часть: K*T*delta(t)';

figure('Name', 'Аналитические переходные характеристики h(t)', ...
    'Position', [50 50 1400 800]);
for i = 1:11
    subplot(3, 4, i);
    h_values = double(subs(h_sym{i}, t, t_anal));
    h_values = real(h_values(:));
    plot(t_anal, h_values, 'LineWidth', 1.3);
    grid on;
    title(names{i}, 'Interpreter', 'none');
    xlabel('t, c');
    ylabel('h(t)');

    if ~isempty(h_singular_note{i})
        text(0.04, 0.94, h_singular_note{i}, ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Interpreter', 'none', ...
            'BackgroundColor', 'w');
    end
end
sgtitle('Переходные характеристики, полученные обратным преобразованием Лапласа', ...
    'FontSize', 13);

figure('Name', 'Аналитические импульсные характеристики w(t)', ...
    'Position', [70 70 1400 800]);
for i = 1:11
    subplot(3, 4, i);
    w_values = double(subs(w_sym{i}, t, t_anal));
    w_values = real(w_values(:));
    plot(t_anal, w_values, 'LineWidth', 1.3);
    grid on;
    title(names{i}, 'Interpreter', 'none');
    xlabel('t, c');
    ylabel('w(t)');

    if ~isempty(w_singular_note{i})
        text(0.04, 0.94, w_singular_note{i}, ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 8, 'Interpreter', 'none', ...
            'BackgroundColor', 'w');
    end
end
sgtitle('Импульсные характеристики, полученные обратным преобразованием Лапласа', ...
    'FontSize', 13);

t_vec = 0:0.01:5;
W_tf = {tf(K, [1 0]), [], tf(K, 1), tf(K, [T 1]), ...
    tf(K, [T2^2, T1a, 1]), tf(K, [T2^2, T1, 1]), tf(K, [T 0 1]), ...
    tf(K, [T 1 0]), tf([K 0], [T 1]), [], tf([K*T K], [1 0])};

figure('Name', 'Переходные характеристики h(t), Control System Toolbox', ...
    'Position', [50 50 1400 800]);
for i = 1:11
    subplot(3, 4, i);

    if isempty(W_tf{i})
        axis([0 5 -1 1]);
        grid on;
        text(2.5, 0, 'step не строит', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        title(names{i}, 'Interpreter', 'none');
        xlabel('t, c'); ylabel('h(t)');
    else
        step(W_tf{i}, t_vec);
        grid on;
        title(names{i}, 'Interpreter', 'none');
        xlabel('t, c'); ylabel('h(t)');
    end
end
sgtitle('Переходные функции h(t): Control System Toolbox', 'FontSize', 14);

figure('Name', 'Импульсные характеристики w(t), Control System Toolbox', ...
    'Position', [50 50 1400 800]);
for i = 1:11
    subplot(3, 4, i);

    if isempty(W_tf{i})
        axis([0 5 -1 1]);
        grid on;
        text(2.5, 0, 'impulse не строит', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        title(names{i}, 'Interpreter', 'none');
        xlabel('t, c'); ylabel('w(t)');
    else
        impulse(W_tf{i}, t_vec);
        grid on;
        title(names{i}, 'Interpreter', 'none');
        xlabel('t, c'); ylabel('w(t)');
    end
end
sgtitle('Импульсные функции w(t): Control System Toolbox', 'FontSize', 14);

figure('Name', 'Влияние K', 'Position', [100 100 800 500]);
hold on; grid on;
K_vals = [0.25, 0.5, 1.0, 1.5, 2.0];
for i = 1:length(K_vals)
    step(tf(K_vals(i), [T 1]), 1);
end
legend('K=0.25', 'K=0.5', 'K=1.0', 'K=1.5', 'K=2.0', ...
    'Location', 'southeast');
title('Звено №4: h(t) при изменении K', 'FontSize', 12);
xlabel('t, c'); ylabel('h(t)');
hold off;

figure('Name', 'Влияние T', 'Position', [150 150 800 500]);
hold on; grid on;
T_vals = [T/5, T/2, T, T*2, T*5];
for i = 1:length(T_vals)
    step(tf(K, [T_vals(i) 1]), 2);
end
legend('T/5', 'T/2', 'T', '2T', '5T', 'Location', 'southeast');
title('Звено №4: h(t) при изменении T', 'FontSize', 12);
xlabel('t, c'); ylabel('h(t)');
hold off;

figure('Name', 'Влияние T1', 'Position', [200 200 800 500]);
hold on; grid on;
T1_vals = [0.3, 0.2, 0.1, 0.05];
for i = 1:length(T1_vals)
    step(tf(K, [T2^2, T1_vals(i), 1]), 35);
end
legend('T1=0.3', 'T1=0.2', 'T1=0.1', 'T1=0.05', ...
    'Location', 'northeast');
title('Звено №6: h(t) при уменьшении T1', 'FontSize', 12);
xlabel('t, c'); ylabel('h(t)');
hold off;

disp('Done');

t_check = (0:0.01:5)';

err_h = nan(11,1);
err_w = nan(11,1);
notes = strings(11,1);

for i = 1:11

    notes(i) = "Сравнение выполнено";

    if isempty(W_tf{i})
        notes(i) = "Численное моделирование невозможно";
        continue;
    end

    [h_num, t_num] = step(W_tf{i}, t_check);

    t_num = t_num(:);
    h_num = squeeze(h_num);
    h_num = h_num(:);

    h_anal = double(subs(h_sym{i}, t, t_num));
    h_anal = real(h_anal(:));
    h_num = real(h_num);

    idx = t_num > 0;

    err_h(i) = max(abs(h_anal(idx) - h_num(idx)));

    if i ~= 3
        [w_num, tw_num] = impulse(W_tf{i}, t_check);

        tw_num = tw_num(:);
        w_num = squeeze(w_num);
        w_num = real(w_num(:));

        w_anal = double(subs(w_sym{i}, t, tw_num));
        w_anal = real(w_anal(:));

        idx_w = tw_num > 0;

        err_w(i) = max(abs(w_anal(idx_w) - w_num(idx_w)));
    end

    if i == 2
        notes(i) = "Идеальное дифференцирующее звено: сингулярные характеристики";

    elseif i == 3
        notes(i) = "Импульсная характеристика содержит только дельта-функцию";

    elseif i == 9 || i == 11
        notes(i) = "Сравниваются регулярные части при t > 0";

    end

end

ComparisonTable = table( ...
    (1:11)', ...
    string(names(:)), ...
    err_h, ...
    err_w, ...
    notes, ...
    'VariableNames', ...
    {'Number','Link','Error_h','Error_w','Comment'});

disp('ТАБЛИЦА СРАВНЕНИЯ РЕЗУЛЬТАТОВ:');
disp(ComparisonTable);

