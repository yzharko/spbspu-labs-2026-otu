clc;
clear;
close all;

s = tf('s');

K  = 0.9;
T  = 0.6;
T1 = 0.8;
T2 = 0.7;
T1_aperiodic = 2.0;

omega = logspace(-2, 2, 2000);

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
         '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
         '9. Реальное дифференцирующее', '10. Форсирующее 1-го порядка', ...
         '11. Изодромное'};

W = {K/s, ...
     K*s, ...
     tf(K), ...
     K/(T*s + 1), ...
     K/(T2^2*s^2 + T1_aperiodic*s + 1), ...
     K/(T2^2*s^2 + T1*s + 1), ...
     K/(T*s^2 + 1), ...
     K/(s*(T*s + 1)), ...
     K*s/(T*s + 1), ...
     K*(T*s + 1), ...
     K*(T*s + 1)/s};

figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for link = 1:numel(W)
    [real_part, imag_part] = nyquist(W{link}, omega);
    real_part = squeeze(real_part);
    imag_part = squeeze(imag_part);

    if link == 7
        resonance_index = find(omega > 1/sqrt(T), 1);
        real_part(resonance_index) = NaN;
        imag_part(resonance_index) = NaN;
    end

    nexttile(link);
    if link == 3
        plot(real_part(1), imag_part(1), 'o');
    else
        plot(real_part, imag_part, 'b', real_part, -imag_part, 'b--');
    end
    if link == 1 || link == 2
        xlim([-1 1]);
    elseif link == 7
        xlim([-5*K 5*K]);
        ylim([-K K]);
    elseif link == 10 || link == 11
        xlim([0 2*max(abs(real_part))]);
    end
    grid on;
    title(names{link}, 'FontSize', 11);
    xlabel('Re W');
    ylabel('Im W');
end

display_figure(W(:), cell(11, 1), names, omega);


K_values = [0.5 0.9 2 3 5];

models_K = cell(11, 5);
labels_K = cell(11, 5);
for link = 1:numel(models_K)
    for k = 1:numel(K_values)
        models_K{link, k} = W{link}*K_values(k)/K;
        labels_K{link, k} = sprintf('K = %.2g', K_values(k));
    end
end
display_figure(models_K, labels_K, names, omega);


T_titles = {'4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
            '5. Апериодическое 2-го порядка', '6. Колебательное', ...
            '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
            '9. Реальное дифференцирующее', '10. Форсирующее 1-го порядка', '11. Изодромное'};
T_names = {'T', 'T1', 'T2', 'T1', 'T2', 'T', 'T', 'T', 'T', 'T'};

T_values = [0.1 0.3 0.6 1.5 3];
T_param = [T_values;
           1.4 1.7 2 3 5;
           0.2 0.3 0.5 0.7 1;
           0.2 0.5 0.8 1.1 1.3;
           0.45 0.6 0.7 1 1.5;
           T_values;
           T_values;
           T_values;
           T_values;
           T_values];

models_T = cell(10, 5);
labels_T = cell(10, 5);
for k = 1:numel(models_T)
    models_T{1, k}  = K/(T_param(1, k)*s + 1);
    models_T{2, k}  = K/(T2^2*s^2 + T_param(2, k)*s + 1);
    models_T{3, k}  = K/(T_param(3, k)^2*s^2 + T1_aperiodic*s + 1);
    models_T{4, k}  = K/(T2^2*s^2 + T_param(4, k)*s + 1);
    models_T{5, k}  = K/(T_param(5, k)^2*s^2 + T1*s + 1);
    models_T{6, k}  = K/(T_param(6, k)*s^2 + 1);
    models_T{7, k}  = K/(s*(T_param(7, k)*s + 1));
    models_T{8, k}  = K*s/(T_param(8, k)*s + 1);
    models_T{9, k}  = K*(T_param(9, k)*s + 1);
    models_T{10, k} = K*(T_param(10, k)*s + 1)/s;
    for row = 1:10
        labels_T{row, k} = sprintf('%s = %.2g', T_names{row}, T_param(row, k));
    end
end
display_figure(models_T, labels_T, T_titles, omega);


corner_5 = sort(-roots([T2^2 T1_aperiodic 1]));
damping_6 = T1/(2*T2);
fprintf('20lg(K) = %.3f дБ\n', 20*log10(K));
fprintf('Звенья 4, 8-11: omega_c = 1/T = %.3f рад/с\n', 1/T);
fprintf('Звено 5: omega_c = %.3f и %.3f рад/с\n', corner_5(1), corner_5(2));
fprintf('Звено 6: zeta = %.3f, omega_0 = %.3f, omega_p = %.3f рад/с\n', ...
    damping_6, 1/T2, sqrt(1 - 2*damping_6^2)/T2);
fprintf('Звено 7: omega_0 = 1/sqrt(T) = %.3f рад/с\n', 1/sqrt(T));


function display_figure(models, labels, titles, omega)
    figure_L = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    figure_phi = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for panel = 1:size(models, 1)
        for curve = 1:size(models, 2)
            [magnitude, phase] = bode(models{panel, curve}, omega);

            figure(figure_L);
            nexttile(panel);
            semilogx(omega, 20*log10(squeeze(magnitude)), 'LineWidth', 1.2);
            hold on;

            figure(figure_phi);
            nexttile(panel);
            semilogx(omega, squeeze(phase), 'LineWidth', 1.2);
            hold on;
        end

        figure(figure_L);
        nexttile(panel);
        grid on;
        title(titles{panel}, 'FontSize', 11);
        xlabel('\omega, рад/с');
        ylabel('L, дБ');

        figure(figure_phi);
        nexttile(panel);
        grid on;
        title(titles{panel}, 'FontSize', 11);
        xlabel('\omega, рад/с');
        ylabel('\psi, град');

        if ~isempty(labels{panel, 1})
            figure(figure_L);
            legend(labels(panel, :), 'Location', 'best', 'FontSize', 8);
            figure(figure_phi);
            legend(labels(panel, :), 'Location', 'best', 'FontSize', 8);
        end
    end
end
