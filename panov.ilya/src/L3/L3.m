clc;
clear;
close all;

s = tf('s');

K0 = 1.1;
T0 = 0.2;
T1 = 0.3;
T2 = 0.5;
T1_5 = 1.2;

omega = logspace(-2, 3, 2000);

W = {
    K0/s
    K0*s
    tf(K0)
    K0/(T0*s + 1)
    K0/(T2^2*s^2 + T1_5*s + 1)
    K0/(T2^2*s^2 + T1*s + 1)
    K0/(T0*s^2 + 1)
    K0/(s*(T0*s + 1))
    K0*s/(T0*s + 1)
    K0*(T0*s + 1)
    K0*(T0*s + 1)/s
};

K_values = [0.5 1 2 3 5];
T_values = [0.1 0.5 1 2 5];

WK = cell(length(W), length(K_values));
labelsK = cell(1, length(K_values));
titles = cell(length(W), 1);

for n = 1:length(W)
    titles{n} = sprintf('W_{%d}', n);
    for i = 1:length(K_values)
        K = K_values(i);
        WK{n, i} = (K/K0)*W{n};
        labelsK{i} = sprintf('K = %.2g', K);
    end
end

time_values = [
    T_values
    1 1.2 1.5 2 3
    0.1 0.2 0.3 0.5 0.6
    0.05 0.15 0.3 0.5 0.8
    0.2 0.3 0.5 1 2
    T_values
    T_values
    T_values
    T_values
    T_values
];

numbersT = [4 5 5 6 6 7 8 9 10 11];
namesT = {'T', 'T1', 'T2', 'T1', 'T2', 'T', 'T', 'T', 'T', 'T'};
WT = cell(length(numbersT), length(T_values));
labelsT = cell(size(WT));
titlesT = cell(length(numbersT), 1);

for i = 1:length(T_values)
    WT{1, i} = K0/(time_values(1, i)*s + 1);
    WT{2, i} = K0/(T2^2*s^2 + time_values(2, i)*s + 1);
    WT{3, i} = K0/(time_values(3, i)^2*s^2 + T1_5*s + 1);
    WT{4, i} = K0/(T2^2*s^2 + time_values(4, i)*s + 1);
    WT{5, i} = K0/(time_values(5, i)^2*s^2 + T1*s + 1);
    WT{6, i} = K0/(time_values(6, i)*s^2 + 1);
    WT{7, i} = K0/(s*(time_values(7, i)*s + 1));
    WT{8, i} = K0*s/(time_values(8, i)*s + 1);
    WT{9, i} = K0*(time_values(9, i)*s + 1);
    WT{10, i} = K0*(time_values(10, i)*s + 1)/s;
    for n = 1:length(numbersT)
        labelsT{n, i} = sprintf('%s = %.2g', namesT{n}, time_values(n, i));
        titlesT{n} = sprintf('W_{%d}: %s', numbersT(n), namesT{n});
    end
end

for group = 1:3
    if group == 1
        models = W;
        numbers = 1:length(W);
        plot_titles = titles;
        group_name = '';
        gain_scale = K0;
    elseif group == 2
        models = WK;
        numbers = 1:length(W);
        plot_titles = titles;
        labels = repmat(labelsK, length(W), 1);
        group_name = ' - Влияние K';
        gain_scale = max(K_values);
    else
        models = WT;
        numbers = numbersT;
        plot_titles = titlesT;
        labels = labelsT;
        group_name = ' - Влияние постоянных времени';
        gain_scale = K0;
    end

    figure_L = figure('Name', ['ЛАЧХ' group_name], 'NumberTitle', 'off', 'Position', [70 70 1200 800]);
    figure_P = figure('Name', ['ЛФЧХ' group_name], 'NumberTitle', 'off', 'Position', [90 90 1200 800]);
    figure_A = figure('Name', ['АФЧХ' group_name], 'NumberTitle', 'off', 'Position', [110 110 1200 800]);
    colors = lines(size(models, 2));

    for n = 1:size(models, 1)
        for i = 1:size(models, 2)
            [mag, phase] = bode(models{n, i}, omega);
            [re, im] = nyquist(models{n, i}, omega);
            mag = mag(:);
            phase = phase(:);
            re = re(:);
            im = im(:);

            if numbers(n) == 7
                if group == 3
                    resonance = 1/sqrt(T_values(i));
                else
                    resonance = 1/sqrt(T0);
                end
                j = find(omega > resonance, 1);
                mag(j) = NaN;
                phase(j) = NaN;
                re(j) = NaN;
                im(j) = NaN;
            end

            figure(figure_L);
            subplot(3, 4, n);
            semilogx(omega, 20*log10(mag), 'Color', colors(i, :));
            hold on;
            xlabel('\omega, рад/с');
            ylabel('L(\omega), дБ');

            figure(figure_P);
            subplot(3, 4, n);
            semilogx(omega, phase, 'Color', colors(i, :));
            hold on;
            xlabel('\omega, рад/с');
            ylabel('\phi(\omega), град');

            figure(figure_A);
            subplot(3, 4, n);
            if numbers(n) == 3
                plot(re(1), im(1), 'o', 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
                hold on;
            else
                plot(re, im, 'Color', colors(i, :));
                hold on;
                plot(re, -im, 'Color', colors(i, :), 'HandleVisibility', 'off');
            end
            xlabel('Re W(j\omega)');
            ylabel('Im W(j\omega)');

            if numbers(n) == 1 || numbers(n) == 2
                xlim([-1 1]);
            elseif numbers(n) == 7
                xlim([-5*gain_scale 5*gain_scale]);
                ylim([-gain_scale gain_scale]);
            elseif numbers(n) == 10 || numbers(n) == 11
                xlim([0 2*max(abs(re))]);
            elseif numbers(n) ~= 8
                axis equal;
            end
        end

        for f = [figure_L figure_P figure_A]
            figure(f);
            subplot(3, 4, n);
            grid on;
            set(gca, 'FontSize', 9);
            title(plot_titles{n});
            if f ~= figure_A
                xlim([omega(1) omega(end)]);
                xticks(10.^(-2:3));
            end
            if group > 1
                legend(labels(n, :), 'Location', 'best', 'FontSize', 6, 'NumColumns', 2);
            end
            hold off;
        end
    end
end

fprintf('\nСопрягающие частоты звена №4 при K = %.2f:\n', K0);
for i = 1:length(T_values)
    T = T_values(i);
    omega_c = 1/T;
    fprintf('T = %.2f c, omega_c = %.3f рад/с\n', T, omega_c);
end
