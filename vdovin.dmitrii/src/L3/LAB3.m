clc;
clear;
close all;

s = tf('s');
K0 = 1.5;
T0 = 0.4;
T1 = 0.6;
T2 = 0.8;
T1_ap = 2*T2 + T1;
omega = logspace(-2, 3, 1000);

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
    '4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
    '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
    '9. Реальное дифференцирующее', '10. Форсирующее', '11. Изодромное'};

W = {K0/s, K0*s, tf(K0, 1), K0/(T0*s+1), ...
    K0/(T2^2*s^2+T1_ap*s+1), K0/(T2^2*s^2+T1*s+1), ...
    K0/(T0*s^2+1), K0/(s*(T0*s+1)), ...
    K0*s/(T0*s+1), K0*(T0*s+1), K0*(T0*s+1)/s};

fig_names = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ', 'Влияние K: ЛАЧХ', ...
    'Влияние T: ЛАЧХ', 'Влияние T: ЛФЧХ'};
for f = 1:6
    figure(f);
    set(gcf, 'Name', fig_names{f}, 'NumberTitle', 'off', ...
        'Position', [50 50 1500 850]);
end

for n = 1:length(W)
    [mag, phase] = bode(W{n}, omega);
    [re, im] = nyquist(W{n}, omega);
    mag = mag(:); phase = phase(:);
    re = re(:); im = im(:);

    if n == 1 || n == 2
        re(:) = 0;
    elseif n == 3
        re(:) = K0; im(:) = 0;
    elseif n == 7
        im(:) = 0;
        j = find(omega > 1/sqrt(T0), 1);
        re(j) = NaN;
    elseif n == 10
        re(:) = K0;
    elseif n == 11
        re(:) = K0*T0;
    end

    figure(1); subplot(3,4,n);
    semilogx(omega, 20*log10(mag)); grid on;
    title(names{n}); ylabel('дБ');

    figure(2); subplot(3,4,n);
    semilogx(omega, phase); grid on;
    title(names{n}); ylabel('град');

    figure(3); subplot(3,4,n);
    if n == 3
        plot(K0, 0, 'bo', 'MarkerFaceColor', 'b');
    else
        plot(re, im, 'b', flipud(re), -flipud(im), 'b');
    end
    grid on;
    if n == 1 || n == 2
        xlim([-0.3 0.3]);
    elseif n == 3
        xlim([0 3]); ylim([-1 1]);
    elseif n == 7
        ylim([-1 1]);
    elseif n == 10
        xlim([0 3]);
    elseif n == 11
        xlim([0 1.2]);
    end
    title(names{n}); xlabel('Re'); ylabel('Im');
end

K_values = [0.5 1 1.5 2 3];
labelsK = cell(1, length(K_values));
for i = 1:length(K_values)
    labelsK{i} = sprintf('K = %.1f', K_values(i));
end

for n = 1:length(W)
    for i = 1:length(K_values)
        WK = (K_values(i)/K0)*W{n};
        mag = bode(WK, omega);

        figure(4); subplot(3,4,n);
        semilogx(omega, 20*log10(mag(:))); hold on;
    end
    figure(4); subplot(3,4,n);
    grid on; title(names{n}); ylabel('дБ');
end

figure(4); subplot(3,4,12); hold on;
for i = 1:length(K_values)
    plot(NaN, NaN);
end
axis off; legend(labelsK, 'Location', 'northwest');

T_values = [0.1 0.2 0.4 0.8 1.6];
T1_ap_values = [1.6 1.8 2.2 2.6 3];
T2_ap_values = [0.2 0.4 0.6 0.8 1];
T1_osc_values = [0.2 0.4 0.6 0.8 1.2];
T2_osc_values = [0.4 0.6 0.8 1 1.2];

WT = cell(10, length(T_values));
for i = 1:length(T_values)
    T = T_values(i);
    WT{1,i} = K0/(T*s+1);
    WT{2,i} = K0/(T2^2*s^2+T1_ap_values(i)*s+1);
    WT{3,i} = K0/(T2_ap_values(i)^2*s^2+T1_ap*s+1);
    WT{4,i} = K0/(T2^2*s^2+T1_osc_values(i)*s+1);
    WT{5,i} = K0/(T2_osc_values(i)^2*s^2+T1*s+1);
    WT{6,i} = K0/(T*s^2+1);
    WT{7,i} = K0/(s*(T*s+1));
    WT{8,i} = K0*s/(T*s+1);
    WT{9,i} = K0*(T*s+1);
    WT{10,i} = K0*(T*s+1)/s;
end

T_all = [T_values; T1_ap_values; T2_ap_values; ...
    T1_osc_values; T2_osc_values; T_values; T_values; ...
    T_values; T_values; T_values];
namesT = {'№4: T', '№5: T1', '№5: T2', '№6: T1', '№6: T2', ...
    '№7: T', '№8: T', '№9: T', '№10: T', '№11: T'};
paramsT = {'T', 'T1', 'T2', 'T1', 'T2', 'T', 'T', 'T', 'T', 'T'};

for n = 1:10
    labelsT = cell(1, length(T_values));
    for i = 1:length(T_values)
        [mag, phase] = bode(WT{n,i}, omega);
        labelsT{i} = sprintf('%s = %.2g', paramsT{n}, T_all(n,i));

        figure(5); subplot(3,4,n);
        semilogx(omega, 20*log10(mag(:))); hold on;

        figure(6); subplot(3,4,n);
        semilogx(omega, phase(:)); hold on;
    end
    figure(5); subplot(3,4,n);
    grid on; title(namesT{n}); ylabel('дБ');
    legend(labelsT, 'Location', 'best', 'FontSize', 7);
    figure(6); subplot(3,4,n);
    grid on; title(namesT{n}); ylabel('град');
    legend(labelsT, 'Location', 'best', 'FontSize', 7);
end

fprintf('\nСдвиг ЛАЧХ при изменении K относительно K = %.1f:\n', K0);
for i = 1:length(K_values)
    fprintf('K = %.1f: %+6.2f дБ, ЛФЧХ без изменений\n', ...
        K_values(i), 20*log10(K_values(i)/K0));
end

fprintf('\nУровень ЛАЧХ при omega = 1 рад/с, дБ:\n');
fprintf('Звено');
for i = 1:length(K_values)
    fprintf('    K=%.1f', K_values(i));
end
fprintf('\n');
for n = 1:length(W)
    mag = bode(W{n}, 1);
    fprintf('%5d', n);
    for i = 1:length(K_values)
        fprintf('   %7.2f', 20*log10(mag(1)*K_values(i)/K0));
    end
    fprintf('\n');
end

fprintf('\nЗвено №4: пересечение с уровнем 0 дБ:\n');
for i = 1:length(K_values)
    K = K_values(i);
    if K > 1
        fprintf('K = %.1f: omega = %.3f рад/с\n', K, sqrt(K^2-1)/T0);
    else
        fprintf('K = %.1f: пересечения при omega > 0 нет\n', K);
    end
end

fprintf('\nЗвено №4: сопрягающая частота при изменении T:\n');
for i = 1:length(T_values)
    fprintf('T = %.2f с: omega_c = %.3f рад/с\n', T_values(i), 1/T_values(i));
end
