clear; clc; close all;

s = tf('s');
Go = 5/(s*(s+0.1));
sigma_req = 10;
t_p_req = 8;

G0 = feedback(Go,1);
S0 = stepinfo(G0);
p = pole(G0);
fprintf('Исходная система: полюса %.3f%+.3fi, %.3f%+.3fi | устойчива: %d\n', ...
        real(p(1)), imag(p(1)), real(p(2)), imag(p(2)), isstable(G0));
fprintf('t_п = %.2f c, sigma = %.2f %%, ошибка = %.4f\n\n', ...
        S0.SettlingTime, S0.Overshoot, dcgain(feedback(1,Go)));

tau = 0.05;
Go_d = Go;  Go_d.InputDelay = tau;
warning('off','all');
[Kcr, ~, Wcg] = margin(Go_d);
warning('on','all');
Tcr = 2*pi/Wcg;
Kp = 0.6*Kcr;  Ti = 0.5*Tcr;  Td = 0.125*Tcr;
fprintf('Зиглер-Николс: K_кр = %.4f, T_кр = %.4f c | Kp = %.4f, Ti = %.4f c, Td = %.4f c\n\n', ...
        Kcr, Tcr, Kp, Ti, Td);

best = [];
for wc_i = 0.5:0.05:3
    for pm_i = 60:2.5:90
        try
            Ct = pidtune(Go,'PID',wc_i,pidtuneOptions('PhaseMargin',pm_i));
        catch
            continue;
        end
        Gt = feedback(Ct*Go,1);
        if ~isstable(Gt), continue; end
        St = stepinfo(Gt);
        if isnan(St.SettlingTime) || isinf(St.SettlingTime), continue; end
        if St.Overshoot <= sigma_req && St.SettlingTime <= t_p_req && ...
           (isempty(best) || (wc_i == best.wc && St.SettlingTime < best.tp))
            best = struct('C',Ct,'wc',wc_i,'pm',pm_i,'tp',St.SettlingTime);
        end
    end
end

C = {pid(Kp, Kp/Ti, Kp*Td), ...
        pidtune(Go,'PID'), ...
        best.C};
name = {'ПИД, Зиглер-Николс','pidtune (по умолчанию)', ...
        sprintf('pidtune (wc = %.2f, PM = %g)', best.wc, best.pm)};

t = 0:0.01:40;
y = cell(1,3);
for k = 1:3
    L = C{k}*Go;
    T = feedback(L,1);
    S = stepinfo(T);
    [Gm, Pm, ~, Wp] = margin(L);
    Cf = pid(C{k}.Kp, C{k}.Ki, C{k}.Kd, C{k}.Kd/C{k}.Kp/10);
    u = max(abs(step(feedback(Cf,Go), 0:0.001:20)));
    fprintf('%s\n  Kp = %.4f, Ki = %.4f, Kd = %.4f\n', name{k}, C{k}.Kp, C{k}.Ki, C{k}.Kd);
    fprintf('  t_п = %.2f c, sigma = %.2f %% | требования выполнены: %d\n', ...
            S.SettlingTime, S.Overshoot, S.SettlingTime <= t_p_req && S.Overshoot <= sigma_req);
    fprintf('  ошибка = %.4f, GM = %.2f дБ, PM = %.2f град (wc = %.3f), пик u = %.3f\n\n', ...
            dcgain(feedback(1,L)), 20*log10(Gm), Pm, Wp, u);
    y{k} = step(T,t);
    figure('Name',['ЛАЧХ/ЛФЧХ: ' name{k}]);  margin(L);  grid on;
end

figure('Name','Переходные характеристики'); hold on; grid on;
plot(t, step(G0,t), 'w--', t, y{1}, 'r', t, y{2}, 'b', t, y{3}, 'g', 'LineWidth', 1.5);
yline(1,'w:');  yline(1+sigma_req/100,'m--');  xline(t_p_req,'m--');
xlabel('t, c');  ylabel('y(t)');  title('Переходные процессы замкнутой системы');
legend('Без регулятора', name{:}, '1', 'Граница \sigma_{доп}', 'Граница t_{п}', 'Location','best');

figure('Name','Крупный план'); hold on; grid on;
plot(t, y{1}, 'r', t, y{3}, 'g', 'LineWidth', 1.5);
yline(1,'w:');  yline(1+sigma_req/100,'m--');  xline(t_p_req,'m--');
xlim([0 20]);  ylim([0 1.8]);
xlabel('t, c');  ylabel('y(t)');  title('Зиглер-Николс и подобранный pidtune');
legend(name{1}, name{3}, '1', 'Граница \sigma_{доп}', 'Граница t_{п}', 'Location','best');
