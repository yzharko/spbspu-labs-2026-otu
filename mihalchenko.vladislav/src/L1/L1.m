clear; clc; close all;

K = 1.1; T = 0.2; T1 = 0.3; T2 = 0.5; T1a = 3*T2;
K_values = K*[0.5 1 2];
T_values = T*[1/5 1/2 1 2 5];
T1_values = T1*[2 1 0.5 0.2 0];
time = (0:0.005:10)'; time_T = (0:0.002:30*T)';
time_osc = (0:0.01:100*T2)'; dt = 0.001*T2;

names = {'Интегрирующее','Дифференцирующее','Усилительное', ...
    'Апериодическое 1-го порядка','Апериодическое 2-го порядка', ...
    'Колебательное','Консервативное','Инерционное интегрирующее', ...
    'Реальное дифференцирующее','Форсирующее','Изодромное'};
n_links = numel(names);
w_labels = {'','w(t) = K\delta''(t)','w(t) = K\delta(t)', ...
    '','','','','','','w(t) = K\delta(t) + KT\delta''(t)',''};
h_sym = cell(1,n_links); w_sym = cell(1,n_links);
first_order_data = zeros(numel(T_values),4);
oscillation_data = zeros(numel(T1_values),6);
checks = NaN(n_links+numel(K_values)+numel(T_values)+numel(T1_values)+1,2);
row = n_links;
visible = 'on';
if ~usejava('desktop'), visible = 'off'; end
out = fullfile(fileparts(mfilename('fullpath')),'results');
if ~isfolder(out), mkdir(out); end
fid = fopen(fullfile(out,'matlab_run.txt'),'w');
assert(fid ~= -1, 'Не удалось открыть файл результатов.');
fclose(fid); diary(fullfile(out,'matlab_run.txt'));
assert(T1a > 2*T2 && all(T1_values >= 0 & T1_values < 2*T2));

a = (T1a+sqrt(T1a^2-4*T2^2))/2; b = T2^2/a;
alpha = T1/(2*T2^2); omega = sqrt(4*T2^2-T1^2)/(2*T2^2);
z = zeros(size(time));
% В h и w записаны регулярные части, без δ(t) и её производной.
h = {K*time,z,K+z,K*(1-exp(-time/T)), ...
    K*(1-(a*exp(-time/a)-b*exp(-time/b))/(a-b)), ...
    K*(1-exp(-alpha*time).*(cos(omega*time)+alpha/omega*sin(omega*time))), ...
    K*(1-cos(time/sqrt(T))),K*(time-T+T*exp(-time/T)), ...
    K/T*exp(-time/T),K+z,K*(T+time)};
w = {K+z,z,z,K/T*exp(-time/T),K*(exp(-time/a)-exp(-time/b))/(a-b), ...
    K/(T2^2*omega)*exp(-alpha*time).*sin(omega*time), ...
    K/sqrt(T)*sin(time/sqrt(T)),K*(1-exp(-time/T)), ...
    -K/T^2*exp(-time/T),z,K+z};

syms s t
W = {K/s,K*s,sym(K),K/(T*s+1),K/(T2^2*s^2+T1a*s+1), ...
    K/(T2^2*s^2+T1*s+1),K/(T*s^2+1),K/(s*(T*s+1)), ...
    K*s/(T*s+1),K*(T*s+1),K*(T*s+1)/s};
U_step = laplace(heaviside(t),t,s);
U_impulse = laplace(dirac(t),t,s);
s_tf = tf('s');
W_tf = {K/s_tf,K*s_tf,tf(K),K/(T*s_tf+1), ...
    K/(T2^2*s_tf^2+T1a*s_tf+1),K/(T2^2*s_tf^2+T1*s_tf+1), ...
    K/(T*s_tf^2+1),K/(s_tf*(T*s_tf+1)), ...
    K*s_tf/(T*s_tf+1),K*(T*s_tf+1),K*(T*s_tf+1)/s_tf};

fprintf('K = %g; T = %g; T1 = %g; T2 = %g; T1a = %g.\n',K,T,T1,T2,T1a);
fprintf('Отклики причинные; dirac(1,t) означает производную dirac(t).\n');
for i = 1:n_links
    Y_step = W{i}*U_step; Y_impulse = W{i}*U_impulse;
    h_sym{i} = simplify(ilaplace(Y_step,s,t));
    w_sym{i} = simplify(ilaplace(Y_impulse,s,t));
    assert(isAlways(simplify(laplace(h_sym{i},t,s)-Y_step) == 0));
    assert(isAlways(simplify(laplace(w_sym{i},t,s)-Y_impulse) == 0));
    fprintf('\n%d. %s\nW = %s\nh = %s\nw = %s\n', ...
        i,names{i},char(W{i}),char(h_sym{i}),char(w_sym{i}));
    % Для звеньев 2 и 10 степень числителя выше степени знаменателя.
    if ~ismember(i,[2 10])
        checks(i,:) = verify(W_tf{i},time,h{i},w{i});
    end
end

fig = figure('Name','Переходные характеристики','Visible',visible, ...
    'Theme','light','Color','w','Position',[50 50 1100 1300]);
for i = 1:n_links
    subplot(4,3,i);
    if i == 2
        axis off;
        text(0.5,0.5,{'h(t) = K\delta(t)','Дельта-функция:', ...
            'обычная кривая не строится.'}, ...
            'Units','normalized','HorizontalAlignment','center','FontSize',11);
    else
        plot(time,h{i},'LineWidth',1.4); grid on;
        xlabel('t, с'); ylabel('h(t)');
        if ismember(i,[4 9]), xlim([0 5*T]); end
        if i == 10
            ylabel('h_{reg}(t)');
            text(0.05,0.85,{'Без импульсной составляющей', ...
                '+ KT\delta(t) при t = 0'},'Units','normalized','FontSize',10);
        end
    end
    set(gca,'FontSize',12);
    title(sprintf('%d. %s',i,names{i}),'FontSize',11);
end
sgtitle(sprintf('Переходные характеристики, K = %g',K));

drawnow;
exportgraphics(fig,fullfile(out,'step_all.png'),'Resolution',200);

fig = figure('Name','Импульсные характеристики','Visible',visible, ...
    'Theme','light','Color','w','Position',[50 50 1100 1300]);
for i = 1:n_links
    subplot(4,3,i);
    if ismember(i,[2 3 10])
        axis off;
        reason = 'Дельта-функция:';
        if i ~= 3, reason = 'Производная дельта-функции:'; end
        text(0.5,0.5,{w_labels{i},reason,'обычная кривая не строится.'}, ...
            'Units','normalized','HorizontalAlignment','center','FontSize',11);
    else
        plot(time,w{i},'LineWidth',1.4); grid on;
        xlabel('t, с'); ylabel('w(t)');
        if ismember(i,[4 9]), xlim([0 5*T]); end
        if ismember(i,[9 11])
            ylabel('w_{reg}(t)');
            label = '+ KT\delta(t)';
            if i == 9, label = '+ (K/T)\delta(t)'; end
            text(0.05,0.15,{'Без импульсной составляющей', ...
                [label ' при t = 0']},'Units','normalized','FontSize',10);
        end
    end
    set(gca,'FontSize',12);
    title(sprintf('%d. %s',i,names{i}),'FontSize',11);
end
sgtitle(sprintf('Импульсные характеристики, K = %g',K));

drawnow;
exportgraphics(fig,fullfile(out,'impulse_all.png'),'Resolution',200);

fig = figure('Name','Влияние K','Visible',visible,...
    'Theme','light','Color','w','Position',[100 100 1100 400]);
for k = K_values
    h_k = h{6}*k/K; w_k = w{6}*k/K;
    row = row+1;
    checks(row,:) = verify(k/(T2^2*s_tf^2+T1*s_tf+1),time,h_k,w_k);
    subplot(1,2,1); hold on; grid on;
    plot(time,h_k,'DisplayName',sprintf('K = %g',k));
    xlabel('t, с'); ylabel('h(t)'); legend('show','Location','best');
    subplot(1,2,2); hold on; grid on;
    plot(time,w_k,'DisplayName',sprintf('K = %g',k));
    xlabel('t, с'); ylabel('w(t)'); legend('show','Location','best');
end

set(findall(fig,'Type','axes'),'FontSize',14);
drawnow;
exportgraphics(fig,fullfile(out,'gain.png'),'Resolution',200);

fig = figure('Name','Влияние T','Visible',visible,...
    'Theme','light','Color','w','Position',[100 100 1100 400]);
j = 0;
for Ti = T_values
    j = j+1; row = row+1;
    h_T = K*(1-exp(-time_T/Ti)); w_T = K/Ti*exp(-time_T/Ti);
    checks(row,:) = verify(K/(Ti*s_tf+1),time_T,h_T,w_T);
    first_order_data(j,:) = [Ti Ti log(9)*Ti -log(0.02)*Ti];
    subplot(1,2,1); hold on; grid on;
    plot(time_T,h_T,'DisplayName',sprintf('T = %g',Ti));
    xlabel('t, с'); ylabel('h(t)'); legend('show','Location','best');
    subplot(1,2,2); hold on; grid on;
    plot(time_T,w_T,'DisplayName',sprintf('T = %g',Ti));
    xlabel('t, с'); ylabel('w(t)'); legend('show','Location','best');
end

set(findall(fig,'Type','axes'),'FontSize',14);
drawnow;
exportgraphics(fig,fullfile(out,'first_order.png'),'Resolution',200);

fig = figure('Name','Влияние T1','Visible',visible,...
    'Theme','light','Color','w','Position',[100 100 1100 400]);
j = 0;
for T1i = T1_values
    j = j+1; row = row+1;
    zeta = T1i/(2*T2); alpha = T1i/(2*T2^2); omega = sqrt(1-zeta^2)/T2;
    h_T1 = K*(1-exp(-alpha*time_osc).* ...
        (cos(omega*time_osc)+alpha/omega*sin(omega*time_osc)));
    w_T1 = K/(T2^2*omega)*exp(-alpha*time_osc).*sin(omega*time_osc);
    checks(row,:) = verify(K/(T2^2*s_tf^2+T1i*s_tf+1),time_osc,h_T1,w_T1);
    sigma = 100*exp(-pi*zeta/sqrt(1-zeta^2)); tau = Inf; ts = Inf;
    if T1i > 0
        tau = 1/alpha;
        % После bound огибающая целиком остаётся в полосе ±2 % от K.
        bound = log(1/(0.02*sqrt(1-zeta^2)))/alpha;
        time_fine = (0:dt:bound+10*T2)';
        h_fine = K*(1-exp(-alpha*time_fine).* ...
            (cos(omega*time_fine)+alpha/omega*sin(omega*time_fine)));
        last = find(abs(h_fine/K-1) > 0.02,1,'last');
        ts = time_fine(last+1);
    end
    oscillation_data(j,:) = [T1i zeta sigma pi/omega tau ts];
    subplot(1,2,1); hold on; grid on;
    plot(time_osc,h_T1,'DisplayName',sprintf('T_1 = %g',T1i));
    xlabel('t, с'); ylabel('h(t)'); xlim([0 20*T2]);
    legend('show','Location','best');
    subplot(1,2,2); hold on; grid on;
    plot(time_osc,w_T1,'DisplayName',sprintf('T_1 = %g',T1i));
    xlabel('t, с'); ylabel('w(t)'); legend('show','Location','best');
end

set(findall(fig,'Type','axes'),'FontSize',14);
drawnow;
exportgraphics(fig,fullfile(out,'oscillatory.png'),'Resolution',200);

fig = figure('Name','Проверка результатов','Visible',visible,...
    'Theme','light','Color','w','Position',[100 100 1100 400]);
subplot(1,2,1);
plot(time,h{4},time,step(W_tf{4},time),'--','LineWidth',1.3); grid on;
xlabel('t, с'); ylabel('h(t)'); xlim([0 5*T]);
legend('Аналитически','step','Location','best');
subplot(1,2,2);
plot(time,w{4},time,impulse(W_tf{4},time),'--','LineWidth',1.3); grid on;
xlabel('t, с'); ylabel('w(t)'); xlim([0 5*T]);
legend('Аналитически','impulse','Location','best');
row = row+1;
checks(row,:) = verify(K/(T2^2*s_tf^2+2*T2*s_tf+1),time, ...
    K*(1-(1+time/T2).*exp(-time/T2)),K*time/T2^2.*exp(-time/T2));

set(findall(fig,'Type','axes'),'FontSize',14);
drawnow;
exportgraphics(fig,fullfile(out,'verification.png'),'Resolution',200);
writetable(array2table(checks,'VariableNames',{'StepError','ImpulseRegularError'}), ...
    fullfile(out,'verification.csv'));
writetable(array2table(first_order_data,'VariableNames', ...
    {'T','Time63','Rise10_90','Settling2'}), ...
    fullfile(out,'first_order.csv'));
writetable(array2table(oscillation_data,'VariableNames', ...
    {'T1','Zeta','Overshoot','PeakTime','EnvelopeTau','Settling2'}), ...
    fullfile(out,'oscillatory.csv'));
save(fullfile(out,'results.mat'),'K','T','T1','T2','T1a','h_sym','w_sym','h','w', ...
    'K_values','T_values','T1_values','checks','first_order_data','oscillation_data');
fprintf('Численных сравнений: %d. Все проверки пройдены.\n',sum(isfinite(checks),'all'));
fprintf('Максимальная ошибка h: %.3e; w: %.3e.\n', ...
    max(checks(:,1),[],'omitnan'),max(checks(:,2),[],'omitnan'));
diary off;

function errors = verify(W,time,h,w)
    h_check = step(W,time); w_check = impulse(W,time);
    errors = [max(abs(h_check-h)) max(abs(w_check-w))];
    assert(all(errors < 1e-8), 'Численная проверка не пройдена.');
end
