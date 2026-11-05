clear; clc; close all;

K = 0.8;
T = 0.5;
T1 = 0.7;
T2 = 0.5;
T1a = 1.25;

gains = K*[0.5 1 2];
tt = linspace(0,10,1001);

syms g s t real
st = tf('s');

titles = {'Интегрирующее','Идеальное дифференцирующее','Усилительное', ...
    'Апериодическое 1-го порядка','Апериодическое 2-го порядка','Колебательное', ...
    'Консервативное','Инерционное интегрирующее','Реальное дифференцирующее', ...
    'Форсирующее','Изодромное'};

W = {@(g,s) g/s, ...
     @(g,s) g*s, ...
     @(g,s) g + 0*s, ...
     @(g,s) g/(T*s+1), ...
     @(g,s) g/(T2^2*s^2+T1a*s+1), ...
     @(g,s) g/(T2^2*s^2+T1*s+1), ...
     @(g,s) g/(T*s^2+1), ...
     @(g,s) g/(s*(T*s+1)), ...
     @(g,s) g*s/(T*s+1), ...
     @(g,s) g*(T*s+1), ...
     @(g,s) g*(T*s+1)/s};

n = numel(W);
hf = cell(n,1);
wf = cell(n,1);
c0 = zeros(n,1);
c1 = zeros(n,1);

for i = 1:n
    Ws = W{i}(g,s);
    a1 = limit(Ws/s,s,inf);
    a0 = limit(Ws - a1*s,s,inf);

    hr = ilaplace((Ws - a1*s)/s,s,t);
    wr = ilaplace(Ws - a0 - a1*s,s,t);

    hf{i} = matlabFunction(hr,'Vars',{g,t});
    wf{i} = matlabFunction(wr,'Vars',{g,t});
    c0(i) = double(subs(a0,g,K));
    c1(i) = double(subs(a1,g,K));

    fprintf('%2d. %s\n', i, titles{i});
    fprintf('    h(t) = %s\n', char(subs(hr + a1*dirac(t),g,K)));
    fprintf('    w(t) = %s\n', char(subs(wr + a0*dirac(t) + a1*diff(dirac(t),t),g,K)));
end

lab = {'h(t)','w(t)'};
fun = {hf,wf};

for m = 1:2
    figure('Color','w','Position',[50 50 900 950]);

    for i = 1:n
        subplot(4,3,i);
        hold on;
        grid on;

        for k = gains
            plot(tt,fun{m}{i}(k,tt) + zeros(size(tt)),'LineWidth',1.2);
        end

        parts = {};
        if m == 1
            if c1(i) ~= 0
                parts{end+1} = sprintf('%g\\delta(t)',c1(i));
            end
        else
            if c0(i) ~= 0
                parts{end+1} = sprintf('%g\\delta(t)',c0(i));
            end
            if c1(i) ~= 0
                parts{end+1} = sprintf('%g\\delta''(t)',c1(i));
            end
        end
        if ~isempty(parts)
            text(0.3,0.85,['+ ' strjoin(parts,' + ')],'Units','normalized','FontSize',9);
        end

        if i == 1
            legend(cellstr(num2str(gains','K = %.1f')),'Location','best');
        end

        title(sprintf('%d. %s',i,titles{i}),'FontSize',9);
        xlabel('t, с');
        ylabel(lab{m});
        xlim([0 10]);
    end
end

err = zeros(n,2);

for i = 1:n
    G = W{i}(K,st);
    Gh = G - c1(i)*st;
    Gw = G - c0(i) - c1(i)*st;

    err(i,1) = max(abs(step(Gh,tt)' - (hf{i}(K,tt) + zeros(size(tt)))));
    err(i,2) = max(abs(impulse(Gw,tt)' - (wf{i}(K,tt) + zeros(size(tt)))));
end

fprintf('\nМаксимальная ошибка (аналитика и step/impulse), K = %.1f\n',K);
fprintf('звено     h(t)        w(t)\n');
fprintf('%4d   %9.2e   %9.2e\n',[1:n; err']);

Tk = T*[1/5 1/2 1 2 5];
t15 = linspace(0,15,1501);

figure('Color','w');
hold on;
grid on;

fprintf('\nЗвено 4, K = %.1f\n',K);
fprintf('  T       t95      tуст\n');

for Tj = Tk
    plot(t15,step(K/(Tj*st+1),t15),'LineWidth',1.3);
    fprintf('%5.2f   %6.3f   %6.3f\n',Tj,-Tj*log(0.05),-Tj*log(0.02));
end

legend(cellstr(num2str(Tk','T = %.2f')),'Location','southeast');
title('Апериодическое звено 1-го порядка');
xlabel('t, с');
ylabel('h(t)');

Tz = T1*[1 1/2 1/5 0];
t30 = 0:0.001:30;

figure('Color','w');
hold on;
grid on;

fprintf('\nЗвено 6, T2 = %.1f\n',T2);
fprintf('  T1    sigma, %%    tmax     tуст\n');

for Tj = Tz
    y = step(K/(T2^2*st^2+Tj*st+1),t30);
    plot(t30,y,'LineWidth',1.3);

    if Tj > 0
        r = stepinfo(y,t30,K);
        fprintf('%5.2f   %6.2f   %6.3f   %6.2f\n',Tj,r.Overshoot,r.PeakTime,r.SettlingTime);
    else
        fprintf('%5.2f   %6.2f   %6.3f      inf\n',Tj,100,pi*T2);
    end
end

legend(cellstr(num2str(Tz','T1 = %.2f')),'Location','north','Orientation','horizontal');
title('Колебательное звено, T2 = 0.5');
xlabel('t, с');
ylabel('h(t)');
ylim([0 2*K+0.4]);
