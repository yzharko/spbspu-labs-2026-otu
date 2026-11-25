function figures = L1()
K = 4;
T = 1.5;
T1 = 1.6;
T2 = 2;
t = (0:0.01:20)';
models = {
    'Интегрирующее', K, [1 0];
    'Идеальное дифференцирующее', [K 0], 1;
    'Усилительное', K, 1;
    'Апериодическое первого порядка', K, [T 1];
    'Апериодическое второго порядка', K, [T2^2 2*T2 1];
    'Колебательное', K, [T2^2 T1 1];
    'Консервативное', K, [T 0 1];
    'Инерционное интегрирующее', K, [T 1 0];
    'Реальное дифференцирующее', [K 0], [T 1];
    'Форсирующее первого порядка', [K*T K], 1;
    'Изодромное', [K*T K], [1 0]
};
reference = regularParts(t,K,T,T1,T2);
maxError = 0;
groups = {1:4, 5:8, 9:11};
labels = {'h(t)', 'w(t)'};
figures = cell(1,5);
for group = 1:numel(groups)
    indices = groups{group};
    figures{group} = figure('Color','w','Position',[50 50 1000 215*numel(indices)]);
    for row = 1:numel(indices)
        i = indices(row);
        for kind = 1:2
            subplot(numel(indices),2,2*(row-1)+kind);
            [y, singular] = response(models{i,2}, models{i,3}, kind, t);
            error = max(abs(y(:)-reference{i,kind}));
            maxError = max(maxError,error);
            plot(t,y,'LineWidth',1.4,'DisplayName','K = 4');
            hold on
            halfGain = response(models{i,2}/2,models{i,3},kind,t);
            plot(t,halfGain,'--','LineWidth',1.2,'DisplayName','K = 2');
            if row == 1 && kind == 1
                legend('Location','northeast');
            end
            if ~isempty(singular)
                text(0.06,0.83,['K = 4: ' singular],'Units','normalized');
            end
            title(sprintf('%d. %s',i,models{i,1}),'FontSize',10);
            if row == numel(indices)
                xlabel('t, с');
            end
            ylabel([labels{kind} ', отн. ед.']);
            grid on
        end
    end
end

fprintf('Максимальная ошибка по 22 характеристикам: %.3e\n',maxError);
assert(maxError < 1e-8,'Проверка аналитических характеристик не пройдена.');

figures{4} = figure('Color','w','Position',[50 50 1000 360]);
subplot(1,2,1);
hold on
for tau = T*[0.2 0.5 1 2 5]
    time = (0:0.01:25)';
    y = step(tf(K,[tau 1]),time);
    plot(time,y,'LineWidth',1.3,'DisplayName',sprintf('T = %g с',tau));
    fprintf('T = %g с, время установления %.2f с\n',tau,tau*log(20));
end
formatPlot('Апериодическое звено');

subplot(1,2,2);
hold on
time = (0:0.001:70)';
for damping = T1*[1 0.5 0.25]
    y = step(tf(K,[T2^2 damping 1]),time);
    outside = find(abs(y-K) > 0.05*K,1,'last');
    peak = 100*(max(y)/K-1);
    fprintf('T1 = %g с, выброс %.1f %%, время установления %.2f с\n', ...
        damping,peak,time(outside+1));
    plot(time(1:10:end),y(1:10:end),'LineWidth',1.2, ...
        'DisplayName',sprintf('T_1 = %g с',damping));
end
formatPlot('Колебательное звено');

figures{5} = figure('Color','w','Position',[50 50 800 370]);
hold on
for gain = K*[0.5 1 2]
    y = step(tf(gain,[T2^2 T1 1]),t);
    plot(t,y,'LineWidth',1.4,'DisplayName',sprintf('K = %g',gain));
    fprintf('K = %g, максимум h(t) = %.3f\n',gain,max(y));
end
formatPlot('Влияние усиления на колебательное звено');
end

function [y, singular] = response(num,den,kind,t)
if kind == 1
    den = [den 0];
end
if numel(num) >= numel(den)
    [q,r] = deconv(num,den);
else
    q = 0;
    r = num;
end
first = find(abs(r) > 1e-12,1);
if isempty(first)
    y = zeros(size(t));
else
    y = impulse(tf(r(first:end),den),t);
end
terms = {};
if q(end) ~= 0
    terms{end+1} = sprintf('%.3gδ',q(end));
end
if numel(q) > 1 && q(end-1) ~= 0
    terms{end+1} = sprintf('%.3gδ′',q(end-1));
end
singular = strjoin(terms,' + ');
end

function formatPlot(name)
title(name);
xlabel('t, с');
ylabel('h(t), отн. ед.');
grid on
legend('Location','northeast');
end

function values = regularParts(t,K,T,T1,T2)
a = T1/(2*T2^2);
b = sqrt(1/T2^2-a^2);
E = exp(-t/T);
F = exp(-t/T2);
D = exp(-a*t);
one = ones(size(t));
zero = zeros(size(t));
values = {
    K*t, K*one;
    zero, zero;
    K*one, zero;
    K*(1-E), K/T*E;
    K*(1-(1+t/T2).*F), K/T2^2*t.*F;
    K*(1-D.*(cos(b*t)+a/b*sin(b*t))), K/(b*T2^2)*D.*sin(b*t);
    K*(1-cos(t/sqrt(T))), K/sqrt(T)*sin(t/sqrt(T));
    K*(t-T*(1-E)), K*(1-E);
    K/T*E, -K/T^2*E;
    K*one, zero;
    K*(T+t), K*one
};
end
