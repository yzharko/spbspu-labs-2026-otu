clearvars;
K = 1.6;
T = 0.8;
T1 = 1.0;
T2 = 1.1;
T1a = 2*T2;                    % Граница апериодичности только для звена 5
gainValues = 1:10;
timeValues = T*[0.2 0.5 1 2 5];
dampingValues = T1*[1 0.75 0.5 0.25 0.1];
time = (0:0.01:20)';
dpi = 300;
out = fullfile(fileparts(mfilename('fullpath')), 'graphs');
if ~isfolder(out)
    mkdir(out);
end
names = {'Интегрирующее','Дифференцирующее','Усилительное',...
    'Апериодическое I порядка','Апериодическое II порядка','Колебательное',...
    'Консервативное','Инерционное интегрирующее','Реальное дифференцирующее',...
    'Форсирующее','Изодромное'};
syms s t
k = sym(K,'r'); tau = sym(T,'r'); a = sym(T1,'r'); b = sym(T2,'r');
W = {k/s,k*s,k,k/(tau*s+1),k/(b^2*s^2+sym(T1a,'r')*s+1),...
    k/(b^2*s^2+a*s+1),k/(tau*s^2+1),k/(s*(tau*s+1)),...
    k*s/(tau*s+1),k*(tau*s+1),k*(tau*s+1)/s};
% Полиномиальная часть образа даёт дельта-функцию и её производные.
% Для строго правильного остатка impulse даёт регулярную часть как h, так и w.
labels = {'h','w'};
for mode = 1:2
    for group = 1:2
        indices = (1:6) + 6*(group-1);
        indices = indices(indices<=numel(W));
        fig = figure('Color','w','Position',[100 100 980 820]);
        tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
        for i = indices
            F = simplify(W{i}/s^(2-mode));
            [num,den] = numden(F);
            [poly,remn] = quorem(num,den,s);
            proper = simplify(remn/den);
            regular = simplify(ilaplace(proper,s,t));
            d0 = subs(poly,s,0);
            d1 = diff(poly,s);
            assert(isequal(simplify(poly-d0-s*d1),sym(0)));
            answer = regular + d0*dirac(t) + d1*dirac(1,t);
            fprintf('%d. %s: W(s) = %s, %s(t) = %s\n',...
                i,names{i},char(W{i}),labels{mode},char(answer));
            if isequal(proper,sym(0))
                y = zeros(size(time));
            else
                [n,d] = numden(proper);
                y = squeeze(impulse(tf(sym2poly(n),sym2poly(d)),time));
            end
            nexttile;
            plot(time,y,'LineWidth',1.5); grid on;
            title(sprintf('%d. %s',i,names{i}));
            xlabel('t, с'); ylabel([labels{mode} '(t)']);
            if ~isequal(poly,sym(0))
                terms = {};
                if d0 ~= 0
                    terms{end+1} = sprintf('%g\\delta(t)',double(d0));
                end
                if d1 ~= 0
                    terms{end+1} = sprintf('%g\\delta''(t)',double(d1));
                end
                text(0.97,0.89,['При t=0: ' strjoin(terms,' + ')],...
                    'Units','normalized','HorizontalAlignment','right');
            end
        end
        print(fig,fullfile(out,sprintf('%s%d.png',labels{mode},group)),...
            '-dpng',sprintf('-r%d',dpi));
    end
end
%% Исследования параметров
p = tf('s');
sets = {gainValues,timeValues,dampingValues};
tags = {'K','T','T1'};
heads = {'Колебательное звено: влияние K','Апериодическое звено: влияние T',...
    'Колебательное звено: влияние T_1'};
for j = 1:3
    fig = figure('Color','w','Position',[100 100 900 510]);
    hold on; grid on;
    for v = sets{j}
        switch j
            case 1
                model = v/(T2^2*p^2+T1*p+1);
            case 2
                model = K/(v*p+1);
                fprintf('T=%g: t95=%.4f с\n',v,v*log(20));
            case 3
                model = K/(T2^2*p^2+v*p+1);
                alpha = v/(2*T2^2);
                omega = sqrt(1/T2^2-alpha^2);
                peakRatio = exp(-pi*alpha/omega);
                n = ceil(log(0.05)/log(peakRatio))-1;
                left = n*pi/omega;
                right = left+(pi/2+atan(alpha/omega))/omega;
                err = @(x) abs(exp(-alpha*x).*(cos(omega*x)+...
                    alpha/omega*sin(omega*x)))-0.05;
                ts = fzero(err,[left right]);
                fprintf('T1=%g: максимум=%.2f %%, t5=%.3f с\n',v,100*peakRatio,ts);
        end
        plot(time,squeeze(step(model,time)),'LineWidth',1.3,...
            'DisplayName',sprintf('%s = %g',tags{j},v));
    end
    xlabel('t, с'); ylabel('h(t)'); title(heads{j});
    legend('Location','eastoutside');
    print(fig,fullfile(out,[tags{j} '.png']),'-dpng',sprintf('-r%d',dpi));
end