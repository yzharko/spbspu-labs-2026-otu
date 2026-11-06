clear; clc; close all;

K0 = 1.2;
T = 0.3;
T1 = 0.5;
T2 = 0.6;
T1a = 1.5;

K_values = [0.6 1.2 2.4];
time = 0:0.01:10;

syms s_sym t K real
s_tf = tf('s');

names = {'Интегрирующее','Идеальное дифференцирующее', ...
    'Усилительное','Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка','Колебательное', ...
    'Консервативное','Инерционное интегрирующее', ...
    'Реальное дифференцирующее','Форсирующее','Изодромное'};

n = numel(names);

W = {K/s_sym,K*s_sym,K,K/(T*s_sym+1), ...
    K/(T2^2*s_sym^2+T1a*s_sym+1), ...
    K/(T2^2*s_sym^2+T1*s_sym+1), ...
    K/(T*s_sym^2+1),K/(s_sym*(T*s_sym+1)), ...
    K*s_sym/(T*s_sym+1),K*(T*s_sym+1), ...
    K*(T*s_sym+1)/s_sym};

U_step = laplace(heaviside(t),t,s_sym);
U_impulse = laplace(dirac(t),t,s_sym);

R = cell(n,2);
notes = cell(n,2);

notes{2,1} = 'K delta(t)';
notes{10,1} = '+ K*T delta(t)';
notes{2,2} = 'K delta''(t)';
notes{3,2} = 'K delta(t)';
notes{9,2} = '+ (K/T) delta(t)';
notes{10,2} = 'K delta(t) + K*T delta''(t)';
notes{11,2} = '+ K*T delta(t)';

for i = 1:n
    Y_step = simplify(W{i}*U_step);
    Y_impulse = simplify(W{i}*U_impulse);

    h = simplify(ilaplace(Y_step,s_sym,t));
    w = simplify(ilaplace(Y_impulse,s_sym,t));

    R{i,1} = simplify(subs(subs(subs( ...
        h,dirac(1,t),0),dirac(t),0),heaviside(t),1));

    R{i,2} = simplify(subs(subs(subs( ...
        w,dirac(1,t),0),dirac(t),0),heaviside(t),1));

    fprintf(['\n%d. %s\nW(s) = %s\nY_step(s) = %s\n' ...
        'h(t) = %s\nY_impulse(s) = %s\nw(t) = %s\n'], ...
        i,names{i},char(W{i}),char(Y_step), ...
        char(h),char(Y_impulse),char(w));
end

for mode = 1:2
    if mode == 1
        label = 'h(t)';
    else
        label = 'w(t)';
    end

    figure('Name',label,'Color','w', ...
        'Position',[50 50 1450 900]);

    for i = 1:n
        subplot(4,3,i);
        hold on;
        grid on;

        singular = (mode == 1 && i == 2) || ...
            (mode == 2 && ismember(i,[2 3 10]));

        if ~singular
            for k = K_values
                y = double(subs(subs(R{i,mode},K,k),t,time));

                if isscalar(y)
                    y = repmat(y,size(time));
                end

                plot(time,y,'LineWidth',1.2, ...
                    'DisplayName',sprintf('K = %.1f',k));
            end
        end

        if ~isempty(notes{i,mode})
            text(0.08,0.85,notes{i,mode}, ...
                'Units','normalized', ...
                'Interpreter','none','FontSize',8);
        end

        title(sprintf('%d. %s',i,names{i}),'FontSize',9);
        xlabel('t, с');
        ylabel(label);
        xlim([0 10]);

        if i == 1
            legend('show','Location','best');
        end
    end
end

Wtf = {K0/s_tf,K0*s_tf,tf(K0),K0/(T*s_tf+1), ...
    K0/(T2^2*s_tf^2+T1a*s_tf+1), ...
    K0/(T2^2*s_tf^2+T1*s_tf+1), ...
    K0/(T*s_tf^2+1),K0/(s_tf*(T*s_tf+1)), ...
    K0*s_tf/(T*s_tf+1),K0*(T*s_tf+1), ...
    K0*(T*s_tf+1)/s_tf};

fprintf('\nПроверка через Control System Toolbox:\n');

step_check = setdiff(1:n,[2 10]);

for i = step_check
    y = double(subs(subs(R{i,1},K,K0),t,time));

    if isscalar(y)
        y = repmat(y,size(time));
    end

    test = step(Wtf{i},time);

    fprintf('Звено %d: ошибка h = %.3g\n', ...
        i,max(abs(test(:)'-y)));
end

impulse_check = [1 4 5 6 7 8 9 11];

for i = impulse_check
    if i == 9
        model = -(K0/T)/(T*s_tf+1);
    elseif i == 11
        model = K0/s_tf;
    else
        model = Wtf{i};
    end

    y = double(subs(subs(R{i,2},K,K0),t,time));

    if isscalar(y)
        y = repmat(y,size(time));
    end

    test = impulse(model,time);

    fprintf('Звено %d: ошибка w = %.3g\n', ...
        i,max(abs(test(:)'-y)));
end

figure('Name','Влияние T','Color','w');
hold on;
grid on;

for Tv = T*[1/5 1/2 1 2 5]
    model = K0/(Tv*s_tf+1);

    plot(time,step(model,time),'LineWidth',1.3, ...
        'DisplayName',sprintf('T = %.2f',Tv));
end

title('Апериодическое звено 1-го порядка');
xlabel('t, с');
ylabel('h(t)');
legend('show');

figure('Name','Влияние T1','Color','w');
hold on;
grid on;

time_long = 0:0.01:40;

for T1v = [T1 T1/2 T1/5 0]
    model = K0/(T2^2*s_tf^2+T1v*s_tf+1);

    plot(time_long,step(model,time_long),'LineWidth',1.3, ...
        'DisplayName',sprintf('T1 = %.2f',T1v));
end

title('Колебательное звено, T2 = 0.6');
xlabel('t, с');
ylabel('h(t)');
legend('show');
