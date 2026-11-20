syms s t K
s_tf = tf('s');
K0 = 1.0; T = 0.2; T1 = 0.4; T2 = 0.5;
T1a = 1.5;   % T1 для апериодического звена 2-го порядка
time = 0:0.01:10;

names = ["Интегрирующее", ...
         "Дифференцирующее идеальное", ...
         "Усилительное", ...
         "Апериодическое 1-го порядка", ...
         "Апериодическое 2-го порядка", ...
         "Колебательное", ...
         "Консервативное", ...
         "Инерционное интегрирующее", ...
         "Реальное дифференцирующее", ...
         "Форсирующее 1-го порядка", ...
         "Изодромное"];

% символьные — для ilaplace
W_sym_dict = dictionary(names, ...
    [K/s, K*s, K, K/(T*s+1), ...
     K/(T2^2*s^2+T1a*s+1), K/(T2^2*s^2+T1*s+1), ...
     K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), ...
     K*(T*s+1), K*(T*s+1)/s]);

% числовые — сразу с K0, для step/impulse
W_tf_dict = dictionary(names, ...
    [K0/s_tf, K0*s_tf, tf(K0,1), K0/(T*s_tf+1), ...
     K0/(T2^2*s_tf^2+T1a*s_tf+1), K0/(T2^2*s_tf^2+T1*s_tf+1), ...
     K0/(T*s_tf^2+1), K0/(s_tf*(T*s_tf+1)), K0*s_tf/(T*s_tf+1), ...
     K0*(T*s_tf+1), K0*(T*s_tf+1)/s_tf]);

key_list = keys(W_sym_dict);

figure(1); hold on; grid on;
figure(2); hold on; grid on;
figure(3); hold on; grid on;
figure(4); hold on; grid on;

U_step    = laplace(1, t, s);         
U_imp = laplace(dirac(t), t, s);  

for i = 1:numel(key_list)
    name  = key_list{i};
    W     = W_sym_dict(name);
    W_check = W_tf_dict(name);

    Y_step    = W * U_step;
    h = ilaplace(Y_step, s, t);
    Y_imp = W * U_imp;
    w = ilaplace(Y_imp, s, t);

    h_k = subs(h, K, K0);
    figure(1); plot(time, double(subs(h_k, t, time)), 'DisplayName', name);
    w_k = subs(w, K, K0);
    figure(2); plot(time, double(subs(w_k, t, time)), 'DisplayName', name);

    try
        figure(3); step(W_check, time);
    catch
        warning('step пропущен: %s', name);
    end
    try
        figure(4); impulse(W_check, time);
    catch
        warning('impulse пропущен: %s', name);
    end
end

figure(1); title('h(t), K=K0'); xlabel('t'); ylabel('h(t)'); legend('show');
figure(2); title('w(t), K=K0'); xlabel('t'); ylabel('w(t)'); legend('show');
figure(3); title('step проверка'); xlabel('t'); ylabel('выход');
figure(4); title('impulse проверка'); xlabel('t'); ylabel('выход');

%% Пункт *:
% Влияние K на интегрирующее звено
K_values = 1:10;

% --- h(t) для разных K ---
figure(5); hold on; grid on;
W_int = W_sym_dict("Интегрирующее");     
h_int = ilaplace(W_int * U_step, s, t); 

for j = 1:numel(K_values)
    Kj = K_values(j);
    y = double(subs(subs(h_int, K, Kj), t, time));
    plot(time, y, 'DisplayName', sprintf('K=%g', Kj));
end
title('Влияние K на h(t) интегрирующего звена');
xlabel('Время t, с'); ylabel('h(t)');
legend('show','Location','northwest');

% --- w(t) для разных K ---
figure(6); hold on; grid on;
w_int = ilaplace(W_int * U_imp, s, t);

for j = 1:numel(K_values)
    Kj = K_values(j);
    y = double(subs(subs(w_int, K, Kj), t, time));
    plot(time, y, 'DisplayName', sprintf('K=%g', Kj));
end
title('Влияние K на w(t) интегрирующего звена');
xlabel('Время t, с'); ylabel('w(t)');
legend('show','Location','northwest');

% Влияние T на апериодическое звено №4
T_nom = 0.2;
T_list = [T_nom, T_nom*2, T_nom*5, T_nom/2, T_nom/5];
labels = {'T=0.2 (номинал)', 'T=0.4 (×2)', 'T=1.0 (×5)', ...
          'T=0.1 (÷2)', 'T=0.04 (÷5)'};

figure(7); hold on; grid on;
for j = 1:numel(T_list)
    Tj = T_list(j);
    W_T = K0 / (Tj*s_tf + 1);
    step(W_T, time);
    L = findobj(gca, 'Type', 'line');
    if ~isempty(L), L(1).DisplayName = labels{j}; end
end
title('Влияние T на h(t)'); xlabel('t'); ylabel('h(t)');
legend('show');

% Влияние T1 на колебательное звено №6 (T2 фиксирована)
T2_fix = 0.5;
T1_list = [0.4, 0.3, 0.2, 0.1, 0.05];
labels6 = {'T1=0.4', 'T1=0.3', 'T1=0.2', 'T1=0.1', 'T1=0.05'};

figure(8); hold on; grid on;
for j = 1:numel(T1_list)
    T1j = T1_list(j);
    W_T1 = K0 / (T2_fix^2*s_tf^2 + T1j*s_tf + 1);
    step(W_T1, time);
    L = findobj(gca, 'Type', 'line');
    if ~isempty(L), L(1).DisplayName = labels6{j}; end
end
title('Влияние T1 на h(t) при T2=0.5'); xlabel('t'); ylabel('h(t)');
legend('show');
