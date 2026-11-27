K = 2.2;
T = 0.2;
T1 = 0.4;
T2 = 0.3;
T1ap = 0.8;
s = tf('s');

W1 = K/s;
W2 = K*s;
W3 = tf(K);
W4 = K/(T*s + 1);
W5 = K/(T2^2*s^2 + T1ap*s + 1);
W6 = K/(T2^2*s^2 + T1*s + 1);
W7 = K/(T*s^2 + 1);
W8 = K/(s*(T*s + 1));
W9 = K*s/(T*s + 1);
W10 = K*(T*s + 1);
W11 = K*(T*s + 1)/s;

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};
headings = {'АФЧХ', 'ЛАЧХ и ЛФЧХ', ...
    'Сравнение K', 'Сравнение K и постоянных времени'};
windowStyle = 'normal';

minFreq = 0.01;
maxFreq = 1000;
w = {minFreq, maxFreq};

W = {W1, W2, W3, W4, W5, W6, W7, W8, W9, W10, W11};

pole = 1/sqrt(T);
gap = 0.001;
beforePole = {minFreq, pole - gap};
afterPole = {pole + gap, maxFreq};
ranges7 = {beforePole, afterPole};
parts7 = {'до полюса', 'после полюса'};

for n = [1:6 8:11]
    model = W{n};
    name = names{n};

    figure('Name', [headings{1} ' - ' name], 'WindowStyle', windowStyle)
    if n == 3
        nyquist(model, 'o', w)
    else
        nyquist(model, w)
    end
    grid on

    figure('Name', [headings{2} ' - ' name], 'WindowStyle', windowStyle)
    bode(model, w)
    grid on
end

for part = 1:2
    name = [names{7} ' ' parts7{part}];
    frequencies = ranges7{part};

    figure('Name', [headings{1} ' - ' name], 'WindowStyle', windowStyle)
    nyquist(W7, frequencies)
    grid on

    figure('Name', [headings{2} ' - ' name], 'WindowStyle', windowStyle)
    bode(W7, frequencies)
    grid on
end

Knew = 4.4;
Tnew = 0.4;
gain = Knew/K;
labelK = ['K = ' num2str(Knew)];
labelT = ['T = ' num2str(Tnew)];

W4new = K/(Tnew*s + 1);
W5new = K/(T2^2*s^2 + 1.2*s + 1);
W6new = K/(T2^2*s^2 + 0.55*s + 1);
W7new = K/(Tnew*s^2 + 1);
W8new = K/(s*(Tnew*s + 1));
W9new = K*s/(Tnew*s + 1);
W10new = K*(Tnew*s + 1);
W11new = K*(Tnew*s + 1)/s;
W5t2 = K/(0.36^2*s^2 + T1ap*s + 1);
W6t2 = K/(0.45^2*s^2 + T1*s + 1);

Wnew = {[], [], [], W4new, W5new, W6new, W7new, W8new, W9new, W10new, W11new};

newPole = 1/sqrt(Tnew);
beforeBoth = {minFreq, newPole - gap};
betweenPoles = {newPole + gap, pole - gap};
afterBoth = {pole + gap, maxFreq};
ranges7new = {beforeBoth, betweenPoles, afterBoth};
parts7new = {'до обоих полюсов', 'между полюсами', 'после обоих полюсов'};

for n = 1:3
    model = W{n};
    figure('Name', [headings{3} ' - ' names{n}], 'WindowStyle', windowStyle)
    bode(model, gain*model, w)
    grid on
    legend('Исходное звено', labelK)
end

for n = [4 8:11]
    model = W{n};
    figure('Name', [headings{4} ' - ' names{n}], 'WindowStyle', windowStyle)
    bode(model, gain*model, Wnew{n}, w)
    grid on
    legend('Исходное звено', labelK, labelT)
end

figure('Name', [headings{4} ' - ' names{5}], 'WindowStyle', windowStyle)
bode(W5, gain*W5, W5new, W5t2, w)
grid on
legend('Исходное звено', labelK, 'T1 = 1.2', 'T2 = 0.36')

figure('Name', [headings{4} ' - ' names{6}], 'WindowStyle', windowStyle)
bode(W6, gain*W6, W6new, W6t2, w)
grid on
legend('Исходное звено', labelK, 'T1 = 0.55', 'T2 = 0.45')

for part = 1:3
    name = [names{7} ' ' parts7new{part}];
    frequencies = ranges7new{part};

    figure('Name', [headings{4} ' - ' name], 'WindowStyle', windowStyle)
    bode(W7, gain*W7, W7new, frequencies)
    grid on
    legend('Исходное звено', labelK, labelT)
end
