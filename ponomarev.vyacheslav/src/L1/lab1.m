
syms s t K;

%% Параметры варианта

variantGain = sym('1.9');
timeConstant = sym('0.7');
firstCoefficient = sym('1.1');
secondCoefficient = sym('1.3');

aperiodicFirstCoefficient = sym('3');

time = 0:0.01:20;
impulseTime = 0:0.001:2;
epsilon = 0.05;

%% Входные воздействия

stepInput = laplace(heaviside(t), t, s);
impulseInput = laplace(dirac(t), t, s);

%% Передаточные функции

linkNames = {
    'Интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1-го порядка';
    'Апериодическое 2-го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1-го порядка';
    'Изодромное'
};

transferFunctions = {
    K/s;
    K*s;
    K;
    K/(timeConstant*s + 1);
    K/(secondCoefficient^2*s^2 + aperiodicFirstCoefficient*s + 1);
    K/(secondCoefficient^2*s^2 + firstCoefficient*s + 1);
    K/(timeConstant*s^2 + 1);
    K/(s*(timeConstant*s + 1));
    K*s/(timeConstant*s + 1);
    K*(timeConstant*s + 1);
    K*(timeConstant*s + 1)/s
};

filterOrders = [0, 2, 1, 0, 0, 0, 0, 0, 1, 2, 1];

%% Аналитические вычисления

for linkIndex = 1:numel(transferFunctions)

    transferFunction = subs(transferFunctions{linkIndex}, K, variantGain);

    stepOutput = simplify(transferFunction * stepInput);
    impulseOutput = simplify(transferFunction * impulseInput);

    stepResponse = simplify(ilaplace(stepOutput, s, t));
    impulseResponse = simplify(ilaplace(impulseOutput, s, t));

    fprintf('\n%d. %s\n', linkIndex, linkNames{linkIndex});

    fprintf('W(s) = ');
    disp(transferFunction);

    fprintf('h(t) = ');
    disp(stepResponse);

    fprintf('w(t) = ');
    disp(impulseResponse);

    %% Определение импульсных составляющих

    [stepNumerator, stepDenominator] = numden(stepOutput);
    [stepPolynomial, ~] = quorem(stepNumerator, stepDenominator, s);

    [impulseNumerator, impulseDenominator] = numden(impulseOutput);
    [impulsePolynomial, ~] = quorem(impulseNumerator, impulseDenominator, s);

    hasStepImpulse = ~isequal(stepPolynomial, sym(0));
    hasImpulseImpulse = ~isequal(impulsePolynomial, sym(0));

    %% Приближённая модель

    stepPlotTime = time;
    impulsePlotTime = time;

    if hasStepImpulse || hasImpulseImpulse

        [numerator, denominator] = numden(transferFunction);

        originalModel = tf(sym2poly(numerator), sym2poly(denominator));

        filterModel = tf([epsilon 1], 1);

        approximateModel = originalModel / filterModel^filterOrders(linkIndex);

    end

    %% Переходная характеристика

    if hasStepImpulse

        stepPlotTime = impulseTime;
        stepValues = squeeze(step(approximateModel, stepPlotTime));

    else

        regularStep = simplify(ilaplace(stepOutput, s, t));
        stepFunction = matlabFunction(regularStep, 'Vars', t);

        stepValues = real(double(stepFunction(stepPlotTime)));

        if isscalar(stepValues)
            stepValues = stepValues * ones(size(stepPlotTime));
        end

    end

    %% Импульсная характеристика

    if hasImpulseImpulse

        impulsePlotTime = impulseTime;
        impulseValues = squeeze(impulse(approximateModel, impulsePlotTime));

    else

        regularImpulse = simplify(ilaplace(impulseOutput, s, t));
        impulseFunction = matlabFunction(regularImpulse, 'Vars', t);

        impulseValues = real(double(impulseFunction(impulsePlotTime)));

        if isscalar(impulseValues)
            impulseValues = impulseValues * ones(size(impulsePlotTime));
        end

    end

    %% Построение графиков

    figure('Name', linkNames{linkIndex}, 'NumberTitle', 'off');

    subplot(1, 2, 1);

    plot(stepPlotTime, stepValues, 'LineWidth', 1.5);
    grid on;

    title('Переходная характеристика h(t)');
    xlabel('Время t, с');
    ylabel('h(t)');

    if hasStepImpulse
        subtitle('Аппроксимация импульсной составляющей');
    end

    subplot(1, 2, 2);

    plot(impulsePlotTime, impulseValues, 'LineWidth', 1.5);
    grid on;

    title('Импульсная характеристика w(t)');
    xlabel('Время t, с');
    ylabel('w(t)');

    if hasImpulseImpulse
        subtitle('Аппроксимация импульсной составляющей');
    end

    sgtitle(sprintf('%d. %s', linkIndex, linkNames{linkIndex}));

end

%% Проверка интегрирующего звена

integratingModel = tf(double(variantGain), [1 0]);

figure('Name', 'Проверка интегрирующего звена', 'NumberTitle', 'off');

subplot(1, 2, 1);

step(integratingModel, time);
grid on;

title('Проверка h(t) через step');
xlabel('Время t, с');
ylabel('h(t)');

subplot(1, 2, 2);

impulse(integratingModel, time);
grid on;

title('Проверка w(t) через impulse');
xlabel('Время t, с');
ylabel('w(t)');

%% Исследование влияния K

gainValues = double(variantGain) * [0.5, 1, 2];

figure('Name', 'Влияние K', 'NumberTitle', 'off');

hold on;
grid on;

for currentGain = gainValues

    response = currentGain * (1 - exp(-time / double(timeConstant)));

    plot(time, response, 'LineWidth', 1.5, 'DisplayName', ['K = ', num2str(currentGain)]);

end

title('Влияние K на апериодическое звено');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

%% Исследование влияния T

timeConstantValues = double(timeConstant) * [1/5, 1/2, 1, 2, 5];

experimentTime = 0:0.01:25;

figure('Name', 'Влияние T', 'NumberTitle', 'off');

hold on;
grid on;

for currentTimeConstant = timeConstantValues

    response = double(variantGain) * (1 - exp(-experimentTime / currentTimeConstant));

    plot(experimentTime, response, 'LineWidth', 1.5, 'DisplayName', ['T = ', num2str(currentTimeConstant)]);

end

title('Апериодическое звено 1-го порядка: влияние T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

%% Исследование влияния T1

firstCoefficientValues = double(firstCoefficient) * [1, 0.75, 0.5, 0.25];

oscillationTime = 0:0.02:80;

overshootValues = zeros(size(firstCoefficientValues));
settlingTimeValues = zeros(size(firstCoefficientValues));

figure('Name', 'Влияние T1', 'NumberTitle', 'off');

hold on;
grid on;

for coefficientIndex = 1:numel(firstCoefficientValues)

    currentFirstCoefficient = firstCoefficientValues(coefficientIndex);

    transferFunction = variantGain / (secondCoefficient^2*s^2 + sym(num2str(currentFirstCoefficient, 16))*s + 1);

    stepOutput = transferFunction * stepInput;
    stepResponse = simplify(ilaplace(stepOutput, s, t));

    responseFunction = matlabFunction(stepResponse, 'Vars', t);

    response = real(double(responseFunction(oscillationTime)));

    plot(oscillationTime, response, 'LineWidth', 1.5, 'DisplayName', ['T1 = ', num2str(currentFirstCoefficient)]);

    model = tf(double(variantGain), [double(secondCoefficient)^2, currentFirstCoefficient, 1]);

    characteristics = stepinfo(model, 'SettlingTimeThreshold', 0.02);

    overshootValues(coefficientIndex) = characteristics.Overshoot;

    settlingTimeValues(coefficientIndex) = characteristics.SettlingTime;

end

title('Колебательное звено: влияние T1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

%% Результаты исследования T1

fprintf('\nИсследование колебательного звена:\n');

fprintf('K = %.1f, T2 = %.1f\n', double(variantGain), double(secondCoefficient));

for coefficientIndex = 1:numel(firstCoefficientValues)

    fprintf('T1 = %.4f: перерегулирование = %.2f%%, время установления = %.2f с\n', firstCoefficientValues(coefficientIndex), overshootValues(coefficientIndex), settlingTimeValues(coefficientIndex));

end
