syms s t K;

%% Параметры

T = sym(1);
T2 = sym(1);
aperiodicT1 = sym(3);
oscillatoryT1 = sym(7) / 10;

gainValues = 1:10;
time = 0:0.01:10;

%% Входные воздействия

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

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
    K / s;
    K * s;
    K;
    K / (T * s + 1);
    K / (T2^2 * s^2 + aperiodicT1 * s + 1);
    K / (T2^2 * s^2 + oscillatoryT1 * s + 1);
    K / (T * s^2 + 1);
    K / (s * (T * s + 1));
    K * s / (T * s + 1);
    K * (T * s + 1);
    K * (T * s + 1) / s
};

%% Аналитические вычисления

for linkIndex = 1:numel(transferFunctions)

    transferFunction = transferFunctions{linkIndex};

    stepOutput = simplify(transferFunction * U_step);
    stepResponse = simplify(ilaplace(stepOutput, s, t));

    impulseOutput = simplify(transferFunction * U_impulse);
    impulseResponse = simplify(ilaplace(impulseOutput, s, t));

    fprintf('\n%d. %s\n', linkIndex, linkNames{linkIndex});

    fprintf('W(s) = ');
    disp(transferFunction);

    fprintf('h(t) = ');
    disp(stepResponse);

    fprintf('w(t) = ');
    disp(impulseResponse);

    %% Отделение импульсных составляющих

    [stepNumerator, stepDenominator] = numden(stepOutput);

    [stepPolynomial, ~] = quorem(stepNumerator, stepDenominator, s);

    regularStep = simplify(ilaplace(stepOutput - stepPolynomial, s, t));

    [impulseNumerator, impulseDenominator] = numden(impulseOutput);

    [impulsePolynomial, ~] = quorem(impulseNumerator, impulseDenominator, s);

    regularImpulse = simplify(ilaplace(impulseOutput - impulsePolynomial, s, t));

    stepFunction = matlabFunction(regularStep, 'Vars', [t, K]);

    impulseFunction = matlabFunction(regularImpulse, 'Vars', [t, K]);

    %% Построение переходных характеристик

    figure('Name', linkNames{linkIndex}, 'NumberTitle', 'off');

    subplot(1, 2, 1);
    hold on;
    grid on;

    for gain = gainValues

        response = real(double(stepFunction(time, gain)));

        if isscalar(response)
            response = response * ones(size(time));
        end

        plot(time, response, 'DisplayName', ['K = ', num2str(gain)]);

    end

    title('Переходная характеристика h(t)');
    xlabel('Время t, с');
    ylabel('h(t)');

    if ~isequal(stepPolynomial, sym(0))
        text(0.03, 0.9, 'Импульсная составляющая при t = 0 не показана', 'Units', 'normalized', 'FontSize', 8);
    end

    legend('show', 'Location', 'best');
    hold off;

    %% Построение импульсных характеристик

    subplot(1, 2, 2);
    hold on;
    grid on;

    for gain = gainValues

        response = real(double(impulseFunction(time, gain)));

        if isscalar(response)
            response = response * ones(size(time));
        end

        plot(time, response, 'DisplayName', ['K = ', num2str(gain)]);

    end

    title('Импульсная характеристика w(t)');
    xlabel('Время t, с');
    ylabel('w(t)');

    if ~isequal(impulsePolynomial, sym(0))
        text(0.03, 0.9, 'Импульсная составляющая при t = 0 не показана', 'Units', 'normalized', 'FontSize', 8);
    end

    legend('show', 'Location', 'best');
    hold off;

    sgtitle(sprintf('%d. %s', linkIndex, linkNames{linkIndex}));

end

%% Проверка через Control System Toolbox

s_tf = tf('s');

figure('Name', 'Проверка интегрирующего звена', 'NumberTitle', 'off');

subplot(1, 2, 1);
hold on;
grid on;

for gain = gainValues

    transferFunctionCheck = gain / s_tf;
    response = squeeze(step(transferFunctionCheck, time));

    plot(time, response, 'DisplayName', ['K = ', num2str(gain)]);

end

title('Проверка h(t) через step');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest');
hold off;

subplot(1, 2, 2);
hold on;
grid on;

for gain = gainValues

    transferFunctionCheck = gain / s_tf;
    response = squeeze(impulse(transferFunctionCheck, time));

    plot(time, response, 'DisplayName', ['K = ', num2str(gain)]);

end

title('Проверка w(t) через impulse');
xlabel('Время t, с');
ylabel('w(t)');
legend('show', 'Location', 'best');
hold off;

%% Влияние постоянной времени T

fixedGain = 2;
experimentTime = 0:0.01:30;
timeConstantValues = [0.2, 0.5, 1, 2, 5];

figure('Name', 'Влияние T', 'NumberTitle', 'off');
hold on;
grid on;

for timeConstant = timeConstantValues

    transferFunction = K / (sym(timeConstant) * s + 1);

    stepOutput = transferFunction * U_step;
    stepResponse = ilaplace(stepOutput, s, t);

    responseFunction = matlabFunction(subs(stepResponse, K, fixedGain), 'Vars', t);

    response = real(double(responseFunction(experimentTime)));

    plot(experimentTime, response, 'DisplayName', ['T = ', num2str(timeConstant)]);

end

title('Апериодическое звено 1-го порядка: влияние T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

%% Влияние T1 на колебательное звено

dampingValues = [1.5, 1, 0.5, 0.1];

figure('Name', 'Влияние T1', 'NumberTitle', 'off');
hold on;
grid on;

for dampingCoefficient = dampingValues

    transferFunction = K / (T2^2 * s^2 + sym(dampingCoefficient) * s + 1);

    stepOutput = transferFunction * U_step;
    stepResponse = ilaplace(stepOutput, s, t);

    responseFunction = matlabFunction(subs(stepResponse, K, fixedGain), 'Vars', t);

    response = real(double(responseFunction(experimentTime)));

    plot(experimentTime, response,'DisplayName', ['T1 = ', num2str(dampingCoefficient)]);

end

title('Колебательное звено: влияние T1 при T2 = 1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;
