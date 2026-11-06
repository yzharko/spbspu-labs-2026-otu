clc;

K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;

syms s t;
syms k tau tau1 tau2 positive;

time = (0:0.01:20)';
longTime = (0:0.005:70)';

names = {
    "Ideal integrator"
    "Ideal differentiator"
    "Proportional link"
    "First-order aperiodic link"
    "Second-order link: original parameters"
    "Oscillatory link"
    "Conservative link"
    "Inertial integrator"
    "Real differentiator"
    "First-order lead link"
    "Isodromic link"
};

W = {
    k / s
    k * s
    k
    k / (tau * s + 1)
    k / (tau2 ^ 2 * s ^ 2 + tau1 * s + 1)
    k / (tau2 ^ 2 * s ^ 2 + tau1 * s + 1)
    k / (tau * s ^ 2 + 1)
    k / (s * (tau * s + 1))
    k * s / (tau * s + 1)
    k * (tau * s + 1)
    k * (tau * s + 1) / s
};

inputs = {
    laplace(heaviside(t), t, s)
    laplace(dirac(t), t, s)
};

gains = K * [1 / 2, 1, 2];
timeConstants = T * [1 / 5, 1 / 2, 1, 2, 5];
T1Values = T1 * [1, 1 / 2, 1 / 5, 0];
aperiodicT1 = 3 * T2;

ordinates = {
    "h(t)"
    "w(t)"
};

tolerance = 1e-8;
maxStepError = 0;
maxImpulseError = 0;

for index = 1:numel(W)
    h = simplify(ilaplace(W{index} * inputs{1}, s, t));
    w = simplify(ilaplace(W{index} * inputs{2}, s, t));

    disp(names{index});
    disp(h);
    disp(w);

    transfer = subs(W{index}, {k, tau, tau1, tau2}, {K, T, T1, T2});
    [numerator, denominator] = numden(transfer);
    [polynomial, remainder] = quorem(numerator, denominator, s);
    model = tf(sym2poly(remainder), sym2poly(denominator));

    stepCheck = step(model, time) + double(subs(polynomial, s, 0));
    impulseCheck = impulse(model, time);

    stepCurve = regularResponse(transfer * inputs{1}, s, t, time);
    impulseCurve = regularResponse(transfer * inputs{2}, s, t, time);

    stepError = max(abs(stepCurve - stepCheck));
    impulseError = max(abs(impulseCurve - impulseCheck));
    assert(stepError < tolerance, 'Step response mismatch for %s.', names{index});
    assert(impulseError < tolerance, 'Impulse response mismatch for %s.', names{index});

    maxStepError = max(maxStepError, stepError);
    maxImpulseError = max(maxImpulseError, impulseError);

    figure('Name', names{index});
    tiledlayout(1, 2);

    for panel = 1:numel(inputs)
        nexttile;
        hold on;

        for gain = gains
            transfer = subs(W{index}, {k, tau, tau1, tau2}, {gain, T, T1, T2});
            curve = regularResponse(transfer * inputs{panel}, s, t, time);
            plot(time, curve, 'DisplayName', sprintf('K = %.2f', gain));
        end

        grid on;
        xlabel('t, s');
        ylabel(ordinates{panel});
        legend('Location', 'best');
    end
end

fprintf('Maximum step error: %.3e\n', maxStepError);
fprintf('Maximum impulse error: %.3e\n', maxImpulseError);
fprintf('Validation passed for %d links.\n', numel(W));

firstOrder = W{4};

figure('Name', 'First-order link: variation of T');
tiledlayout(1, 2);

for panel = 1:numel(inputs)
    nexttile;
    hold on;

    for timeConstant = timeConstants
        transfer = subs(firstOrder, {k, tau}, {K, timeConstant});
        curve = regularResponse(transfer * inputs{panel}, s, t, time);
        plot(time, curve, 'DisplayName', sprintf('T = %.2f', timeConstant));
    end

    grid on;
    xlabel('t, s');
    ylabel(ordinates{panel});
    legend('Location', 'best');
end

secondOrder = W{6};

figure('Name', 'Second-order link: variation of T1');
tiledlayout(1, 2);

for panel = 1:numel(inputs)
    nexttile;
    hold on;

    for T1Value = T1Values
        transfer = subs(secondOrder, {k, tau1, tau2}, {K, T1Value, T2});
        curve = regularResponse(transfer * inputs{panel}, s, t, longTime);
        plot(longTime, curve, 'DisplayName', sprintf('T1 = %.2f', T1Value));
    end

    grid on;
    xlabel('t, s');
    ylabel(ordinates{panel});
    legend('Location', 'best');
end

transfer = subs(secondOrder, {k, tau1, tau2}, {K, aperiodicT1, T2});

figure('Name', 'Second-order aperiodic example');
tiledlayout(1, 2);

for panel = 1:numel(inputs)
    nexttile;

    curve = regularResponse(transfer * inputs{panel}, s, t, time);
    plot(time, curve);

    grid on;
    xlabel('t, s');
    ylabel(ordinates{panel});
end

function curve = regularResponse(image, s, t, time)
    [numerator, denominator] = numden(image);
    [~, remainder] = quorem(numerator, denominator, s);

    regular = ilaplace(remainder / denominator, s, t);

    response = matlabFunction(regular, 'Vars', t);
    curve = response(time) + zeros(size(time));

    assert(all(isfinite(curve)), 'Response contains nonfinite values.');
    assert(max(abs(imag(curve))) < 1e-12 * max(1, max(abs(curve))), 'Response contains a nonzero imaginary part.');

    curve = real(curve);
end
