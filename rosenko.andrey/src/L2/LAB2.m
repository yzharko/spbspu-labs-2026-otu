
clc;
clear;
close all;

% W(s) = (s^3 + 4*s^2 + 50*s + 47) / (s^3 + 4*s^2 + 110*s + 147)
num = [1 4 50 47];
den = [1 4 110 147];

fprintf('\nLaboratory Work 2. Variant 10\n');
fprintf('Equivalent transfer function:\n');
fprintf('W(s) = (s^3 + 4*s^2 + 50*s + 47) / (s^3 + 4*s^2 + 110*s + 147)\n\n');
fprintf('Numerator:   ');
disp(num);
fprintf('Denominator: ');
disp(den);
fprintf('W(0) = %.9f\n', num(end) / den(end));
fprintf('System poles:\n');
disp(roots(den));

modelName = 'scheme_L2_variant10_R2018a';
modelFile = [modelName '.slx'];

if bdIsLoaded(modelName)
    close_system(modelName, 0);
end

if exist(modelFile, 'file') == 2
    delete(modelFile);
end

new_system(modelName);

set_param(modelName, ...
    'StopTime', '10', ...
    'SolverType', 'Variable-step', ...
    'Solver', 'ode45', ...
    'RelTol', '1e-7', ...
    'AbsTol', '1e-9', ...
    'MaxStep', '0.001');

addVariant10Lane(modelName, 'Step', ...
    'simulink/Sources/Step', 60, num, den);
addVariant10Lane(modelName, 'Pulse', ...
    'simulink/Sources/Pulse Generator', 280, num, den);
addVariant10Lane(modelName, 'Sine', ...
    'simulink/Sources/Sine Wave', 500, num, den);

set_param(modelName, 'ZoomFactor', 'FitSystem');
save_system(modelName, modelFile);

open_system(modelName);
open_system([modelName '/Step_Scope']);
open_system([modelName '/Pulse_Scope']);
open_system([modelName '/Sine_Scope']);

fprintf('\nModel created: %s\n', modelFile);
fprintf('Running simulation over the interval 0...10 s.\n');
sim(modelName);

fprintf('\nDone. Three scopes are open: Step_Scope, Pulse_Scope, Sine_Scope.\n');
fprintf('Each Scope should display two overlapping signals:\n');
fprintf('the original system output and the equivalent transfer function output.\n');
fprintf('Save screenshots of all three plots for the report.\n\n');

function addVariant10Lane(modelName, laneName, sourceLibraryPath, y, num, den)

    inputBlock = [laneName '_Input'];
    sum1 = [laneName '_Sum1'];
    w1 = [laneName '_W1'];
    sum2 = [laneName '_Sum2'];
    w2 = [laneName '_W2'];
    sum3 = [laneName '_Sum3'];
    w3 = [laneName '_W3'];
    w4 = [laneName '_W4'];
    equivalent = [laneName '_Equivalent'];
    mux = [laneName '_Mux'];
    scope = [laneName '_Scope'];

    add_block(sourceLibraryPath, [modelName '/' inputBlock], ...
        'Position', [30 y+45 65 y+75]);

    if strcmp(laneName, 'Step')
        set_param([modelName '/' inputBlock], ...
            'Time', '1', 'Before', '0', 'After', '1', 'SampleTime', '0');
    elseif strcmp(laneName, 'Pulse')
        set_param([modelName '/' inputBlock], ...
            'PulseType', 'Time based', ...
            'Amplitude', '1', ...
            'Period', '10', ...
            'PulseWidth', '5', ...
            'PhaseDelay', '0');
    elseif strcmp(laneName, 'Sine')
        set_param([modelName '/' inputBlock], ...
            'Amplitude', '1', ...
            'Bias', '0', ...
            'Frequency', '1', ...
            'Phase', '0', ...
            'SampleTime', '0');
    end

    add_block('simulink/Math Operations/Sum', [modelName '/' sum1], ...
        'Inputs', '+-', 'Position', [115 y+45 135 y+75]);
    add_block('simulink/Continuous/Transfer Fcn', [modelName '/' w1], ...
        'Numerator', '[1]', 'Denominator', '[1 1]', ...
        'Position', [180 y+45 235 y+75]);
    add_block('simulink/Math Operations/Sum', [modelName '/' sum2], ...
        'Inputs', '++-', 'Position', [270 y+45 290 y+75]);
    add_block('simulink/Math Operations/Gain', [modelName '/' w2], ...
        'Gain', '2', 'Position', [330 y+45 370 y+75]);
    add_block('simulink/Math Operations/Sum', [modelName '/' sum3], ...
        'Inputs', '++-', 'Position', [405 y+45 425 y+75]);
    add_block('simulink/Continuous/Transfer Fcn', [modelName '/' w3], ...
        'Numerator', '[5]', 'Denominator', '[1 3 2]', ...
        'Position', [465 y+45 525 y+75]);
    add_block('simulink/Math Operations/Gain', [modelName '/' w4], ...
        'Gain', '4', 'Position', [330 y+125 370 y+155]);

    add_block('simulink/Continuous/Transfer Fcn', [modelName '/' equivalent], ...
        'Numerator', mat2str(num), ...
        'Denominator', mat2str(den), ...
        'Position', [250 y-55 365 y-20]);

    add_block('simulink/Signal Routing/Mux', [modelName '/' mux], ...
        'Inputs', '2', 'Position', [570 y+48 575 y+92]);
    add_block('simulink/Sinks/Scope', [modelName '/' scope], ...
        'Position', [625 y+48 665 y+92]);

    add_line(modelName, [inputBlock '/1'], [sum1 '/1'], 'autorouting', 'on');
    add_line(modelName, [inputBlock '/1'], [equivalent '/1'], 'autorouting', 'on');

    add_line(modelName, [sum1 '/1'], [w1 '/1'], 'autorouting', 'on');
    add_line(modelName, [sum1 '/1'], [sum2 '/1'], 'autorouting', 'on');
    add_line(modelName, [sum1 '/1'], [sum3 '/1'], 'autorouting', 'on');
    hOriginal = add_line(modelName, [sum1 '/1'], [mux '/1'], 'autorouting', 'on');
    set_param(hOriginal, 'Name', 'Original');

    add_line(modelName, [w1 '/1'], [sum2 '/2'], 'autorouting', 'on');
    add_line(modelName, [sum2 '/1'], [w2 '/1'], 'autorouting', 'on');
    add_line(modelName, [w2 '/1'], [sum3 '/2'], 'autorouting', 'on');
    add_line(modelName, [sum3 '/1'], [w3 '/1'], 'autorouting', 'on');

    add_line(modelName, [w3 '/1'], [w4 '/1'], 'autorouting', 'on');
    add_line(modelName, [w3 '/1'], [sum3 '/3'], 'autorouting', 'on');
    add_line(modelName, [w4 '/1'], [sum1 '/2'], 'autorouting', 'on');
    add_line(modelName, [w4 '/1'], [sum2 '/3'], 'autorouting', 'on');

    hEquivalent = add_line(modelName, [equivalent '/1'], [mux '/2'], 'autorouting', 'on');
    set_param(hEquivalent, 'Name', 'Equivalent');
    add_line(modelName, [mux '/1'], [scope '/1'], 'autorouting', 'on');
end
