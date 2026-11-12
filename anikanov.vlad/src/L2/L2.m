syms s;

W = {
    1 / (s + 1)
    sym(2)
    5 / (s ^ 2 + 3 * s + 2)
};
W{end + 1} = simplify(W{2} * (W{1} + W{3}) / (1 + W{1} + W{2}));

[numerators, denominators] = cellfun(@numden, W, 'UniformOutput', false);
numerators = cellfun(@sym2poly, numerators, 'UniformOutput', false);
denominators = cellfun(@sym2poly, denominators, 'UniformOutput', false);

disp('Equivalent transfer function:');
disp(W{end});
disp('Numerator:');
disp(numerators{end});
disp('Denominator:');
disp(denominators{end});

model = 'scheme_L2';
open_system(fullfile(fileparts(mfilename('fullpath')), [model, '.slx']));
workspace = get_param(model, 'ModelWorkspace');
assignin(workspace, 'numerators', numerators);
assignin(workspace, 'denominators', denominators);

equivalent = tf(numerators{end}, denominators{end});
original = linearize(model, getlinio(model));
assert(norm(minreal(original - equivalent), inf) < 1e-9);
