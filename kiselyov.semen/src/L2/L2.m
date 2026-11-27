function [t, curves, errors] = L2()
W1 = tf(1,[1 1]);
W2 = 2;
W3 = tf(5,[1 3 2]);
W4 = 4;
inner = feedback(W3,1+W2*W4);
forward = (1+W2*(1+W1))*inner;
equivalent = minreal(feedback(tf(1),W4*forward));
disp(equivalent);

original = ss([-1 -4 0; 0 0 1; 10 -107 -3],[1;0;15],[0 -4 0],1);
t = (0:0.001:12)';
names = {'Скачок','Прямоугольный импульс','Синусоида'};
curves = zeros(numel(t),2,3);
models = {original,equivalent};
for k = 1:2
    h = step(models{k},t);
    on = zeros(size(t));
    off = zeros(size(t));
    on(1001:end) = h(1:end-1000);
    off(1051:end) = h(1:end-1050);
    curves(:,k,1) = on;
    curves(:,k,2) = 20*(on-off);
    curves(:,k,3) = lsim(models{k},sin(t),t);
end

errors = zeros(3,1);
figure('Color','w');
for k = 1:3
    subplot(3,1,k);
    plot(t,curves(:,1,k),'-',t,curves(:,2,k),'--','LineWidth',1.3);
    grid on
    title(names{k});
    xlabel('Время, с');
    ylabel('Выход, отн. ед.');
    legend('Исходная схема','Эквивалентное звено');
    errors(k) = max(abs(curves(:,1,k)-curves(:,2,k)));
    fprintf('%s: %.3e\n',names{k},errors(k));
end
assert(all(errors < 1e-8),'Модели дают разные отклики.');

if exist('Simulink.SimulationInput','class')
    model = 'scheme_L2';
    load_system(fullfile(fileparts(mfilename('fullpath')),[model '.slx']));
    amplitudes = [1 0 0; 20 -20 0; 0 0 1];
    for k = 1:3
        cfg = Simulink.SimulationInput(model);
        cfg = cfg.setBlockParameter([model '/StepOn'],'After',num2str(amplitudes(k,1)));
        cfg = cfg.setBlockParameter([model '/StepOff'],'After',num2str(amplitudes(k,2)));
        cfg = cfg.setBlockParameter([model '/Sine'],'Amplitude',num2str(amplitudes(k,3)));
        cfg = cfg.setModelParameter('StopTime','12','Solver','ode4','FixedStep','0.001');
        out = sim(cfg);
        pair = out.get('outputs');
        difference = max(abs(pair.Data(:,1)-pair.Data(:,2)));
        fprintf('Simulink, %s: %.3e\n',names{k},difference);
        assert(difference < 1e-7,'Расхождение выходов в Simulink.');
    end
end
end
