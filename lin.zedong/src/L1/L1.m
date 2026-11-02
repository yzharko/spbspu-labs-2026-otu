%% Лабораторная работа №1. Временные характеристики динамических звеньев
% Запуск: кнопка «Выполнить» или команда L1.
% Результаты: results_L1/base и results_L1/extra рядом с этим файлом.

if isempty(ver('symbolic')) || exist('sym','file') ~= 2
    error('Требуется Symbolic Math Toolbox.');
end
if isempty(ver('control')) || exist('tf','file') ~= 2
    error('Требуется Control System Toolbox.');
end

p = struct('K',4,'T',1.5,'T1',1.6,'T2',2,'T1_link5',5);
assert(p.T1_link5 >= 2*p.T2,'Для звена 5 требуется T1 >= 2*T2.');
assert(p.T1 > 0 && p.T1 < 2*p.T2,'Для звена 6 требуется 0 < T1 < 2*T2.');
rootDir = fullfile(fileparts(mfilename('fullpath')),'results_L1');
baseDir = fullfile(rootDir,'base');
extraDir = fullfile(rootDir,'extra');
figDir = fullfile(baseDir,'figures');
smoothDir = fullfile(baseDir,'smooth_input_check');
extraFigDir = fullfile(extraDir,'figures');
folders = {figDir,smoothDir,extraFigDir};
for folderIndex=1:numel(folders)
    if ~exist(folders{folderIndex},'dir'), mkdir(folders{folderIndex}); end
end

[base,s,t] = l1_models(p,true);
comparisonLabels = {'Аналитическое вычисление','Численное моделирование'};
U_step = laplace(heaviside(t),t,s);
U_impulse = laplace(dirac(t),t,s);
fprintf('Лабораторная работа №1: K=%g; T=%g; T1=%g; T2=%g.\n',p.K,p.T,p.T1,p.T2);
fprintf('Звено 5: T1=%g >= 2*T2=%g; звено 6: 0 < T1=%g < 2*T2=%g.\n', ...
    p.T1_link5,2*p.T2,p.T1,2*p.T2);
fprintf('U_step(s)=%s; U_impulse(s)=%s.\n',char(U_step),char(U_impulse));
fprintf('Символьные вычисления могут занять некоторое время.\n');

%% 1. Аналитические отклики и независимая численная проверка
results = struct('parameters',p,'U_step',U_step,'U_impulse',U_impulse);
numbers = (1:11).';
names = cell(11,1);
stepError = zeros(11,1);
impulseError = zeros(11,1);
singularStep = zeros(11,1);
singularImpulseDelta = zeros(11,1);
singularImpulseDerivative = zeros(11,1);
textPath = fullfile(baseDir,'formulas.txt');
fid = fopen(textPath,'w','n','UTF-8');
if fid<0, error('Не удалось создать файл %s.',textPath); end
try
    fprintf(fid,'ЛАБОРАТОРНАЯ РАБОТА №1. Временные характеристики динамических звеньев\n');
    fprintf(fid,'Линь Цзэдун, 5130904/30106\n');
    fprintf(fid,'Личные параметры: K=%g; T=%g; T1=%g; T2=%g.\n',p.K,p.T,p.T1,p.T2);
    fprintf(fid,'Звено 5: T1=%g; T2=%g; T1 >= 2*T2.\n',p.T1_link5,p.T2);
    fprintf(fid,'Звено 6: T1=%g; T2=%g; 0 < T1 < 2*T2.\n',p.T1,p.T2);
    fprintf(fid,'U_step(s)=%s; U_impulse(s)=%s.\n\n',char(U_step),char(U_impulse));
    fprintf(fid,'Нулевые начальные условия; до t=0 все отклики равны нулю.\n');
    fprintf(fid,'Обычные части даны для t>0; на графиках при t=0 принят правый предел.\n');
    fprintf(fid,'delta имеет единичную площадь; delta'' является производной распределения.\n');
    fprintf(fid,'Сингулярные члены не являются конечными пиками. Ошибки в CSV относятся к обычным частям.\n\n');

    for j=1:11
        m = base(j);
        fprintf('[%02d/11] %s\n',j,m.name);
        Y_step = simplify(m.Wsym*U_step);
        h = simplify(ilaplace(Y_step,s,t));
        Y_impulse = simplify(m.Wsym*U_impulse);
        w = simplify(ilaplace(Y_impulse,s,t));

        % Полные h и w содержат распределения; для кривых выделяем обычные части.
        R_symbolic = simplify(m.Wsym-sym(m.a1,'r')*s-sym(m.a0,'r'));
        hRegular = simplify(sym(m.a0,'r')+ilaplace(R_symbolic/s,s,t));
        wRegular = simplify(ilaplace(R_symbolic,s,t));
        if ismember(j,[1 4 8 9 11])
            time = (0:0.01:10).';
        elseif j==7
            time = (0:0.01:25).';
        elseif ismember(j,[5 6])
            time = (0:0.01:40).';
        else
            time = (0:0.01:5).';
        end
        hAnalytic = l1_evaluate(hRegular,t,time);
        wAnalytic = l1_evaluate(wRegular,t,time);
        [hNumeric,wNumeric] = l1_response(m,time);
        stepError(j) = max(abs(hAnalytic-hNumeric));
        impulseError(j) = max(abs(wAnalytic-wNumeric));
        names{j} = m.name;
        singularStep(j) = m.a1;
        singularImpulseDelta(j) = m.a0;
        singularImpulseDerivative(j) = m.a1;

        fprintf(fid,'[%02d] %s\n',j,m.name);
        fprintf(fid,'Параметры: K=%g; T=%g; T1=%g; T2=%g.\n', ...
            m.parameters.K,m.parameters.T,m.parameters.T1,m.parameters.T2);
        if ~isempty(m.conditionNote), fprintf(fid,'%s\n',m.conditionNote); end
        fprintf(fid,'W(s)=%s\nY_step(s)=W(s)/s=%s\nh(t)=%s\n',char(m.Wsym),char(Y_step),char(h));
        fprintf(fid,'Y_impulse(s)=W(s)=%s\nw(t)=%s\n',char(Y_impulse),char(w));
        fprintf(fid,'h_regular(t)=%s; коэффициент при delta=%.12g\n',char(hRegular),m.a1);
        fprintf(fid,'w_regular(t)=%s; коэффициент при delta=%.12g; при delta''=%.12g\n', ...
            char(wRegular),m.a0,m.a1);
        fprintf(fid,'Максимальные абсолютные расхождения обычных частей: h=%.6g; w=%.6g\n\n', ...
            stepError(j),impulseError(j));

        fig = l1_figure(sprintf('%02d. %s',j,m.name),'base');
        l1_plot(subplot(1,2,1,'Parent',fig),time,[hAnalytic hNumeric],comparisonLabels, ...
            'h(t), отн. ед.','Переходная характеристика',11,m.a1,0,true);
        l1_plot(subplot(1,2,2,'Parent',fig),time,[wAnalytic wNumeric],comparisonLabels, ...
            'w(t), отн. ед./с','Импульсная характеристика',11,m.a0,m.a1,true);
        header = sprintf('%02d. %s | K=%g; T=%g; T1=%g; T2=%g', ...
            j,m.name,m.parameters.K,m.parameters.T,m.parameters.T1,m.parameters.T2);
        if ~isempty(m.conditionNote), header=sprintf('%s\n%s',header,m.conditionNote); end
        l1_export(fig,fullfile(figDir,m.fileName),header);
        results.links(j) = struct('number',j,'name',m.name,'parameters',m.parameters, ...
            'W',m.Wsym,'Y_step',Y_step,'h',h,'Y_impulse',Y_impulse,'w',w, ...
            'hRegular',hRegular,'wRegular',wRegular,'a1',m.a1,'a0',m.a0,'time',time, ...
            'hAnalytic',hAnalytic,'hNumeric',hNumeric,'wAnalytic',wAnalytic, ...
            'wNumeric',wNumeric,'conditionNote',m.conditionNote);
    end
    verification = table(numbers,names,stepError,impulseError,singularStep, ...
        singularImpulseDelta,singularImpulseDerivative,'VariableNames', ...
        {'Link','Name','MaxAbsError_h_regular','MaxAbsError_w_regular', ...
         'h_delta_coefficient','w_delta_coefficient','w_delta_derivative_coefficient'});
    writetable(verification,fullfile(baseDir,'verification.csv'));
    results.verification = verification;

    %% 2. Вспомогательная проверка звеньев с сингулярными откликами
    % Вход сглаживается фильтром Q; идеальные отклики остаются в results.links.
    epsilon = 0.05;
    epsilonSym = sym(epsilon,'r');
    Qsym = 1/(epsilonSym*s+1)^3;
    Qnum = tf(1,conv(conv([epsilon 1],[epsilon 1]),[epsilon 1]));
    smoothTime = (0:0.001:6).';
    special = [2 3 9 10 11];
    smoothStepError = zeros(numel(special),1);
    smoothImpulseError = zeros(numel(special),1);
    fprintf(fid,'ПРОВЕРКА СО СГЛАЖЕННЫМ ВХОДОМ\n');
    fprintf(fid,'Q(s)=1/(epsilon*s+1)^3; epsilon=%g с.\n',epsilon);
    fprintf(fid,'Аппроксимация импульса: q(t)=t^2*exp(-t/epsilon)/(2*epsilon^3), t>=0; интеграл=1.\n');
    fprintf(fid,'Аппроксимация ступени: 1-exp(-t/epsilon)*(1+t/epsilon+t^2/(2*epsilon^2)).\n\n');
    for r=1:numel(special)
        j = special(r);
        m = base(j);
        fprintf('[Сглаженный вход %d/%d] Звено %d, epsilon=%g с\n',r,numel(special),j,epsilon);
        WsmoothSymbolic = simplify(m.Wsym*Qsym);
        hSmooth = simplify(ilaplace(WsmoothSymbolic/s,s,t));
        wSmooth = simplify(ilaplace(WsmoothSymbolic,s,t));
        hA = l1_evaluate(hSmooth,t,smoothTime);
        wA = l1_evaluate(wSmooth,t,smoothTime);
        WsmoothNumeric = m.Wnum*Qnum;
        hN = reshape(step(WsmoothNumeric,smoothTime),[],1);
        wN = reshape(impulse(WsmoothNumeric,smoothTime),[],1);
        smoothStepError(r) = max(abs(hA-hN));
        smoothImpulseError(r) = max(abs(wA-wN));
        fig = l1_figure(sprintf('Сглаженный вход: звено %02d',j),'base');
        ah = subplot(1,2,1,'Parent',fig);
        l1_plot(ah,smoothTime,[hA hN],comparisonLabels, ...
            'y(t), отн. ед.','Отклик на сглаженную ступень',11,0,0,true);
        aw = subplot(1,2,2,'Parent',fig);
        l1_plot(aw,smoothTime,[wA wN],comparisonLabels, ...
            'y(t), отн. ед./с','Отклик на импульс единичной площади',11,0,0,true);
        if ismember(j,[2 3 10]), xlim(ah,[0 0.7]); xlim(aw,[0 0.7]); end
        header = sprintf(['Звено %02d: %s\n' ...
            'Вспомогательная проверка: Q(s)=1/(epsilon*s+1)^3, epsilon=%g с; сглаженный вход'],j,m.name,epsilon);
        l1_export(fig,fullfile(smoothDir,['smooth_' m.fileName]),header);
        fprintf(fid,'[%02d] h_smoothed(t)=%s\nw_smoothed(t)=%s\n',j,char(hSmooth),char(wSmooth));
        fprintf(fid,'Максимальные абсолютные расхождения: h=%.6g; w=%.6g\n\n', ...
            smoothStepError(r),smoothImpulseError(r));
        results.smooth(r) = struct('number',j,'epsilon',epsilon,'time',smoothTime, ...
            'hAnalytic',hA,'hNumeric',hN,'wAnalytic',wA,'wNumeric',wN);
    end
    smoothVerification = table(special(:),smoothStepError,smoothImpulseError, ...
        'VariableNames',{'Link','MaxAbsError_h_smoothed','MaxAbsError_w_smoothed'});
    writetable(smoothVerification,fullfile(baseDir,'smooth_verification.csv'));
    results.smoothVerification = smoothVerification;
    save(fullfile(baseDir,'results.mat'),'results');
catch caughtError
    fclose(fid);
    rethrow(caughtError);
end
fclose(fid);

fprintf('Проверка обычных частей откликов: максимальные абсолютные ошибки.\n');
for row=1:height(verification)
    fprintf('%02d. %s: h=%.6g; w=%.6g.\n',verification.Link(row), ...
        verification.Name{row},verification.MaxAbsError_h_regular(row),verification.MaxAbsError_w_regular(row));
end
fprintf('Проверка со сглаженным входом: максимальные абсолютные ошибки.\n');
for row=1:height(smoothVerification)
    fprintf('Звено %02d: h=%.6g; w=%.6g.\n',smoothVerification.Link(row), ...
        smoothVerification.MaxAbsError_h_smoothed(row),smoothVerification.MaxAbsError_w_smoothed(row));
end

%% 3. Задание со звёздочкой: влияние K на все 11 звеньев
fprintf('Исследование параметров: K, T, T1 и T2.\n');
extra = struct('baseline',p);
Kvalues = p.K*[0.5 1 2];
Kfig = [l1_figure('Влияние K: переходные характеристики','grid'), ...
        l1_figure('Влияние K: импульсные характеристики','grid')];
for j=1:11
    if ismember(j,[1 8 11]), time=(0:0.02:10).'; else, time=(0:0.02:40).'; end
    h=zeros(numel(time),3); w=h;
    stepD=zeros(1,3); impulseD=zeros(1,3); impulseDp=zeros(1,3);
    for r=1:numel(Kvalues)
        pr=p; pr.K=Kvalues(r); models=l1_models(pr,false);
        [h(:,r),w(:,r)]=l1_response(models(j),time);
        stepD(r)=models(j).a1; impulseD(r)=models(j).a0; impulseDp(r)=models(j).a1;
    end
    legends=arrayfun(@(v)sprintf('K=%g',v),Kvalues,'UniformOutput',false);
    l1_plot(subplot(4,3,j,'Parent',Kfig(1)),time,h,legends, ...
        'h(t), отн. ед.',sprintf('%02d. %s',j,base(j).name),9,stepD,zeros(1,3));
    l1_plot(subplot(4,3,j,'Parent',Kfig(2)),time,w,legends, ...
        'w(t), отн. ед./с',sprintf('%02d. %s',j,base(j).name),9,impulseD,impulseDp);
    extra.K(j)=struct('number',j,'baselineParameters',base(j).parameters, ...
        'values',Kvalues,'time',time,'hRegular',h,'wRegular',w, ...
        'hDelta',stepD,'wDelta',impulseD,'wDeltaDerivative',impulseDp);
end
Knotes={{'Влияние K: обычные части','переходных характеристик.', ...
    'Члены delta указаны в панелях.','Коэффициенты перечислены','в порядке легенды.'}, ...
    {'Влияние K: обычные части','импульсных характеристик.', ...
    'Члены delta и delta'' указаны','в панелях; их нельзя заменять','нулевой линией.'}};
Kstems={'01_K_step_all_links','02_K_impulse_all_links'};
for panel=1:2
    note=[Knotes{panel},{sprintf('Звено 5: T1=%g; T2=%g.',p.T1_link5,p.T2), ...
        sprintf('Звено 6: T1=%g; T2=%g.',p.T1,p.T2)}];
    ax=subplot(4,3,12,'Parent',Kfig(panel)); axis(ax,'off');
    text(ax,0,0.9,note,'Units','normalized','Interpreter','none', ...
        'FontSize',11,'VerticalAlignment','top','Color','k');
    l1_export(Kfig(panel),fullfile(extraFigDir,Kstems{panel}),'','grid');
end

%% 4. Постоянная времени звена первого порядка
Tvalues=p.T*[1/5 1/2 1 2 5];
time=(0:0.01:8*max(Tvalues)).';
h=zeros(numel(time),numel(Tvalues)); w=h;
for r=1:numel(Tvalues)
    pr=p; pr.T=Tvalues(r); models=l1_models(pr,false);
    [h(:,r),w(:,r)]=l1_response(models(4),time);
end
fig=l1_figure('Звено 4: влияние постоянной времени');
legends=arrayfun(@(v)sprintf('T=%g с',v),Tvalues,'UniformOutput',false);
l1_plot(subplot(1,2,1,'Parent',fig),time,h,legends,'h(t), отн. ед.','Ступенчатое воздействие',11);
l1_plot(subplot(1,2,2,'Parent',fig),time,w,legends,'w(t), отн. ед./с','Импульсное воздействие',11);
l1_export(fig,fullfile(extraFigDir,'03_first_order_T'), ...
    sprintf('Звено 4: K=%g; T/5, T/2, T, 2T, 5T; базовое T=%g с',p.K,p.T));
firstOrder=table(Tvalues(:),Tvalues(:),-log(0.05)*Tvalues(:),-log(0.02)*Tvalues(:), ...
    'VariableNames',{'T_seconds','Time_63p2_percent_seconds','Time_95_percent_seconds','Settling_2percent_seconds'});
writetable(firstOrder,fullfile(extraDir,'first_order_metrics.csv'));
extra.firstOrder=struct('Tvalues',Tvalues,'time',time,'h',h,'w',w,'metrics',firstOrder);

%% 5. Параметры T1 и T2 звеньев второго порядка
time=(0:0.02:100).';
parameterNames={'T1','T2'};
for j=[5 6]
    fig=l1_figure(sprintf('Звено %d: влияние T1 и T2',j),'four');
    if j==5
        groupName='aperiodic'; fields={'T1_link5','T2'};
        scanValues={[4 5 8],[1 2 2.5]};
        condition='T1 >= 2*T2 (включая критический случай)';
    else
        groupName='oscillatory'; fields={'T1','T2'};
        scanValues={p.T1*[0.5 1 2],p.T2*[0.5 1 2]};
        condition='0 < T1 < 2*T2';
    end
    extra.secondOrder.(groupName).number=j;
    extra.secondOrder.(groupName).baselineParameters=base(j).parameters;
    for col=1:2
        parameter=parameterNames{col}; values=scanValues{col};
        h=zeros(numel(time),numel(values)); w=h;
        curveParameters=repmat(base(j).parameters,1,numel(values));
        legends=cell(1,numel(values));
        for r=1:numel(values)
            pr=p; pr.(fields{col})=values(r); models=l1_models(pr,false);
            curveParameters(r)=models(j).parameters;
            actual=curveParameters(r);
            if j==5
                assert(actual.T1>=2*actual.T2,'Звено 5: требуется T1 >= 2*T2.');
            else
                assert(actual.T1>0 && actual.T1<2*actual.T2,'Звено 6: требуется 0 < T1 < 2*T2.');
            end
            % step/impulse корректно обрабатывают и повторные полюса при критическом демпфировании.
            [h(:,r),w(:,r)]=l1_response(models(j),time);
            legends{r}=sprintf('%s=%g',parameter,values(r));
            if j==5 && actual.T1==2*actual.T2, legends{r}=[legends{r} ' (крит.)']; end
        end
        fixedName=parameterNames{3-col}; fixedValue=base(j).parameters.(fixedName);
        l1_plot(subplot(2,2,col,'Parent',fig),time,h,legends,'h(t), отн. ед.', ...
            sprintf('Ступень: изменение %s (%s=%g)',parameter,fixedName,fixedValue),10);
        l1_plot(subplot(2,2,2+col,'Parent',fig),time,w,legends,'w(t), отн. ед./с', ...
            sprintf('Импульс: изменение %s (%s=%g)',parameter,fixedName,fixedValue),10);
        extra.secondOrder.(groupName).(parameter)=struct('values',values, ...
            'parameters',curveParameters,'time',time,'h',h,'w',w);
    end
    l1_export(fig,fullfile(extraFigDir,['04_second_order_' groupName '_T1_T2']), ...
        sprintf('Звено %d: %s | K=%g; %s',j,base(j).name,p.K,condition),'four');
end

%% 6. Влияние T на звенья 7–11
affected=[7 8 9 10 11];
Tvalues=p.T*[0.5 1 2];
for q=1:numel(affected)
    j=affected(q);
    if ismember(j,[8 11]), time=(0:0.01:20).'; else, time=(0:0.01:40).'; end
    h=zeros(numel(time),3); w=h;
    stepD=zeros(1,3); impulseD=zeros(1,3); impulseDp=zeros(1,3);
    for r=1:3
        pr=p; pr.T=Tvalues(r); models=l1_models(pr,false);
        [h(:,r),w(:,r)]=l1_response(models(j),time);
        stepD(r)=models(j).a1; impulseD(r)=models(j).a0; impulseDp(r)=models(j).a1;
    end
    fig=l1_figure(sprintf('Звено %d: влияние T',j));
    legends=arrayfun(@(v)sprintf('T=%g',v),Tvalues,'UniformOutput',false);
    l1_plot(subplot(1,2,1,'Parent',fig),time,h,legends, ...
        'h(t), отн. ед.','Ступенчатое воздействие',11,stepD,zeros(1,3));
    l1_plot(subplot(1,2,2,'Parent',fig),time,w,legends, ...
        'w(t), отн. ед./с','Импульсное воздействие',11,impulseD,impulseDp);
    l1_export(fig,fullfile(extraFigDir,sprintf('05_T_link_%02d',j)), ...
        sprintf('Звено %d: %s | K=%g; изменение T при остальных постоянных параметрах',j,base(j).name,p.K));
    extra.otherT(q)=struct('number',j,'values',Tvalues,'time',time,'hRegular',h, ...
        'wRegular',w,'hDelta',stepD,'wDelta',impulseD,'wDeltaDerivative',impulseDp);
end

%% 7. Уменьшение демпфирования и предел T1=0
T1values=[p.T1 p.T1/2 p.T1/5 0];
time=(0:0.01:150).';
h=zeros(numel(time),numel(T1values)); w=h;
zeta=T1values(:)/(2*p.T2);
overshoot=zeros(numel(T1values),1);
peakTime=zeros(numel(T1values),1);
settling2=nan(numel(T1values),1);
envelopeBound=nan(numel(T1values),1);
states=cell(numel(T1values),1);
for r=1:numel(T1values)
    pr=p; pr.T1=T1values(r); models=l1_models(pr,false);
    [h(:,r),w(:,r)]=l1_response(models(6),time);
    wd=sqrt(1-zeta(r)^2)/p.T2;
    peakTime(r)=pi/wd;
    overshoot(r)=100*exp(-pi*zeta(r)/sqrt(1-zeta(r)^2));
    if T1values(r)>0
        alpha=T1values(r)/(2*p.T2^2);
        envelopeBound(r)=-log(0.02*sqrt(1-zeta(r)^2))/alpha;
        settling2(r)=l1_settling(time,h(:,r),p.K,0.02);
        states{r}='Затухает';
    else
        states{r}='Незатухающие колебания; установления нет';
    end
end
fig=l1_figure('Звено 6: уменьшение T1');
legends=arrayfun(@(v)sprintf('T1=%g',v),T1values,'UniformOutput',false);
l1_plot(subplot(1,2,1,'Parent',fig),time,h,legends,'h(t), отн. ед.','Переходные характеристики',11);
l1_plot(subplot(1,2,2,'Parent',fig),time,w,legends,'w(t), отн. ед./с','Импульсные характеристики',11);
l1_export(fig,fullfile(extraFigDir,'06_oscillatory_T1_decrease'), ...
    sprintf('Звено 6: K=%g; T2=%g (фиксировано); T1=0 показан отдельно как предел без затухания',p.K,p.T2));
oscillatory=table(T1values(:),repmat(p.T2,numel(T1values),1),zeta,overshoot,peakTime, ...
    settling2,envelopeBound,states,'VariableNames', ...
    {'T1','T2','DampingRatio','PeakExcess_percent_of_K','FirstPeak_seconds', ...
     'Settling_2percent_seconds','EnvelopeBound_2percent_seconds','Behavior'});
writetable(oscillatory,fullfile(extraDir,'oscillatory_metrics.csv'));
extra.oscillatory=struct('T1values',T1values,'time',time,'h',h,'w',w,'metrics',oscillatory);
save(fullfile(extraDir,'extra_results.mat'),'extra');

%% 8. Параметры рисунков и итоговые показатели
textPath=fullfile(extraDir,'parameter_notes.txt');
fid=fopen(textPath,'w','n','UTF-8');
if fid<0, error('Не удалось создать файл %s.',textPath); end
try
    fprintf(fid,'ПАРАМЕТРЫ ИССЛЕДОВАНИЯ\nK=%g; T=%g; T1=%g; T2=%g.\n',p.K,p.T,p.T1,p.T2);
    fprintf(fid,'Звено 5: T1=%g; T2=%g; T1 >= 2*T2.\n',p.T1_link5,p.T2);
    fprintf(fid,'Звено 6: T1=%g; T2=%g; 0 < T1 < 2*T2.\n\n',p.T1,p.T2);
    fprintf(fid,'01_K_step_all_links и 02_K_impulse_all_links: K=[%s].\n',num2str(extra.K(1).values));
    for j=1:11
        fixed=extra.K(j).baselineParameters;
        fprintf(fid,'  Звено %d: T=%g; T1=%g; T2=%g; t=[0,%g] с.\n', ...
            j,fixed.T,fixed.T1,fixed.T2,extra.K(j).time(end));
    end
    fprintf(fid,'03_first_order_T: K=%g; T=[%s]; t=[0,%g] с.\n', ...
        p.K,num2str(extra.firstOrder.Tvalues),extra.firstOrder.time(end));
    groups={'aperiodic','oscillatory'};
    for q=1:numel(groups)
        group=extra.secondOrder.(groups{q});
        fprintf(fid,'04_second_order_%s_T1_T2: звено %d.\n',groups{q},group.number);
        for col=1:2
            study=group.(parameterNames{col});
            fprintf(fid,'  Изменение %s; t=[0,%g] с:\n',parameterNames{col},study.time(end));
            for r=1:numel(study.values)
                actual=study.parameters(r);
                fprintf(fid,'    K=%g; T1=%g; T2=%g; T1-2*T2=%g.\n', ...
                    actual.K,actual.T1,actual.T2,actual.T1-2*actual.T2);
            end
        end
    end
    for q=1:numel(extra.otherT)
        study=extra.otherT(q);
        fprintf(fid,'05_T_link_%02d: K=%g; T=[%s]; t=[0,%g] с.\n', ...
            study.number,p.K,num2str(study.values),study.time(end));
    end
    fprintf(fid,'06_oscillatory_T1_decrease: K=%g; T2=%g; T1=[%s]; t=[0,%g] с.\n\n', ...
        p.K,p.T2,num2str(extra.oscillatory.T1values),extra.oscillatory.time(end));
    notes={ ...
        'Звено 4: время до 63.2% равно T; до 95% равно -T*log(0.05); установление в полосе 2% равно -T*log(0.02).', ...
        'Звено 5: T1=2*T2 соответствует критическому демпфированию; эти случаи вычислены функциями step и impulse.', ...
        'Звено 6: zeta=T1/(2*T2); частота собственных колебаний 1/T2; показатель экспоненциальной огибающей T1/(2*T2^2).', ...
        'На рисунках с delta и delta'' показаны обычные части; коэффициенты распределений приведены отдельно.', ...
        'В oscillatory_metrics.csv: Settling_2percent_seconds — последнее пересечение границы +/-2% от K с линейной интерполяцией сетки 0.01 с.', ...
        'EnvelopeBound_2percent_seconds — верхняя оценка по огибающей; она отличается от последнего пересечения границы.', ...
        'При T1=0 установления нет, соответствующие времена равны NaN. Превышение 100% вычислено относительно K.'};
    for noteIndex=1:numel(notes), fprintf(fid,'%s\n\n',notes{noteIndex}); end
catch caughtError
    fclose(fid);
    rethrow(caughtError);
end
fclose(fid);

fprintf('Звено первого порядка: времена достижения уровня и установления.\n');
for row=1:height(firstOrder)
    fprintf('T=%.6g с; до 63,2%%: %.6g с; до 95%%: %.6g с; установление в полосе +/-2%%: %.6g с.\n', ...
        firstOrder.T_seconds(row),firstOrder.Time_63p2_percent_seconds(row), ...
        firstOrder.Time_95_percent_seconds(row),firstOrder.Settling_2percent_seconds(row));
end
fprintf('Колебательное звено: характеристики при уменьшении T1.\n');
for row=1:height(oscillatory)
    fprintf('T1=%.6g; T2=%.6g; относительное затухание=%.6g; превышение над K=%.6g%%; первый максимум: %.6g с.\n', ...
        oscillatory.T1(row),oscillatory.T2(row),oscillatory.DampingRatio(row), ...
        oscillatory.PeakExcess_percent_of_K(row),oscillatory.FirstPeak_seconds(row));
    fprintf('  Установление в полосе +/-2%%: %.6g с; оценка по огибающей: %.6g с; %s.\n', ...
        oscillatory.Settling_2percent_seconds(row),oscillatory.EnvelopeBound_2percent_seconds(row),oscillatory.Behavior{row});
end
fprintf('Выполнено. Результаты:\n%s\n',rootDir);
fprintf('Открыть пример рисунка:\nopenfig(fullfile(''%s'',''05_second_order_aperiodic.fig''),''visible'');\n',figDir);

%% Локальные функции
function [models,s,t]=l1_models(p,makeSymbolic)
% Независимые коэффициенты tf и символические выражения для всех 11 звеньев.
K=p.K; T=p.T; T1=p.T1; T2=p.T2; T1_link5=p.T1_link5;
names={'Интегрирующее (идеальное)','Дифференцирующее (идеальное)', ...
    'Усилительное (безынерционное)','Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка','Колебательное','Консервативное', ...
    'Инерционное интегрирующее','Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка','Изодромное'};
fileNames={'01_integrator','02_ideal_differentiator','03_gain', ...
    '04_first_order','05_second_order_aperiodic','06_oscillatory', ...
    '07_conservative','08_inertial_integrator','09_real_differentiator','10_forcing','11_isodromic'};
Wnum={tf(K,[1 0]),tf([K 0],1),tf(K,1),tf(K,[T 1]), ...
    tf(K,[T2^2 T1_link5 1]),tf(K,[T2^2 T1 1]),tf(K,[T 0 1]), ...
    tf(K,[T 1 0]),tf([K 0],[T 1]),tf([K*T K],1),tf([K*T K],[1 0])};
% W=a1*s+a0+R; R строго правильная или нулевая. h содержит a1*delta, w — a0*delta+a1*delta'.
a1=[0 K 0 0 0 0 0 0 0 K*T 0];
a0=[0 0 K 0 0 0 0 0 K/T K K*T];
Rnum={tf(K,[1 0]),tf(0),tf(0),tf(K,[T 1]), ...
    tf(K,[T2^2 T1_link5 1]),tf(K,[T2^2 T1 1]),tf(K,[T 0 1]), ...
    tf(K,[T 1 0]),tf(-K/T,[T 1]),tf(0),tf(K,[1 0])};
Wsym=cell(1,11); s=[]; t=[];
if makeSymbolic
    % Без предположения t>0: сингулярные члены при t=0 входят в полные формулы.
    syms s t
    k=sym(K,'r'); tau=sym(T,'r'); tau1=sym(T1,'r');
    tau2=sym(T2,'r'); tau1_link5=sym(T1_link5,'r');
    Wsym={k/s,k*s,k,k/(tau*s+1), ...
        k/(tau2^2*s^2+tau1_link5*s+1),k/(tau2^2*s^2+tau1*s+1), ...
        k/(tau*s^2+1),k/(s*(tau*s+1)),k*s/(tau*s+1),k*(tau*s+1),k*(tau*s+1)/s};
end
for j=1:11
    models(j).name=names{j}; models(j).fileName=fileNames{j};
    models(j).Wnum=Wnum{j}; models(j).Rnum=Rnum{j}; models(j).Wsym=Wsym{j};
    models(j).a1=a1(j); models(j).a0=a0(j);
    models(j).zeroRegular=ismember(j,[2 3 10]);
    models(j).parameters=struct('K',K,'T',T,'T1',T1,'T2',T2);
    models(j).conditionNote='';
    if j==5
        models(j).parameters.T1=T1_link5;
        models(j).conditionNote=sprintf('Звено 5: T1=%g; 2*T2=%g; условие T1 >= 2*T2.',T1_link5,2*T2);
    end
end
end

function [h,w]=l1_response(m,time)
% Обычные части откликов; коэффициенты распределений хранятся в a0 и a1.
if m.a1==0
    h=reshape(step(m.Wnum,time),[],1);
else
    h=m.a0+reshape(step(m.Rnum,time),[],1);
end
if m.zeroRegular, w=zeros(size(time)); else, w=reshape(impulse(m.Rnum,time),[],1); end
end

function y=l1_evaluate(expression,t,time)
% Правый предел при t=0 относится только к обычной части графика.
expression=subs(expression,heaviside(t),sym(1));
f=matlabFunction(expression,'Vars',t);
y=double(f(time));
if isscalar(y), y=repmat(y,size(time)); end
y=reshape(y,[],1);
if any(~isfinite(y)) || any(abs(imag(y))>1e-8*max(1,max(abs(y))))
    error('В обычной части отклика получено бесконечное, неопределённое или комплексное значение.');
end
y=real(y);
end

function fig=l1_figure(name,layout)
if nargin<2, layout='pair'; end
switch layout
    case 'grid', position=[50 50 1680 1400];
    case 'four', position=[60 60 1440 1000];
    case 'base', position=[100 100 1440 620];
    otherwise, position=[70 70 1440 620];
end
fig=figure('Name',name,'NumberTitle','off','Visible','off','Color','w','Position',position);
end

function l1_plot(ax,time,Y,legends,yLabel,caption,fontSize,delta,deltaPrime,isComparison)
% Сопоставление двух методов или семейство параметрических кривых.
if nargin<8, delta=0; end
if nargin<9, deltaPrime=0; end
if nargin<10, isComparison=false; end
if isComparison
    plot(ax,time,Y(:,1),'-','Color',[0 0.35 0.75],'LineWidth',1.9);
    hold(ax,'on');
    plot(ax,time,Y(:,2),'--','Color',[0.9 0.3 0.05],'LineWidth',1.2);
else
    set(ax,'ColorOrder',[0 0.35 0.75;0.9 0.3 0.05;0.25 0.6 0.15;0.6 0.2 0.65;0.2 0.65 0.7]);
    plot(ax,time,Y,'LineWidth',1.5);
end
grid(ax,'on');
set(ax,'Color','w','FontName','Arial','FontSize',fontSize,'XColor','k','YColor','k');
xlabel(ax,'Время t, с','Interpreter','none','Color','k');
ylabel(ax,yLabel,'Interpreter','none','Color','k');
if isComparison
    title(ax,caption,'Interpreter','none','FontWeight','normal','Color','k');
    legend(ax,legends,'Location','best','Interpreter','none','TextColor','k','Color','w');
    hold(ax,'off');
else
    title(ax,caption,'Interpreter','none','FontSize',fontSize,'FontWeight','normal','Color','k');
    legend(ax,legends,'Location','best','Interpreter','none','FontSize',fontSize,'TextColor','k','Color','w');
end
xlim(ax,[time(1) time(end)]);
if ~any(delta) && ~any(deltaPrime), return; end
limits=ylim(ax);
ylim(ax,[limits(1),limits(2)+0.45*(limits(2)-limits(1))]);
legend(ax,'Location','southeast');
parts={};
if isComparison
    if delta~=0, parts{end+1}=sprintf('%g*delta(t)',delta); end
    if deltaPrime~=0, parts{end+1}=sprintf('%g*delta''(t)',deltaPrime); end
    note=sprintf('%s(t) = обычная часть + %s\nПри t=0: распределение, не конечный пик.',yLabel(1),strjoin(parts,' + '));
    margin=5;
else
    if any(delta), parts{end+1}=sprintf('Коэф. delta: [%s]',strtrim(sprintf('%g ',delta))); end
    if any(deltaPrime), parts{end+1}=sprintf('Коэф. delta'': [%s]',strtrim(sprintf('%g ',deltaPrime))); end
    note=sprintf('%s: обычная часть + члены при t=0\n%s',yLabel(1),strjoin(parts,sprintf('\n')));
    margin=3;
end
text(ax,0.03,0.96,note,'Units','normalized','VerticalAlignment','top', ...
    'Interpreter','none','BackgroundColor',[1 0.96 0.85],'Color',[0.25 0.15 0], ...
    'FontSize',min(fontSize,10),'Margin',margin);
end

function l1_export(fig,stem,header,layout)
% Заголовок, PNG и редактируемый FIG для каждого рисунка.
if nargin<4, layout='pair'; end
if ~isempty(header)
    if strcmp(layout,'four')
        annotation(fig,'textbox',[0.05 0.94 0.9 0.05],'String',header, ...
            'Interpreter','none','HorizontalAlignment','center','EdgeColor','none','FontSize',12,'Color','k');
    else
        annotation(fig,'textbox',[0.04 0.89 0.92 0.095],'String',header, ...
            'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none', ...
            'EdgeColor','none','FontSize',12,'FontWeight','bold','Color','k');
        axs=findall(fig,'Type','axes');
        for k=1:numel(axs)
            pos=axs(k).Position; pos(2)=0.14; pos(4)=0.65; axs(k).Position=pos;
        end
    end
end
set(fig,'PaperPositionMode','auto','InvertHardcopy','off');
savefig(fig,[stem '.fig']);
print(fig,[stem '.png'],'-dpng','-r150');
close(fig);
end

function ts=l1_settling(time,y,target,relativeTolerance)
% Последнее пересечение границы полосы с линейной интерполяцией временной сетки.
tolerance=relativeTolerance*abs(target);
err=abs(y-target);
lastOutside=find(err>tolerance,1,'last');
if isempty(lastOutside)
    ts=time(1);
elseif lastOutside==numel(time)
    ts=NaN;
else
    i=lastOutside;
    fraction=(err(i)-tolerance)/(err(i)-err(i+1));
    ts=time(i)+fraction*(time(i+1)-time(i));
end
end
