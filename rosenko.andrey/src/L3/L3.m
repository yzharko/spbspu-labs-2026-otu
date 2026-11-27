clc;
clear;
close all;

K0 = 0.8;
T0 = 0.5;
T1 = 0.7;
T2 = 0.5;
T1a = 1.25;
K_values = [0.4 0.8 1.6];
T_values = [0.1 0.25 0.5 1.0 2.5];
w = logspace(-2,3,900);
s = tf('s');

names = {'Интегрирующее','Идеальное дифференцирующее','Усилительное', ...
    'Апериодическое 1-го порядка','Апериодическое 2-го порядка', ...
    'Колебательное','Консервативное','Инерционное интегрирующее', ...
    'Реальное дифференцирующее','Форсирующее 1-го порядка','Изодромное'};
outDir = 'L3_results';
if ~exist(outDir,'dir'), mkdir(outDir); end

baseTitles = {'ЛАЧХ типовых звеньев','ЛФЧХ типовых звеньев','АФЧХ типовых звеньев'};
gainTitles = {'Влияние усиления на ЛАЧХ','Влияние усиления на ЛФЧХ'};
timeTitles = {'Влияние постоянных времени на ЛАЧХ','Влияние постоянных времени на ЛФЧХ'};
files = {};

baseSheets = makeSheets(baseTitles);
for i = 1:11
    W = makeW(i,K0,T0,T1,T2,T1a,s);
    [ww,L,phase] = bodeData(W,w,i,T0);
    ax = subplot(4,3,i,'Parent',baseSheets(1));
    semilogx(ax,ww,L,'LineWidth',1.1); formatAxis(ax,w,'L(\omega), дБ',i,names);
    ax = subplot(4,3,i,'Parent',baseSheets(2));
    semilogx(ax,ww,phase,'LineWidth',1.1); formatAxis(ax,w,'\phi(\omega), °',i,names);
    if i == 7, markPhaseJump(ax,T0,[0 .32 .60]); end
    ax = subplot(4,3,i,'Parent',baseSheets(3));
    plotLocus(ax,W,i,K0,T0,w); title(ax,sprintf('%d. %s',i,names{i}),'FontSize',8);
end
for h = 1:3
    addPageLabels(baseSheets(h),h);
    files{end+1} = saveSheet(baseSheets(h),outDir,sprintf('base%d.pdf',h-1));
end

gainSheets = makeSheets(gainTitles);
crossings = cell(11,numel(K_values));
for i = 1:11
    for j = 1:numel(K_values)
        k = K_values(j);
        W = makeW(i,k,T0,T1,T2,T1a,s);
        [ww,L,phase] = bodeData(W,w,i,T0);
        ax = subplot(4,3,i,'Parent',gainSheets(1));
        if j == 1, hold(ax,'on'); end
        semilogx(ax,ww,L,'LineWidth',1.1,'DisplayName',sprintf('K = %g',k));
        formatAxis(ax,w,'L(\omega), дБ',i,names);
        ax = subplot(4,3,i,'Parent',gainSheets(2));
        if j == 1, hold(ax,'on'); end
        semilogx(ax,ww,phase,'LineWidth',1.1,'DisplayName',sprintf('K = %g',k));
        formatAxis(ax,w,'\phi(\omega), °',i,names);
        crossings{i,j} = zeroDbFrequencies(W,ww,L);
    end
end
for h = 1:2
    for i = 1:11
        ax = subplot(4,3,i,'Parent',gainSheets(h));
        legend(ax,'show','FontSize',5,'Location','best');
    end
    addPageLabels(gainSheets(h),h);
    files{end+1} = saveSheet(gainSheets(h),outDir,sprintf('gain%d.pdf',h-1));
end

fprintf('\nЧастоты пересечения уровня 0 дБ, рад/с\n');
fprintf('%-36s %14s %14s %14s\n','Звено','K = 0.4','K = 0.8','K = 1.6');
for i = 1:11
    fprintf('%-36s',names{i});
    for j = 1:3
        v = crossings{i,j};
        if isempty(v), value = '-'; else, value = sprintf('%.3f; ',v); value = value(1:end-2); end
        fprintf('%14s',value);
    end
    fprintf('\n');
end

timeSheets = makeSheets(timeTitles);
commonLinks = [4 7 8 9 10 11];
commonSlots = [1 6 7 8 9 10];
timeLabels = {'4. T','5. T_{1a}','5. T_2','6. T_1','6. T_2', ...
    '7. T','8. T','9. T','10. T','11. T'};

for j = 1:numel(commonLinks)
    i = commonLinks(j);
    slot = commonSlots(j);
    for T = T_values
        W = makeW(i,K0,T,T1,T2,T1a,s);
        [ww,L,phase] = bodeData(W,w,i,T);
        ax = subplot(4,3,slot,'Parent',timeSheets(1));
        if T == T_values(1), hold(ax,'on'); end
        semilogx(ax,ww,L,'LineWidth',1.1,'DisplayName',sprintf('T = %g',T));
        formatAxis(ax,w,'L(\omega), дБ',slot,timeLabels);
        ax = subplot(4,3,slot,'Parent',timeSheets(2));
        if T == T_values(1), hold(ax,'on'); end
        semilogx(ax,ww,phase,'LineWidth',1.1,'DisplayName',sprintf('T = %g',T));
        formatAxis(ax,w,'\phi(\omega), °',slot,timeLabels);
        if i == 7 && T == T_values(1), markPhaseJump(ax,T,[0 .32 .60]); end
    end
end

for i = [5 6]
    if i == 5
        firstValues = [1.0 1.25 1.5 2.0 2.5];
        secondValues = [.25 .4 .5 .6];
        firstBase = T1a;
    else
        firstValues = [.1 .3 .5 .7];
        secondValues = [.4 .5 .7 1.0];
        firstBase = T1;
    end
    slot = 2*(i-4);
    for v = firstValues
        W = makeW(i,K0,T0,v,T2,v,s);
        [ww,L,phase] = bodeData(W,w,i,T0);
        ax = subplot(4,3,slot,'Parent',timeSheets(1));
        if v == firstValues(1), hold(ax,'on'); end
        semilogx(ax,ww,L,'LineWidth',1.1,'DisplayName',sprintf('T_1 = %g',v));
        formatAxis(ax,w,'L(\omega), дБ',slot,timeLabels);
        ax = subplot(4,3,slot,'Parent',timeSheets(2));
        if v == firstValues(1), hold(ax,'on'); end
        semilogx(ax,ww,phase,'LineWidth',1.1,'DisplayName',sprintf('T_1 = %g',v));
        formatAxis(ax,w,'\phi(\omega), °',slot,timeLabels);
    end
    for v = secondValues
        W = makeW(i,K0,T0,firstBase,v,firstBase,s);
        [ww,L,phase] = bodeData(W,w,i,T0);
        ax = subplot(4,3,slot+1,'Parent',timeSheets(1));
        if v == secondValues(1), hold(ax,'on'); end
        semilogx(ax,ww,L,'LineWidth',1.1,'DisplayName',sprintf('T_2 = %g',v));
        formatAxis(ax,w,'L(\omega), дБ',slot+1,timeLabels);
        ax = subplot(4,3,slot+1,'Parent',timeSheets(2));
        if v == secondValues(1), hold(ax,'on'); end
        semilogx(ax,ww,phase,'LineWidth',1.1,'DisplayName',sprintf('T_2 = %g',v));
        formatAxis(ax,w,'\phi(\omega), °',slot+1,timeLabels);
    end
end

for h = 1:2
    for i = 1:11
        ax = subplot(4,3,i,'Parent',timeSheets(h));
        if ~isempty(findobj(ax,'Type','line')), legend(ax,'show','FontSize',5,'Location','best'); end
    end
    addPageLabels(timeSheets(h),h);
    files{end+1} = saveSheet(timeSheets(h),outDir,sprintf('time%d.pdf',h-1));
end
zip('L3_graphs.zip',files);

function f = makeSheets(titles)
f = gobjects(1,numel(titles));
for k = 1:numel(titles)
    f(k) = figure('Name',titles{k},'Visible','off','Color','w', ...
        'Units','inches','Position',[.3 .3 7.1 5.14], ...
        'PaperUnits','inches','PaperSize',[7.1 5.14], ...
        'PaperPosition',[0 0 7.1 5.14]);
end
end

function W = makeW(i,K,T,T1,T2,T1a,s)
switch i
    case 1, W = K/s;
    case 2, W = K*s;
    case 3, W = tf(K);
    case 4, W = K/(T*s+1);
    case 5, W = K/(T2^2*s^2+T1a*s+1);
    case 6, W = K/(T2^2*s^2+T1*s+1);
    case 7, W = K/(T*s^2+1);
    case 8, W = K/(s*(T*s+1));
    case 9, W = K*s/(T*s+1);
    case 10, W = K*(T*s+1);
    case 11, W = K*(T*s+1)/s;
    otherwise, error('Неизвестный номер звена: %d.',i);
end
end

function [w,L,phase] = bodeData(W,grid,i,T)
w = grid(:);
if i == 7
    w0 = 1/sqrt(T);
    d = logspace(-3,-1,16)';
    w = unique([w; w0*(1-d); w0; w0*(1+d)]);
    ok = w ~= w0;
    [mag,p] = bode(W,w(ok));
    L = nan(size(w)); phase = nan(size(w));
    L(ok) = 20*log10(squeeze(mag));
    phase(ok) = -abs(squeeze(p));
else
    [mag,p] = bode(W,w);
    L = 20*log10(squeeze(mag));
    phase = squeeze(p);
end
end

function formatAxis(ax,w,ylabelText,i,names)
grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log','XLim',w([1 end]),'FontSize',7);
xlabel(ax,'\omega, рад/с'); ylabel(ax,ylabelText);
if iscell(names) && numel(names) == 11
    title(ax,sprintf('%d. %s',i,names{i}),'FontSize',8);
else
    title(ax,names{i},'FontSize',8);
end
end

function plotLocus(ax,W,i,K,T,w)
[re,im] = nyquist(W,w);
re = re(:); im = im(:);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
switch i
    case {1,2}
        re(:) = 0; xlim(ax,[-1 1]);
    case 3
        plot(ax,K,0,'o','MarkerFaceColor',[0 .32 .60]);
        xlim(ax,[0 1.5*K]); ylim(ax,[-1 1]);
        xlabel(ax,'Re W'); ylabel(ax,'Im W'); return
    case 7
        lower = w(:) < 1/sqrt(T); upper = w(:) > 1/sqrt(T);
        plot(ax,re(lower),zeros(nnz(lower),1),'LineWidth',1.4);
        plot(ax,re(upper),zeros(nnz(upper),1),'LineWidth',1.4);
        xlim(ax,[-3*K 4*K]); ylim(ax,[-1 1]);
        xlabel(ax,'Re W'); ylabel(ax,'Im W'); return
    case 8
        extent = 1.06*max(abs(im));
        xlim(ax,[-1.16 .16]*K*T); ylim(ax,[-extent extent]);
        line(ax,[-K*T -K*T],[-extent extent],'Color',[.35 .35 .35],'LineStyle',':');
    case 10
        re(:) = K; xlim(ax,[0 1.6*K]);
    case 11
        re(:) = K*T; xlim(ax,[0 1.6*K*T]);
end
plot(ax,re,im,'LineWidth',1.2);
plot(ax,re,-im,'--','LineWidth',1.0);
if any(i == [4 5 9]), axis(ax,'equal'); end
xlabel(ax,'Re W'); ylabel(ax,'Im W');
end

function markPhaseJump(ax,T,color)
w0 = 1/sqrt(T);
line(ax,[w0 w0],[-180 0],'LineStyle',':','Color',color,'LineWidth',1.1, ...
    'HandleVisibility','off');
set(ax,'YLim',[-195 15],'YTick',[-180 -90 0]);
end

function addPageLabels(f,h)
if h == 1
    x = '\omega, рад/с'; y = 'L(\omega), дБ';
else
    x = '\omega, рад/с'; y = '\phi(\omega), °';
end
annotation(f,'textbox',[.4 .005 .3 .025],'String',x,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',8);
annotation(f,'textbox',[.005 .4 .025 .2],'String',y,'EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',8,'Rotation',90);
end

function c = zeroDbFrequencies(W,w,L)
idx = find(L(1:end-1).*L(2:end) < 0);
c = zeros(size(idx));
for j = 1:numel(idx)
    k = idx(j);
    c(j) = fzero(@(x) abs(evalfr(W,1i*x))-1,[w(k) w(k+1)]);
end
c = unique(c(:))';
end

function name = saveSheet(f,outDir,file)
name = fullfile(outDir,file);
print(f,name,'-dpdf','-painters');
close(f);
end
