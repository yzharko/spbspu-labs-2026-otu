clc;
clear;
close all;
K = 1.6;
T = 0.8;
T1 = 1.0;
T2 = 1.1;
K_values = [0.5 0.9 K 3.2];
scales = [0.5 1 2];
varyA = {2*T2*[1 1.25 2], T1*[0.5 1 1.5]};
varyB = {T2*[0.5 0.75 1], T2*[0.75 1 1.5]};
omega = logspace(-2,2,1201);
dpi = 300;
scriptDir = fileparts(mfilename('fullpath'));
out = fullfile(scriptDir,'graphs');
if ~isfolder(out), mkdir(out); end
names = {'Интегрирующее','Дифференцирующее','Усилительное',...
    'Апериодическое I порядка','Апериодическое II порядка','Колебательное',...
    'Консервативное','Инерционное интегрирующее','Реальное дифференцирующее',...
    'Форсирующее','Изодромное'};
s = tf('s');
base = links(s,K,T,T1,T2);
colors = [0.58 0.33 0.16; 0.13 0.50 0.37; 0.06 0.29 0.51; 0.56 0.25 0.54];
for i = 1:11
    fig = figure('Color','w','Position',[80 80 1050 720]);
    lay = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    title(lay,sprintf('%d. %s',i,names{i}));
    axK = nexttile; hold(axK,'on');
    axT = nexttile; hold(axT,'on');
    axPhi = nexttile; hold(axPhi,'on');
    axN = nexttile; hold(axN,'on');
    for j = 1:numel(K_values)
        k = K_values(j);
        W = (k/K)*base{i};
        [w,L] = response(W,i,omega,T);
        semilogx(axK,w,L,'Color',colors(j,:),'LineWidth',1.5,...
            'DisplayName',sprintf('K = %g',k));
        cross = crossings(W);
        fprintf('Звено %d, K=%g: частоты 0 дБ, рад/с: %s\n',i,k,mat2str(cross,5));
    end
    La = asymptote(i,omega,K,T,T2);
    semilogx(axK,omega,La,'k:','LineWidth',1.2,'DisplayName','Асимптота, исходные');
    for j = 1:numel(scales)
        q = scales(j);
        if i<=3 && q~=1, continue; end
        group = links(s,K,q*T,q*T1,q*T2);
        [w,L,phase] = response(group{i},i,omega,q*T);
        name = sprintf('q = %g',q);
        semilogx(axT,w,L,'Color',colors(j+1,:),'LineWidth',1.5,'DisplayName',name);
        semilogx(axPhi,w,phase,'Color',colors(j+1,:),'LineWidth',1.5,'DisplayName',name);
    end
    [~,~,~,z] = response(base{i},i,omega,T);
    drawLocus(axN,z,i);
    marker = evalfr(base{i},1i);
    plot(axN,real(marker),imag(marker),'o','Color','k','MarkerFaceColor','w',...
        'DisplayName','\omega = 1');
    for ax = [axK axT axPhi]
        set(ax,'XScale','log'); grid(ax,'on'); xlim(ax,omega([1 end]));
        xlabel(ax,'\omega, рад/с');
    end
    ylabel(axK,'L, дБ'); ylabel(axT,'L, дБ'); ylabel(axPhi,'\phi, град');
    yline(axK,0,'Color',[.6 .6 .6],'HandleVisibility','off');
    yline(axT,0,'Color',[.6 .6 .6],'HandleVisibility','off');
    title(axK,'ЛАЧХ при изменении K'); title(axT,'ЛАЧХ при изменении времени');
    title(axPhi,'ЛФЧХ при изменении времени'); title(axN,'АФЧХ исходного варианта');
    for ax = [axK axT axN], legend(ax,'Location','best','FontSize',8); end
    print(fig,fullfile(out,sprintf('link_%02d.png',i)),'-dpng',sprintf('-r%d',dpi));
end

fig = figure('Color','w','Position',[80 80 1050 420]);
for j = 1:2
    ax = subplot(1,2,j); hold(ax,'on');
    if j==1, vals=K_values; else, vals=T*scales; end
    for n=1:numel(vals)
        if j==1
            W=vals(n)/(T*s+1); label=sprintf('K = %g',vals(n));
        else
            W=K/(vals(n)*s+1); label=sprintf('T = %g',vals(n));
        end
        [re,im]=nyquist(W,omega);
        plot(ax,re(:),im(:),'Color',colors(n,:),'LineWidth',1.5,'DisplayName',label);
        z=evalfr(W,1i);
        plot(ax,real(z),imag(z),'o','Color',colors(n,:),'HandleVisibility','off');
    end
    axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'Re W'); ylabel(ax,'Im W');
    title(ax,'АФЧХ при \omega > 0, метка при \omega = 1');
    legend(ax,'Location','best');
end
print(fig,fullfile(out,'first_order_locus.png'),'-dpng',sprintf('-r%d',dpi));

for n=1:2
    fig=figure('Color','w','Position',[80 80 1050 700]);
    for col=1:2
        axL=subplot(2,2,col); hold(axL,'on');
        axP=subplot(2,2,col+2); hold(axP,'on');
        if col==1, vals=varyA{n}; else, vals=varyB{n}; end
        for j=1:numel(vals)
            a=T1; b=T2;
            if n==1, a=2*T2; end
            if col==1, a=vals(j); else, b=vals(j); end
            W=K/(b^2*s^2+a*s+1);
            [w,L,ph]=response(W,n+4,omega,T);
            label=sprintf('T%d = %.3g с',col,vals(j));
            semilogx(axL,w,L,'Color',colors(j,:),'DisplayName',label);
            semilogx(axP,w,ph,'Color',colors(j,:));
        end
        for ax=[axL axP]
            set(ax,'XScale','log'); grid(ax,'on'); xlabel(ax,'\omega, рад/с');
        end
        ylabel(axL,'L, дБ'); ylabel(axP,'Фаза, град');
        title(axL,sprintf('Звено %d: изменение T%d',n+4,col));
        legend(axL,'Location','best');
    end
    print(fig,fullfile(out,sprintf('independent_%02d.png',n+4)),'-dpng',sprintf('-r%d',dpi));
end
fprintf('Графики сохранены в %s\n',out);

function W = links(s,k,t,a,b)

W = {k/s,k*s,tf(k),k/(t*s+1),k/(b*s+1)^2,k/(b^2*s^2+a*s+1),...
    k/(t*s^2+1),k/(s*(t*s+1)),k*s/(t*s+1),k*(t*s+1),k*(t*s+1)/s};
end

function [w,L,phase,z] = response(W,i,omega,t)

segments = {omega};
if i==7
    p=1/sqrt(t);
    segments={sort(unique([omega(omega<p*0.999),p*0.999])),...
        sort(unique([p*1.001,omega(omega>p*1.001)]))};
end
w=[]; L=[]; phase=[]; z=[];
for j=1:numel(segments)
    v=segments{j};
    [mag,phi]=bode(W,v);
    [re,im]=nyquist(W,v);
    h=re(:)+1i*im(:);
    assert(max(abs(mag(:)-abs(h))./max(1,abs(h)))<1e-8);
    phi=mod(phi(:)+180,360)-180;
    if j>1, w(end+1)=NaN; L(end+1)=NaN; phase(end+1)=NaN; z(end+1)=NaN; end
    w=[w(:); v(:)]; L=[L(:); 20*log10(mag(:))]; 
    phase=[phase(:); phi]; z=[z(:); h]; 
end
end

function L=asymptote(i,w,k,t,b)
x=20*log10(w); a=20*log10(max(1,t*w));
switch i
    case 1, shift=-x;
    case 2, shift=x;
    case 3, shift=zeros(size(w));
    case 4, shift=-a;
    case {5,6}, shift=-40*log10(max(1,b*w));
    case 7, shift=-40*log10(max(1,sqrt(t)*w));
    case 8, shift=-x-a;
    case 9, shift=x-a;
    case 10, shift=a;
    case 11, shift=a-x;
end
L=20*log10(k)+shift;
end

function c=crossings(W)

[n,d]=tfdata(W,'v');
n=n.*(1i).^(numel(n)-1:-1:0); d=d.*(1i).^(numel(d)-1:-1:0);
a=real(conv(n,conj(n))); b=real(conv(d,conj(d)));
m=max(numel(a),numel(b));
p=[zeros(1,m-numel(a)) a]-[zeros(1,m-numel(b)) b];
if all(p==0), c=Inf; return; end
r=roots(p); c=sort(real(r(abs(imag(r))<1e-7 & real(r)>1e-8))).';
end

function drawLocus(ax,z,i)
plot(ax,real(z),-imag(z),'--','Color',[.55 .55 .55],'DisplayName','\omega < 0');
plot(ax,real(z),imag(z),'Color',[.06 .29 .51],'LineWidth',1.7,'DisplayName','\omega > 0');
axis(ax,'equal'); grid(ax,'on'); xlabel(ax,'Re W'); ylabel(ax,'Im W');
switch i
    case {1,2}, xlim(ax,[-2 2]); ylim(ax,[-5 5]);
    case 7, xlim(ax,[-6 6]); ylim(ax,[-1.5 1.5]);
    case {8,10,11}, xlim(ax,[-2 3]); ylim(ax,[-5 5]);
    case 3, xlim(ax,[0 2]); ylim(ax,[-0.7 0.7]);
end
end
