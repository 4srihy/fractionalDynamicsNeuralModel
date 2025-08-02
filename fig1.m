%% the code to generate alpha lambda matrix and plot fig1
% D^alpha (x) + lambda*x = B*randn(1)
figure
clear

slopeFlag = 1;

nTrials = 200;
numericalFlag = 0; % 0 : analytical method, 1: numerical method

dt = 0.1;

if numericalFlag
    tmax  = 225;
    tstart = 1001;%rejecting the y values before these time points
else 
    precision0 = 1;
    tmax = 125;
    tstart  =1;
    t1 = dt:dt:tmax;
end
tdur = 0.5; % time duration for which we want to normalize t;
Fs = round(((tmax/dt)-(tstart-1))/tdur);
alpha =1.05:0.05:1.95;
lambda0 = 0.1:0.1:1;

settings = struct();
settings.max_n_peaks = 1;
settings.peak_width_limits = [4,80];

params.tapers   = [1 1];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 round(Fs/10)];
params.trialave = 1;

for iLambda = 1:length(lambda0)
    lambda = lambda0(iLambda);
    if slopeFlag
        f1= @(t,x)-lambda*x+randn(1);%x.*(4-y);
    else
        f1= @(t,x)-lambda*x+2*sin(0.8*t)+randn(1);%x.*(4-y);
    end
    for iAlpha = 1:length(alpha)
        if ~numericalFlag
        alpha1 = alpha(iAlpha);
        precision(iAlpha,iLambda) = precision0;
        z1(iAlpha,iLambda,:) = t1.^(alpha1-1).*mlf(alpha1,alpha1,-lambda.*t1.^alpha1,precision(iAlpha,iLambda));

% if z1 has diverging behaviour, set the precision using the following
%          while max(abs(z1(iAlpha,iLambda,:)))>5 || abs(z1(iAlpha,iLambda,end))>0.01
%             precision(iAlpha,iLambda) = precision(iAlpha,iLambda)-2;
%             z1(iAlpha,iLambda,:) = t1.^(alpha1-1).*mlf(alpha1,alpha1,-lambda.*t1.^alpha1,precision(iAlpha,iLambda));
%         end
        end
        for j = 1:nTrials
            if numericalFlag
                [t1,y1(iAlpha,iLambda,j,:)] = fde12(alpha(iAlpha),f1,dt,tmax,[0 0],dt);
            else
                 z2 = conv(randn(1,length(t1)),squeeze(z1(iAlpha,iLambda,:)));
                y1(iAlpha,iLambda,j,:) = z2(1:length(t1));            
            end
            hfd0(iAlpha,iLambda,j) = HigFracDimV2(squeeze(y1(iAlpha,iLambda,j,tstart:end)),5,0.05,1,0,0);
            %   hurst(i) = genhurst(squeeze(y1(i,:)),1,30);
        end
        Power1(iAlpha,iLambda,:,:)= mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
        %Power1(iAlpha,iLambda,:) = squeeze(mean(abs(fft(squeeze(y1(iAlpha,iLambda,:,tstart:end))'))'.^2));
       
    end
end
 hfd = mean(hfd0,3);
[~,f] = mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
%Power1 =10.^(log10(Power1)-log10(Power1(:,:,1))); %offset correction

if tstart == 1
    t2 = t1.*tdur./tmax;
else
    t2 = (t1(tstart:end)-t1(tstart)).*tdur./(t1(end)-t1(tstart));
end

% if Power is calculated as 'fft', then calculate frequency by this
% df = 1/tdur;
% f = df:df:length(t2)*df;
%
% if length(Power1)~=length(f)
%     f = df:df:(length(t2)+1)*df;
% end

%% computation of parameters
cf = zeros(length(alpha),length(lambda0));
ph = cf;
pw = cf;
cf00 = cf;
fooofRange = [64,140];
for iAlpha = 1:length(alpha)
    
    for iLambda = 1:length(lambda0)
        [ph0(iAlpha,iLambda),PosPh0] = max(log10(Power1(iAlpha,iLambda,1:end-5))); %removing the initial high power:adjust the no. 5 accordingly
        cf00(iAlpha,iLambda) = f(PosPh0(1));
        fd = fooof(f,squeeze(Power1(iAlpha,iLambda,:)),fooofRange,settings,true);
        errfit(iAlpha,iLambda) = fd.error;
        [exponent(iAlpha,iLambda),~,~,cf0{iAlpha,iLambda},ph10{iAlpha,iLambda},pw0{iAlpha,iLambda}] = fooof_get_params(fd,squeeze(Power1(iAlpha,iLambda,:)),settings);
        if ~isempty(cf0{iAlpha,iLambda})
            cf(iAlpha,iLambda) = cf0{iAlpha,iLambda}{1};
            ph1(iAlpha,iLambda) = ph10{iAlpha,iLambda}{1};
            pw(iAlpha,iLambda) = pw0{iAlpha,iLambda}{1};
        end
    end
end


%% plot figure
% hf = figure;
% hf.Position = [1.4877e+03 41.6667 560 599.3333];
numRows = 5;
numCols = 4;
hPlotA{1} = getPlotHandles(2,2,[0.09 0.61 0.28 0.38], 0.06,0.02,0);
hPlotA{2} = getPlotHandles(2,2,[0.09 0.1 0.28 0.38], 0.06,0.02,0);
hPlot{1} = getPlotHandles(4,2,[0.52 0.1 0.24 0.88],0.01,0.03,0);
hPlot{2} = getPlotHandles(4,1,[0.83 0.1 0.15 0.88],0,0.03,0);
fontSize = 14;      fontWeight = 'bold';    lineWidth = 2;
lineWidthFooof = 2; %for showing fooof range
colorNames = turbo(10);

param{1} = exponent;    paramLabel{1} = 'Slope';    clim{1} = [2  4];
param{2} = cf00;        paramLabel{2} = 'CF';       clim{2} = [10 40];
param{3} = ph0;         paramLabel{3} = 'PH';       clim{3} = '';
param{4} = hfd;         paramLabel{4} = 'HFD';      clim{4} = [1 1.5];

for iParam = 1:length(param)
    h1 = hPlot{2}(iParam);
    axes(h1);
    pcolor(alpha,lambda0,param{iParam}');
    colormap(jet); shading('flat');
    c = colorbar;
    c.Label.String = paramLabel{iParam};
    if ~isempty(clim{iParam})
        caxis(clim{iParam});
    end
    if iParam==length(param)
        xlabel(h1,'\alpha','FontWeight','bold','FontSize',fontSize);
    else
        set(h1,'Xtick',[]);
    end
    % if iParam == 1
    ylabel(h1,'\lambda','FontWeight','bold','FontSize',fontSize);
    %end
    set(h1,'FontWeight','bold','FontSize',fontSize,'lineWidth',lineWidth,'TickLength',[0.02 0.02],'TickDir','out','box','off');
end

k = 1:4; iTer = 1:2;
xAlphaPos{1} = 5.*k-1;  yAlphaPos{1} = 5.*iTer-2;
xAlphaPos{2} = 15.*iTer-12;  yAlphaPos{2} = 2.*k;

for i = 1:2 % used for alpha and lambda
    if i==1
        for iTer = 1:2
            hTime(i,iTer) = hPlotA{i}(iTer);%hPlot(2*i-1,2*iTer-1);
            hPSD(i,iTer) = hPlotA{i}(2+iTer);%hPlot(2*i-1,2*iTer);

            axes(hTime(i,iTer))
            for k = 1:4
                plot( t2,squeeze(y1(xAlphaPos{i}(k),yAlphaPos{i}(iTer),5,tstart:end)),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
            end

            if iTer  == 1
                set(hTime(i,iTer),'XTickLabel',[]);
            else
                xlabel('Time (s)');
            end
            ylabel('Voltage');
            set(hTime(i,iTer),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            text(-0.25,-4,['\lambda = ' num2str(lambda0(yAlphaPos{i}(iTer)))],'Color',colorNames(5*iTer-1,:,:),'FontSize',fontSize+3,'FontWeight','bold','Rotation',90);

            axes(hPSD(i,iTer))
            for k = 1:4
                plot(f,log10(squeeze(Power1(xAlphaPos{i}(k),yAlphaPos{i}(iTer),:))),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
                label0{k} = [num2str(alpha(xAlphaPos{i}(k))) '(' num2str(hfd(xAlphaPos{i}(k),yAlphaPos{i}(iTer)),4) ')'];
                %title(['\lambda =' num2str(lambda0(2*j))]);
            end
            for k=1:4
                cfPos = f==cf00(xAlphaPos{i}(k),yAlphaPos{i}(iTer));
                plot(f(cfPos),squeeze(log10(Power1(xAlphaPos{i}(k),yAlphaPos{i}(iTer),cfPos))),'.k','MarkerSize',10);
            end
            xline(fooofRange(1),'--','lineWidth',lineWidthFooof);
            xline(fooofRange(2),'--','lineWidth',lineWidthFooof);
            xlim([4 200]);
            [l1,l2] =  legend(label0,'Location','southeast','box','off','FontSize',fontSize-3);
            for k = 1:4
                l2(2*k+3).XData(1) = 0.2;
            end
            l1.Position(1) = l1.Position(1)+0.08;
            ylabel('log_{10}Power')
            if iTer  == 1
                set(hPSD(i,iTer),'XTick',[],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            else
                xlabel('Frequency (Hz)');
                set(hPSD(i,iTer),'XTick',[10 100],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            end

        end
        for iParam  = 1:length(param)
            hPlotParam(i,iParam) = hPlot{1}(iParam,i);
            axes(hPlotParam(i,iParam))
            for iTer = 1:2
                yy1  = smooth(param{iParam}(1:end,yAlphaPos{i}(iTer)));
                plot(alpha(1:end),yy1,'LineWidth',2.5,'color',colorNames(5*iTer-1,:,:)); hold on; %need to change color
                label1{iTer} = num2str(lambda0(yAlphaPos{i}(iTer)));
            end
            if iParam==4
                hfd2 = 2-alpha./2;
                plot(alpha,hfd2,'--','lineWidth',lineWidth,'Color','k');
                
                [l3,l4] = legend(label1,'Location','southwest','box','off','FontSize',fontSize-3);
                for iTer = 1:2
                    l4(2*iTer+1).XData(1) = 0.2;
                end
            end
            if iParam == length(param)
                xlabel('\alpha');
            else
                set(hPlotParam(i,iParam),'XTickLabel',[]);
            end
            ylabel(paramLabel{iParam});
            set(hPlotParam(i,iParam),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
        end
        a = annotation('textbox',[0.3813    0.9414    0.1    0.0579],'String','\alpha (HFD)','FontSize',fontSize-2,'EdgeColor','none','FontWeight','bold');

    else
        for iTer = 1:2
            hTime(i,iTer) = hPlotA{i}(iTer);%hPlot(2*i-1,2*iTer-1);
            hPSD(i,iTer) = hPlotA{i}(2+iTer);%hPlot(2*i-1,2*iTer);

            axes(hTime(i,iTer))
            for k = 1:4
                plot( t2, squeeze(y1(xAlphaPos{i}(iTer),yAlphaPos{i}(k),5,tstart:end)),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
            end
            if iTer  == 1
                set(hTime(i,iTer),'XTick',[]);
            else
                xlabel('Time (s)');
            end
            ylabel('Voltage');
            set(hTime(i,iTer),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            text(-0.25,-4,['\alpha = ' num2str(alpha(xAlphaPos{i}(iTer)))],'Color',colorNames(5*iTer-1,:,:),'FontSize',fontSize+3,'FontWeight','bold','Rotation',90);

            axes(hPSD(i,iTer))
            for k = 1:4
                plot(f,log10(squeeze(Power1(xAlphaPos{i}(iTer),yAlphaPos{i}(k),:))),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
                label0{k} = [num2str(lambda0(yAlphaPos{i}(k))) '(' num2str(hfd(xAlphaPos{i}(iTer),yAlphaPos{i}(k)),4) ')'];
                %title(['\lambda =' num2str(lambda0(2*j))]);
            end
            for k=1:4
                cfPos = f==cf00(xAlphaPos{i}(iTer),yAlphaPos{i}(k));
                plot(f(cfPos),squeeze(log10(Power1(xAlphaPos{i}(iTer),yAlphaPos{i}(k),cfPos))),'.k','MarkerSize',10);
            end
            xline(fooofRange(1),'--','lineWidth',lineWidthFooof);
            xline(fooofRange(2),'--','lineWidth',lineWidthFooof);
            xlim([4 200]);
            [l1,l2] =  legend(label0,'Location','southeast','box','off','FontSize',fontSize-3);

            l1.Position(1) = l1.Position(1)+0.08;
            for k = 1:4
                l2(2*k+3).XData(1) = 0.2;
            end
            ylabel('log_{10}Power');
            if iTer  == 1
                set(hPSD(i,iTer),'XTick',[],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            else
                xlabel('Frequency (Hz)');
                set(hPSD(i,iTer),'XTick',[10 100],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            end

        end
        for iParam  = 1:length(param)
            hPlotParam(i,iParam) = hPlot{1}(iParam,i);
            axes(hPlotParam(i,iParam))
            for iTer = 1:2
                yy1  = smooth(param{iParam}(xAlphaPos{i}(iTer),:));
                plot(lambda0,yy1,'LineWidth',2.5,'color',colorNames(5*iTer-1,:,:)); hold on; %need to change color
                label1{iTer} = num2str(alpha(xAlphaPos{i}(iTer)));
            end
            if iParam==4
                [l3,l4] = legend(label1,'Location','southwest','box','off','FontSize',fontSize-3);
                for iTer = 1:2
                    l4(2*iTer+1).XData(1) = 0.2;
                end
            end
            if iParam==length(param)
                xlabel('\lambda');
            else
                set(hPlotParam(i,iParam),'XTickLabel',[]);
            end
            %ylabel(paramLabel{iParam});
            set(hPlotParam(i,iParam),'YTick',[],'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
        end
        annotation('textbox',[0.3813    0.4310    0.1    0.0579],'String','\lambda (HFD)','FontSize',fontSize-2,'EdgeColor','none','FontWeight','bold');
    end
end

linkaxes(hTime);
linkaxes(hPSD);
for iParam = 1:length(param)
    linkaxes(hPlotParam(:,iParam),'y');
end
a = annotation('textbox',[0.008    0.9524    0.0228    0.0579],'String','A','FontSize',18,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[0.008    0.4793    0.0228    0.0579],'String','B','FontSize',18,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[0.4470    0.9524    0.0228    0.0579],'String','C','FontSize',18,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[0.7693   0.9524    0.0228    0.0579],'String','D','FontSize',18,'EdgeColor','none','FontWeight','bold');

