%% generates figure for AutoRegressive model (Damped Harmonic Oscillator)
clear
slopeFlag = 1;

nTrials = 200;
dt = 0.1;
tmax  = 125; %keep it 125 only to get similar results as that of FDE
tstart = 1001;  %rejecting the y values before these time points
tdur = 0.5; % time duration for which we want to normalize t;
settings = struct();
settings.max_n_peaks = 1;
settings.peak_width_limits = [4,80];

Fs = round(((tmax/dt)-(tstart-1))/tdur);
params.tapers   = [1 1];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 200];
params.trialave = 1;
%% AR(2) Modifications
alpha = 1.8:0.01:1.99;%1.77:0.01:1.99; %==b1;
lambda0 = -(1:-0.001:0.9890);%==b2
N = tmax/dt;
x0 = [0 0];%[0.1 -0.5]; %initial conditions

b = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',1,'PassbandFrequency2',90,'PassbandRipple',0.2,'SampleRate',tmax/(dt*tdur));

for iLambda = 1:length(lambda0)
    for iAlpha = 1:length(alpha)
        for j = 1:nTrials
            y1(iAlpha,iLambda,j,:) = ar2([alpha(iAlpha),lambda0(iLambda)],x0,N);%ar2([alpha(iAlpha),lambda0(iLambda)],x0,N);
            hfd0(iAlpha,iLambda,j) = HigFracDimV2(squeeze(y1(iAlpha,iLambda,j,tstart:end)),10,0.05,1,0,0);
            %             y2 = filtfilt(b,squeeze(y1(iAlpha,iLambda,j,tstart:end)));
            %             hfd0(iAlpha,iLambda,j) = HigFracDimV2(y2,10,0.05,1,0,0);
            %   hurst(i) = genhurst(squeeze(y1(i,:)),1,30);
        end
        Power1(iAlpha,iLambda,:,:)= mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
        %Power1(iAlpha,iLambda,:) = squeeze(mean(abs(fft(squeeze(y1(iAlpha,iLambda,:,tstart:end))'))'.^2));
        hfd = mean(hfd0,3);
    end
end
[~,f] = mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
%Power1 = Power1-Power1(:,:,1); %offset correction

t1 = dt:dt:tmax;
t2 = (t1(tstart:end)-t1(tstart)).*tdur./(t1(end)-t1(tstart));%setdiff(t1,t1(1:tstart-1)).*tdur./t1(end);

% 
% df = 1/tdur;
% f = 0:df:(length(t2)-1)*df;
% 
% if length(Power1)~=length(f)
%     f = 0:df:(length(t2))*df;
% end


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
numRows = 5;
numCols = 4;
%hPlot = getPlotHandles(4,5,[0.1 0.07 0.86 0.87],0.06,0.06,0);
hPlotA{1} = getPlotHandles(2,2,[0.09 0.61 0.28 0.38], 0.06,0.02,0);
hPlotA{2} = getPlotHandles(2,2,[0.09 0.1 0.28 0.38], 0.06,0.02,0);
hPlot{1} = getPlotHandles(4,2,[0.51 0.1 0.24 0.88],0.02,0.03,0);
hPlot{2} = getPlotHandles(4,1,[0.83 0.1 0.15 0.88],0,0.03,0);
fontSize = 14;      fontWeight = 'bold';    lineWidth = 2;
lineWidthFooof = 2; %for showing fooof range
colorNames = turbo(10);


param{1} = cf00;        paramLabel{1} = 'CF';       clim{1} = [10 40];
param{2} = ph0;         paramLabel{2} = 'PH';       clim{2} = '';
param{3} = exponent;    paramLabel{3} = 'Slope';    clim{3} = [2  4];
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
        xlabel(h1,'\beta_1','FontWeight','bold','FontSize',fontSize);
    else
        set(h1,'XtickLabel',[]);
    end
    % if iParam == 1
    ylabel(h1,'\beta_2','FontWeight','bold','FontSize',fontSize);
    %end
    set(h1,'FontWeight','bold','FontSize',fontSize)
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
                plot( t2,squeeze(y1(xAlphaPos{i}(k),yAlphaPos{i}(iTer),1,tstart:end)),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
            end
            
            if iTer  == 1
                set(hTime(i,iTer),'XTickLabel',[]);
            else
                xlabel('Time (s)');
            end
            ylabel('Voltage');
            set(hTime(i,iTer),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth,'TickLength',[0.01 0.01],'TickDir','out');
            text(-0.3,-200,['\beta_2 = ' num2str(lambda0(yAlphaPos{i}(iTer)))],'Color',colorNames(5*iTer-1,:,:),'FontSize',fontSize+3,'FontWeight','bold','Rotation',90);

            axes(hPSD(i,iTer))
            for k = 1:4
                plot(f,log10(squeeze(Power1(xAlphaPos{i}(k),yAlphaPos{i}(iTer),:))),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
                label0{k} = [num2str(alpha(xAlphaPos{i}(k))) '(' num2str(hfd(xAlphaPos{i}(k),yAlphaPos{i}(iTer)),4) ')'];
                %title(['\beta_2 =' num2str(lambda0(2*j))]);
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
                set(hPSD(i,iTer),'XTick',[],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
            else
                xlabel('Frequency (Hz)');
                set(hPSD(i,iTer),'XTick',[10 100],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
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
               [l3,l4] = legend(label1,'Location','southwest','box','off','FontSize',fontSize-3); 
               for iTer = 1:2
                    l4(2*iTer+1).XData(1) = 0.2;
               end
           end
            if iParam == length(param)
                xlabel('\beta_1');
            else
                set(hPlotParam(i,iParam),'XTickLabel',[]);
            end
            ylabel(paramLabel{iParam});
            set(hPlotParam(i,iParam),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
        end
        a = annotation('textbox',[0.3813    0.9414    0.1    0.0579],'String','\beta_1 (HFD)','FontSize',fontSize-2,'EdgeColor','none','FontWeight','bold');

    else
        for iTer = 1:2
            hTime(i,iTer) = hPlotA{i}(iTer);%hPlot(2*i-1,2*iTer-1);
            hPSD(i,iTer) = hPlotA{i}(2+iTer);%hPlot(2*i-1,2*iTer);

            axes(hTime(i,iTer))
            for k = 1:4
                plot( t2, squeeze(y1(xAlphaPos{i}(iTer),yAlphaPos{i}(k),1,tstart:end)),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
            end
            if iTer  == 1
                set(hTime(i,iTer),'XTick',[]);
            else
                xlabel('Time (s)');
            end
            ylabel('Voltage');
            set(hTime(i,iTer),'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
            text(-0.3,-200,['\beta_1 = ' num2str(alpha(xAlphaPos{i}(iTer)))],'Color',colorNames(5*iTer-1,:,:),'FontSize',fontSize+3,'FontWeight','bold','Rotation',90);

            axes(hPSD(i,iTer))
            for k = 1:4
                plot(f,log10(squeeze(Power1(xAlphaPos{i}(iTer),yAlphaPos{i}(k),:))),'LineWidth',lineWidth,'color',colorNames(k+5*(iTer-1),:,:)); hold on;
                label0{k} = [num2str(lambda0(yAlphaPos{i}(k))) '(' num2str(hfd(xAlphaPos{i}(iTer),yAlphaPos{i}(k)),4) ')'];
                %title(['\beta_2 =' num2str(lambda0(2*j))]);
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
                set(hPSD(i,iTer),'XTick',[],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
            else
                xlabel('Frequency (Hz)');
                set(hPSD(i,iTer),'XTick',[10 100],'XScale','log','fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
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
                xlabel('\beta_2');
            else
                set(hPlotParam(i,iParam),'XTickLabel',[]);
            end
            %ylabel(paramLabel{iParam});
            set(hPlotParam(i,iParam),'YTick',[],'fontsize',fontSize,'FontWeight',fontWeight,'box','off','LineWidth',lineWidth);
        end
        annotation('textbox',[0.3813    0.4310    0.1    0.0579],'String','\beta_2 (HFD)','FontSize',fontSize-2,'EdgeColor','none','FontWeight','bold');
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

%% Ar(2) function
function x = ar2(b,x0,N)
x = x0;
for i=3:N
    x(i) = b(1)*x(i-1) + b(2)*x(i-2)+randn(1);
end
end

