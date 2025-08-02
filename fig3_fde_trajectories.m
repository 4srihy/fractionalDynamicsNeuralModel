%% the code to generate alpha lambda matrix
% D^alpha (x) + lambda*x = B*randn(1)
figure
clear
%lambda = .9;
slopeFlag = 1;
numericalFlag = 0; % 0 : analytical method, 1: numerical method
nTrials = 200;
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

alpha = 1.1:0.04:1.9;
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

%         while max(abs(z1(iAlpha,iLambda,:)))>5 || abs(z1(iAlpha,iLambda,end))>0.01
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
             hfd0(iAlpha,iLambda,j) = HigFracDimV2(squeeze(y1(iAlpha,iLambda,j,tstart:end)),10,0.05,1,0,0);
          %   hurst(i) = genhurst(squeeze(y1(i,:)),1,30);
        end
         Power1(iAlpha,iLambda,:,:)= mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
       % Power1(iAlpha,iLambda,:) = squeeze(mean(abs(fft(squeeze(y1(iAlpha,iLambda,:,tstart:end))'))'.^2));
        
    end
end
hfd = mean(hfd0,3);
[~,f] = mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);

if tstart == 1
    t2 = t1.*tdur./tmax;
else
    t2 = (t1(tstart:end)-t1(tstart)).*tdur./(t1(end)-t1(tstart));
end
% df = 1/tdur;
% f = 0:df:(length(t2)-1)*df;
% 
% if length(Power1)~=length(f)
%     f = df:df:(length(t2))*df;
% end

%% computation of parameters
cf = zeros(length(alpha),length(lambda0));
ph = cf;
pw = cf;
cf00 = cf;
fooofRange = [64,140];
for iAlpha = 1:length(alpha)
    for iLambda = 1:length(lambda0)
        [ph0(iAlpha,iLambda),PosPh0{iAlpha,iLambda}] = max(log10(Power1(iAlpha,iLambda,1:end-5))); %removing the initial high power:adjust the no. 5 accordingly
        cf00(iAlpha,iLambda) = f(PosPh0{iAlpha,iLambda}(1)+0); %addition of 4 bcz max of power is taking starting from 5
        fd = fooof(f,squeeze(Power1(iAlpha,iLambda,:)),fooofRange,settings,true);
        [exponent(iAlpha,iLambda),~,~,cf0{iAlpha,iLambda},ph10{iAlpha,iLambda},pw0{iAlpha,iLambda}] = fooof_get_params(fd,squeeze(Power1(iAlpha,iLambda,:)),settings);
        if ~isempty(cf0{iAlpha,iLambda})
            cf(iAlpha,iLambda) = cf0{iAlpha,iLambda}{1};
            ph1(iAlpha,iLambda) = ph10{iAlpha,iLambda}{1};
            pw(iAlpha,iLambda) = pw0{iAlpha,iLambda}{1};
        end
    end
end

%% plots

hPlot = getPlotHandles(2,4,[0.06 0.1 0.85 0.83],0.07,0.16,0);
colorMarker = [0.72 0.27 1]; % for plots with markers and fitted line
colorLine = [0.44 0.15 0.60];
lineWidth = 2;
lineWidthDashed = 3;

%h1 = subplot('Position',[0.0800 0.8250 0.1625 0.1250]);
h1 = subplot('Position',[0.06500 0.5850 0.1625 0.3450]);
%h1 = hPlot(1);%subplot(4,4,1);
pcolor(h1,alpha,lambda0,hfd'); hold on;
shading('interp');
c = colorbar;
c.Label.String = 'HFD';  c.Label.FontWeight = 'bold'; c.FontSize = 12;
colormap(jet);
set(h1,'FontSize',14,'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off');
xlabel(h1,'\alpha','FontWeight','bold','FontSize',18);
ylabel(h1,'\lambda','FontWeight','bold','FontSize',18);
% t = title(h1,'FDE','Rotation',90,'FontSize',24);
% t.Position = [0.7720 -0.2707 0];


% 1. Baseline Aging
alphaPtsAge0 = [1.1 1.2 1.3 1.4];
alphaPtsAge = 1:2:11;%findClosestValue(alphaPtsAge0,alpha);
lambdaPtsAge = 2;%round(length(lambda0)/3);
dataPtsAge = [alphaPtsAge',repmat(lambdaPtsAge,[1,length(alphaPtsAge)])'];
colorNames2 = copper(2*length(alphaPtsAge));
labelData = 80:-15:0;

h3 = hPlot(3);%subplot(4,4,2);
makePredefinedTrajectories(h1, h3,f, Power1, alpha, lambda0, dataPtsAge,colorNames2,'^');
xlim([4 300]);
xline(h3,fooofRange(1),'--','LineWidth',2,'color',[0.3 0.3 0.3])
xline(h3,fooofRange(2),'--','LineWidth',2,'color',[0.3 0.3 0.3])
xlabel(h3,'Frequency (Hz)','FontWeight','bold');
ylabel(h3,'log_{10}Power','FontWeight','bold','FontSize',14);
set(h3,'XScale','log','XTick',[10 100],'FontSize',14,'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);
title(h3,'Baseline Activity','FontSize',16);

%slope vs hfd
for i = 1:length(dataPtsAge)
    exp1(i) = exponent(dataPtsAge(i,1),dataPtsAge(i,2));
    hfd1(i) = hfd(dataPtsAge(i,1),dataPtsAge(i,2));
end
mdl0 = fitlm(exp1,hfd1);
[c1,p1] = corrcoef(exp1,hfd1);
h31 = hPlot(2);%subplot(4,4,5);
axes(h31);
plot(exp1,hfd1,'x','MarkerSize',10,'LineStyle','none','LineWidth',2,'color',colorMarker); hold on;
plot(exp1,exp1.*mdl0.Coefficients.Estimate(2)+mdl0.Coefficients.Estimate(1),'LineWidth',lineWidthDashed,'LineStyle','--','color',colorLine);
xlabel('Slope','FontWeight','bold','FontSize',14);
ylabel('HFD','FontWeight','bold','FontSize',14);
set(h31,'FontSize',14,'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);
% h31.Position(1) = 0.08;
% h31.Position(3) = 0.18;
text(2.08,1.3,['\rho = ' num2str(c1(2),'%.2f')],'Fontsize',12,'FontWeight','bold');
text(2.08,1.2,['p = ' num2str(p1(2),'%.2e')],'Fontsize',12,'FontWeight','bold');

%age vs slope
mdl02 = fitlm(labelData,exp1);
h32 = hPlot(4);%subplot(4,4,6);
axes(h32);
plot(labelData,exp1,'x','MarkerSize',10,'LineStyle','none','LineWidth',2,'color',colorMarker); hold on;
plot(labelData,labelData.*mdl02.Coefficients.Estimate(2)+mdl02.Coefficients.Estimate(1),'LineWidth',lineWidthDashed,'LineStyle','--','color',colorLine);
xlabel('Age (years)','FontWeight','bold','FontSize',14);
ylabel('Slope','FontWeight','bold','FontSize',14);
set(h32,'FontSize',14,'XTick',[0 20 40 60 80],'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);
text(5,2.4,['\beta = ' num2str(mdl02.Coefficients.Estimate(2),3)],'Fontsize',12,'FontWeight','bold');
text(5,2.2,['p = ' num2str(mdl02.Coefficients.pValue(2),'%.2e')],'Fontsize',12,'FontWeight','bold');

%%%%%%%%%%%%%%% Stimulus Aging %%%%%%%%%%%%%%%%%%%%%%
alphaPtsAge0St = [1.1 1.2 1.3 1.4]+0.5;%+0.05*(1:4);
alphaPtsAgeSt = [11 15 20 21 21 21];
lambdaPtsAge0St = round(length(lambda0)/3);
lambdaPtsAgeSt = [3:6 8 10];
dataPtsAgeSt = [alphaPtsAgeSt', lambdaPtsAgeSt'];
colorNames2 = copper(2*length(alphaPtsAgeSt));

clear z1
for i = 1:length(dataPtsAgeSt)
    cf1(i) = cf00(dataPtsAgeSt(i,1),dataPtsAgeSt(i,2));
    ptsPower(i) = PosPh0{dataPtsAgeSt(i,1),dataPtsAgeSt(i,2)}(1);
    peakPowerChange(i) = sum(10*squeeze(log10(Power1(dataPtsAgeSt(i,1),dataPtsAgeSt(i,2),ptsPower(i)))-log10(Power1(alphaPtsAge(i),lambdaPtsAge,ptsPower(i)))));
    dhfdSt(i) = hfd(dataPtsAgeSt(i,1),dataPtsAgeSt(i,2))-hfd1(i);
    z1(i,:) = 10*squeeze(log10(Power1(dataPtsAgeSt(i,1),dataPtsAgeSt(i,2),:))-log10(Power1(alphaPtsAge(i),lambdaPtsAge,:)));
end

%purple [0.44 0.16 0.49]
h4 = hPlot(5);%subplot(4,4,3);

makePredefinedTrajectories(h1, h4,f, Power1, alpha, lambda0, dataPtsAgeSt, colorNames2,'d')
axes(h4);
for i=1:length(dataPtsAgeSt)
    plot(cf1(i),squeeze(log10(Power1(dataPtsAgeSt(i,1),dataPtsAgeSt(i,2),ptsPower(i)))),'Marker','.','MarkerSize',15,'Color','k'); hold on;
end
xlim([4 300]);
xlabel(h4,'Frequency (Hz)','FontWeight','bold','FontSize',14);
ylabel(h4,'log_{10}Power','FontWeight','bold','FontSize',14);
set(h4,'XScale','log','XTick',[10 100],'FontSize',14,'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);
title(h4,'Stimulus Activity','FontSize',16)

linkaxes([h3 h4]);

% change in Power
h5 = hPlot(7);%subplot(4,4,4);
axes(h5);
for i = 1:size(z1,1)
    plot(h5,f,z1(i,:),'lineWidth',lineWidth,'color',colorNames2(2*i,:,:)); hold on;
    plot(h5,cf1(i),peakPowerChange(i),'Marker','.','MarkerSize',15,'Color','k');
end
ylabel('\DeltaPower (dB)','FontSize',14,'FontWeight','bold');
xlabel('Frequency (Hz)','FontWeight','bold','FontSize',14);
title('Stimulus-Baseline','FontSize',16)
xlim([4 300]);
%legend(h4,num2str(labelData'),'Position',[0.9145 0.7757 0.0497 0.1490],'box','off');
legend(h4,num2str(labelData'),'Position',[0.9145 0.6757 0.0497 0.1490],'box','off');
set(h5,'XScale','log','XTick',[10 100],'FontSize',14,'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);


% Age vs Cf
mdl10 = fitlm(labelData,cf1);
[c2,p2] = corrcoef(labelData,cf1);
h41 = hPlot(6);%subplot(4,4,7);
axes(h41);
plot(labelData,cf1,'x','MarkerSize',10,'LineStyle','none','LineWidth',2,'color',colorMarker); hold on;
plot(labelData,labelData.*mdl10.Coefficients.Estimate(2)+mdl10.Coefficients.Estimate(1),'LineWidth',lineWidthDashed,'LineStyle','--','color',colorLine);
xlabel('Age (years)','FontWeight','bold','FontSize',14);
ylabel('CF (Hz)','FontWeight','bold','FontSize',14);
set(h41,'FontSize',14,'XTick',[0 20 40 60 80],'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);
text(5,25,['\rho = ' num2str(c2(2),2)],'Fontsize',12,'FontWeight','bold');
text(5,20,['p = ' num2str(p2(2),'%.2e')],'Fontsize',12,'FontWeight','bold');

%Age vs change in power and HFD
h42 = hPlot(8);%subplot(4,4,8);
axes(h42);
yyaxis left
plot(labelData,peakPowerChange,'x','MarkerSize',10,'LineStyle','none','LineWidth',2,'color',colorMarker); hold on;
mdl = fitlm(labelData,peakPowerChange,'poly2','Weights',[5 5 4 10 5 1]);
x = labelData;
plot(x,x.^2.*mdl.Coefficients.Estimate(3)+x.*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1),'lineWidth',lineWidthDashed,'LineStyle','--','color',colorLine);
xlabel('Age (years)','FontWeight','bold','FontSize',14);
ylabel('\Delta Power (dB)','FontWeight','bold','FontSize',14);
set(h42,'YColor',colorLine);
%h42.YAxis(1).Color = colorLine;
yyaxis right
plot(labelData,dhfdSt,'.','MarkerSize',20,'MarkerFaceColor',[1 0.5 0.5]); hold on;
mdl2 = fitlm(labelData,dhfdSt,'poly2');
plot(x,x.^2.*mdl2.Coefficients.Estimate(3)+x.*mdl2.Coefficients.Estimate(2)+mdl2.Coefficients.Estimate(1),'lineWidth',lineWidth,'LineStyle','-','color',[1 0.41 0.16]);
ylabel('\DeltaHFD');
%h42.YAxis(2).Color = [1 0.5 0.5];
set(h42,'FontSize',14,'XTick',[0 20 40 60 80],'FontWeight','bold','TickLength',[0.015 0.015],'TickDir','out','box','off','lineWidth',lineWidth);

a = annotation('textbox',[0.0276    0.9390    0.0228    0.0579],'String','A','FontSize',18,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[0.2564 0.9390 0.0228 0.0579],'String','B','FontSize',18,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.4732 0.9390 0.0228 0.0579],'String','C','FontSize',18,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[.6932 0.9390 0.0228 0.0579],'String','D','FontSize',18,'EdgeColor','none','FontWeight','bold');
e = annotation('textbox',[0.9238    0.9102    0.0228    0.0579],'String','Age (years)','FontSize',16,'EdgeColor','none','FontWeight','bold');

%%%%%%%%% Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makePredefinedTrajectories(h1, h2,f, Power1, alpha, lambda0, dataPts,colorNames,marker)
alphaVal = alpha(dataPts(:,1));
lambdaVal = lambda0(dataPts(:,2));

axes(h1); hold on;
for i = 1:size(dataPts,1)
    plot(alphaVal(i),lambdaVal(i),'Marker',marker,'MarkerEdgeColor','k','MarkerFaceColor',colorNames(2*i,:,:),'MarkerSize',6);
end

axes(h2); hold on;
for i = 1:size(dataPts,1)
    plot(f,log10(squeeze(Power1(dataPts(i,1),dataPts(i,2),:))),'LineWidth',2,'Color',colorNames(2*i,:,:)); hold on;
end

end

function index = findClosestValue(x,y)
%x: array you give
%y: array you want to find closest values in

for i = 1: length(x)
    [~,index(i)] = min(abs(y-x(i)));
end
end