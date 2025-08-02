% can be used to visualise fig 2 _ data and supp Fig _all electrode data
%%%%%% TLSAEEG data
elecGroup  = 'HP'; % can be 'HP' or 'all' for all electrodes

tData = load('TLSAEEGdata_eyesClosed_slope64to140.mat');
ageTLSA = tData.ageList;
if strcmp(elecGroup,'HP')
    slopeTLSA = tData.slope_HighPriority;
elseif strcmp(elecGroup,'All')
     slopeTLSA = tData.slope_allElec;
end
%%%%%%%% Voytek's Data
vData = load('VoytekData.mat');
ageVoytek = vData.age1;
slopeVoytek = vData.slope2';


hf = figure;
hf.Position = [360.3333 321.6667 818.6667 296.6667];
axTopo = [];    axParam=[]; axPSD = [];

age{1}  = ageTLSA;      slope{1} = slopeTLSA;   xPos(1) = 75;   yPos(1) = 1.9; dy(1) = 0.15;
age{2} = ageVoytek;     slope{2} = slopeVoytek; xPos(2) = 40;   yPos(2) = 2.4; dy(2) = 0.1;
labelType = [{'EEG'} {'ECoG'}];

for iType  = 1:2
        hPlot  = subplot(1,2,iType);
        scatterPlot(hPlot,age{iType},slope{iType},1,xPos(iType),yPos(iType),dy(iType)); hold on;
        xlabel(hPlot,'Age (years)','FontWeight','bold','FontSize',12);
        yLab  = ylabel(hPlot,'Slope','FontWeight','bold','FontSize',12);
        title(hPlot,labelType{iType})
        set(hPlot,'LineWidth',1,'FontSize',12,'FontWeight','bold','box','off');
end
a = annotation('textbox',[0.0264    0.9603    0.0228    0.0579],'String','A','FontSize',22,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[0.4771    0.9603    0.0228    0.0579],'String','B','FontSize',22,'EdgeColor','none','FontWeight','bold');


function scatterPlot(hPlot,x,y,fitFlag,xPos,yPos,dy)
if ~exist('xPos','var');    xPos = 75;  end
if ~exist('yPos','var');    yPos = 3;  end
if ~exist('dy','var');      dy = 0.2;  end
mdl = fitlm(x,y);
plot(hPlot,x,y,'.','Markersize',15,'LineStyle','none','color',[0.72 0.27 1]); hold on;

if fitFlag
    yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*x;
    plot(hPlot,x,yFit,'LineWidth',2,'color',[0.5 0.5 0.5]);
    text(xPos,yPos,['N = ' num2str(length(x))]);
    text(xPos,yPos-dy,['\beta = ' num2str(mdl.Coefficients.Estimate(2),3)]);
    if mdl.Coefficients.pValue(2)>0.001
        text(xPos,yPos-2*dy,['p-value: ' num2str(mdl.Coefficients.pValue(2)/2,'%.3f')]);
    else
        text(xPos,yPos-2*dy,['p-value: ' num2str(mdl.Coefficients.pValue(2)/2,'%.2e')]);
    end

     fit(x',y,'poly1');
end
end

