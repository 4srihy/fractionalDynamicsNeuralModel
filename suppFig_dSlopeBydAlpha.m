% Sup fig : showing the need of optimum range that satisfies theory

%% the code to generate alpha lambda matrix
% D^alpha (x) + lambda*x = B*randn(1)
%figure
clear
%clf
%lambda = .9;
alpha = 1.1:0.018:1.5;
lambda0 = 0.1:0.1:1;
slopeFlag = 1;
TLSAFlag = 1; %if ~TLSAFlag, then Voytek
numericalFlag = 0; % 0 : analytical method, 1: numerical method

nTrials = 200;
%Fs = 250; Fs = tmax/(dt*tdur);
dt = 0.1;

if numericalFlag
    tmax  = 225;
    tstart = 1001;%round(0.2*(dt*tmax))+900;  l%rejecting the y values before these time points
else
    precision0 = 1;
    tmax = 125;
    tstart  =1;
    t1 = dt:dt:tmax;
end

if TLSAFlag
    %lambda0 = 0.2512;
    tmax  = 125+tstart-1;
    yPos = 0.62;
    fooofRange = [64, 140];
else  %for Voytek Data with Fs 1000
    % lambda0 = 0.9;
    tmax = 50+tstart-1;%300;
    yPos = 0.12;
    fooofRange = [75, 150];
end

tdur = 0.5; % time duration for which we want to normalize t;
settings = struct();
settings.max_n_peaks = 1;
settings.peak_width_limits = [4,80];

Fs = round(((tmax/dt)-(tstart-1))/tdur);
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
        end

        for j = 1:nTrials
            if numericalFlag
                [t1,y1(iAlpha,iLambda,j,:)] = fde12(alpha(iAlpha),f1,dt,tmax,[0 0],dt);
            else
                z2 = conv(randn(1,length(t1)),squeeze(z1(iAlpha,iLambda,:)));
                y1(iAlpha,iLambda,j,:) = z2(1:length(t1));
            end
        
        end
        Power1(iAlpha,iLambda,:,:)= mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
        %Power1(iAlpha,iLambda,:) = squeeze(mean(abs(fft(squeeze(y1(iAlpha,iLambda,:,tstart:end))'))'.^2));
     
    end
end
[~,f] = mtspectrumc(squeeze(y1(iAlpha,iLambda,:,tstart:end))',params);
t2 = (t1(tstart:end)-t1(tstart)).*tdur./(t1(end)-t1(tstart));

% df = 1/tdur;
% f = 0:df:(length(t2)-1)*df;
%
% if length(Power1)~=length(f)
%     f = 0:df:(length(t2))*df;
% end

centreFreq = 20:10:120;
minFreqVal = 4; maxFreqVal = 160;
freqRangeWidth = 80;
for jFreq = 1:length(centreFreq)
    freq_range{jFreq} = [centreFreq(jFreq)-freqRangeWidth/2,centreFreq(jFreq)+freqRangeWidth/2];
    if freq_range{jFreq}(1)<minFreqVal
        freq_range{jFreq}(1) = minFreqVal;
    end

    if freq_range{jFreq}(2)>maxFreqVal
        freq_range{jFreq}(2) = maxFreqVal;
    end
end

for iAlpha = length(alpha) %at alpha 1.5
    for iLambda = 1:length(lambda0)
        [ph0(iAlpha,iLambda),PosPh0] = max(log10(Power1(iAlpha,iLambda,1:end-5))); %removing the initial high power:adjust the no. 5 accordingly
        cf00(iAlpha,iLambda) = f(PosPh0(1));
        for jFreq = 1:length(centreFreq)
            fd = fooof(f,squeeze(Power1(iAlpha,iLambda,:)), freq_range{jFreq},settings,true);
            exponent(iAlpha,iLambda,jFreq) = fooof_get_params(fd,squeeze(Power1(iAlpha,iLambda,:)),settings);
        end
    end
end

for iAlpha=1:length(alpha)
    for iLambda = 1:length(lambda0)
        fd = fooof(f,squeeze(Power1(iAlpha,iLambda,:)), fooofRange,settings,true);
        exponent1(iAlpha,iLambda) = fooof_get_params(fd,squeeze(Power1(iAlpha,iLambda,:)),settings);
    end
end

for iLambda = 1:length(lambda0)
    for jFreq = 1:length(centreFreq)
        mdl = fitlm(alpha,squeeze(exponent(:,iLambda,jFreq)));
        slope(iLambda,jFreq) = mdl.Coefficients.Estimate(2);
    end
end

%%
figure
contourf(lambda0,centreFreq,slope',[1.8 2.2],'ShowText','on');hold on;
colormap(bone);
plot(lambda0,cf00(end,:),'.','MarkerSize',20);
xlabel('\lambda');
ylabel('Centre Frequency for Slope (Hz)');
set(gca,'FontWeight','bold','FontSize',14,'LineWidth',2,'box','off')
text(0.5,40,'<1.8','FontSize',14,'FontWeight','bold');
text(0.5,80,'~2','FontSize',14,'FontWeight','bold');
% 


