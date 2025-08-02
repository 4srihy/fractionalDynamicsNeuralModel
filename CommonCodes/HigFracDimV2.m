function [HFD,rmse,kmaxFinal,minuslnK,lnL] = HigFracDimV2(data, kmax, maxRmse,fastFlag,optimiseRmseFlag,plotFlag)
%fit function changes by Srishty 29-05-23

%fastFlag: if fastFlag, it would not calculate rmse
%optimiseRmseFlag: rmse is optimised ony if it is on
%plotFlag: plots log(L(k) vs Log(1/k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Higuchi Fractal Dimension of a time series %
% The Higuchi Fractal Dimension describes the degree of self-similarity %
% along a 1D time series, according to a sequencing integer, k. The     %
% function returns the dimensionality value, D, where 1 <= D <= 2. A    %
% higher D corresponds to higher complexity. It is a time measure of    %
% complexity.                                                           %    
%                                                                       %
% Higuchi, T. (1988). Approach to an Irregular Time Series on the Basis %
%         of the Fractal Theory. Physica D, 31, 277-283.                %
% Accardo, A., Affinito, M., Carrozzi, M., & Bouquet, F. (1997). Use of %
%         the fractal dimension for the analysis of                     %
%         electroencephalographic time series. Biol Cybern, 77, 339-350.% 
%         doi:10.1007/s004220050394                                     %
%                                                                       %
% INPUTS                                                                %
%  Parameters:                                                          %
%   data = time series                                                  %                                               %
%   kmax = maximum size of iteration                                    %
%                                                                       %
%   Standard value is kmax = 8                                          %
%                                                                       %
% OUTPUTS                                                               %
%   result = Higuchi fractal dimension value                            %
%                                                                       %
%                       Created by Brian Lord                           %
%                       University of Arizona                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Higuchi Fractal Dimension calculation
if ~exist('kmax','var')     kmax=8;             end
if ~exist('maxRmse','var')  maxRmse = 0.05;     end
if ~exist('fastFlag','var') fastFlag = 0;       end
if ~exist("optimiseRmseFlag",'var') optimiseRmseFlag = 1; end
if ~exist('plotFlag','var') plotFlag = 0;       end
% Initialize components
N = length(data);
L = zeros(1,kmax);
x = zeros(1,kmax);
y = zeros(1,kmax);

for k = 1:kmax
    for m = 1:k
        norm_factor = (N-1)/(round((N-m)/k)*k); % normalization factor
        X = sum(abs(diff(data(m:k:N)))); % generate iterative sum difference
        L(m)=X*norm_factor/k; % mean of normalized differences
    end

    y(k)=log(sum(L)/k); % ln(L(k))
    x(k)=log(1/k); % ln(1/k)
end
iRed = 0;   %final kmax = kmax-iRed
if fastFlag
    D = polyfit(x,y,1);
    HFD = D(1);
    kmaxFinal = kmax;
    rmse = nan;
else
    [D,E] = fit(x',y','poly1'); % linear regression fit
    
        if optimiseRmseFlag
        while E.rmse>maxRmse
            iRed = iRed+1;
            [D,E] = fit(x(1:kmax-iRed)',y(1:kmax-iRed)','poly1');
    
            if (kmax-iRed)<3    break;      end   
        end
        end
    
   
    rmse = E.rmse;
    HFD = D.p1; % return slope
    kmaxFinal = kmax-iRed;
end
minuslnK = x(1:kmax-iRed);
lnL = y(1:kmax-iRed);
 if plotFlag plot(minuslnK,lnL,'-*');hold on;   end
end