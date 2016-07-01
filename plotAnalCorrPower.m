function power = plotAnalCorrPower(rho)
% USAGE: power = plotAnalCorrPower(rho)
% output
%   power [length(alphas) length(tails) length(Ns)]
% input 
%   rho expected correlation
% make a plot of power vs. sample size
% using the analytical power calculation for a correlation analysis
% show results for different criteria and tails
% 07/01/2016 Julien Dubois from scratch
dbstop if error

% expected correlation
if nargin<1
    rho = 0.2;
end

% criterion for significance
alphas = [0.05 0.01 0.001];
c = {'k','m','r'}; % corresponding colors for plotting
% one-sided or two-sided?
tails  = {'right','both'};
s = {'o-','o--'}; % corresponding line styles for plotting
% sample sizes
Ns = [5:50 55:5:100 110:10:200 220:20:500 550:50:1000];

% DO THE WORK
fprintf('Computing power at all sample sizes\n');
tic
power = zeros(length(alphas),length(tails),length(Ns));
for itail = 1:length(tails)
    tail = tails{itail};
    fprintf('\t\t tail = %s\n',tail);
    for ialpha = 1:length(alphas)
        alpha = alphas(ialpha);
        fprintf('\t alpha = %0.3f\n',alpha);
        for iN = 1:length(Ns),
            N = Ns(iN);
            power(ialpha,itail,iN) = analCorrPower(rho,alpha,N,tail);
        end
    end
end
elapsed = toc;
fprintf('done in %.1fs\n',elapsed);

% plot result
figure;hold on;
cell4legend = cell(1,1);
for itail = 1:length(tails),
    for ialpha = 1:length(alphas)
        plot(Ns,squeeze(power(ialpha,itail,:)),[c{ialpha},s{itail}]);
        cell4legend{(itail-1)*length(alphas)+ialpha} = sprintf('\\alpha = %.3f (%s)',alphas(ialpha),tails{itail});
    end
end
legend(cell4legend,'Location','SouthEast');
ylabel('Statistical Power (analytical)');
xlabel('Sample size');
title(sprintf('Expected effect size: \\rho = %.3f',rho)); 

keyboard
 
