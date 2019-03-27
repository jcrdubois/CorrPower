function power = plotAnalCorrPower(r,targetPower,alphas,tail,Nmax)
% USAGE: power = plotAnalCorrPower(r,targetPower,alphas,tail,Nmax)
% output
%   power [length(alphas),length(Ns)]
% input 
%   r: expected correlation (effect size)
%   targetPower: how much power you want your experiment to have
%   alphas: significance thresholds of interest
%   tail: 'right' or 'both", depends on your hypoethesis
%   Nmax: until what sample size to compute the iwer

% make a plot of power vs. sample size
% using the analytical power calculation for a correlation analysis
% show results for different criteria and tails
% 07/01/2016 Julien Dubois from scratch
% 03/27/2019 Julien Dubois -- added input arguments, and plotting/output of
% targetN
dbstop if error

% expected correlation
if nargin<1
    r = 0.4;
end
if nargin<2
    targetPower = .99;
end
if nargin<3
    % criterion for significance
    alphas = [0.05 0.01 0.001];
end
if nargin<4
    tail = 'right'; % may be both
end
if nargin<5
    % sample sizes
    Nmax = 200;
end
c = {'k','m','r'}; % corresponding colors for plotting

Ns = 5:Nmax;
% DO THE WORK
fprintf('Computing power at all sample sizes\n');
tic
power = zeros(length(alphas),length(Ns));
fprintf('\t\t tail = %s\n',tail);
for ialpha = 1:length(alphas)
    alpha = alphas(ialpha);
    fprintf('\t alpha = %0.3f\n',alpha);
    power(ialpha,:) = analCorrPower(r,alpha,Ns,tail);
end
elapsed = toc;
fprintf('done in %.1fs\n',elapsed);

% plot result
figure;hold on;
cell4legend = cell(1,1);
for ialpha = 1:length(alphas)
    plot(Ns,squeeze(power(ialpha,:)),c{ialpha});
    cell4legend{ialpha} = sprintf('\\alpha = %.3f',alphas(ialpha));
end
legend(cell4legend,'Location','SouthEast');
ylabel('Statistical Power (analytical)');
xlabel('Sample size');
title(sprintf('Expected effect size: r = %.3f, tail = %s',r, tail)); 

targetN = zeros(1,length(alphas));
for ialpha=1:length(alphas)
    [~,indmini] = min(abs(power(ialpha,:)-targetPower));
    targetN(ialpha) = Ns(indmini);
    fprintf('alpha=%.3f: N for power=%.2f is %d\n',alphas(ialpha),targetPower,targetN(ialpha));
    plot(Ns(indmini)*ones(1,2),[0,1],[c{ialpha},':']);
end
plot([Ns(1),Ns(end)],targetPower*ones(1,2),'k:');



