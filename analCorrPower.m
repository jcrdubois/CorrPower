function power = analCorrPower(rho,alpha,N,tail)
% USAGE: beta = analCorrPower(rho,alpha,N)
% inputs
% rho      = expected correlation
% alpha    = threshold probability for rejecting the null hypothesis 
%            (Type I error rate) 
% N        = sample size
% tail     = 'left', 'right' or 'both' (mirroring Matlab's nomenclature)
%
% 07/01/2016. Julien Dubois, from scratch
% Ref: Hulley SB, Cummings SR, Browner WS, Grady D, Newman TB. 
% Designing clinical research : an epidemiologic approach. 4th ed. 
% Philadelphia, PA: Lippincott Williams & Wilkins; 2013. 
% Appendix 6C, page 79

if nargin<4
    tail   = 'right'; % one-sided test
end

C = 0.5 * log((1+rho)./(1-rho));

switch tail
    case {'left','right'}
        z_alpha = norminv(1-alpha);
    case 'both'
        z_alpha = norminv(1-alpha/2);
end

z_beta = C * sqrt(N-2) - z_alpha;

power = normcdf(z_beta);