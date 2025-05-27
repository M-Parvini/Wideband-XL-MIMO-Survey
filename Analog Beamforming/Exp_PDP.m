function [pathGains, pathDelays] = Exp_PDP(BW, delay_spread, L)

%%% Exponential power delay profile
Tm = (1/BW)/delay_spread;
pathDelays = (1/BW)*([0:L-1]);
% PDP_EXP_norm = exp((-1*Tm)*([0:params.sys.L-1]));
PDP_EXP = Tm*exp((-1*Tm)*([0:L-1]));
PDP_EXP_norm = reNormalize(PDP_EXP);
pathGains=pow2db(PDP_EXP_norm);

end

function y = reNormalize(x)
% Renormalization to 1
Powersum = sum(x);
y = (1/Powersum)*x;
% y = pow2db(y);
end