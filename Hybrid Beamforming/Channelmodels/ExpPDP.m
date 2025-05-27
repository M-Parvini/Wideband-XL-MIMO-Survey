function pathGains = ExpPDP(Fs, Tm, L)
%%%
% Exponential power delay profile
%%%
Tm = (1/Fs)/Tm;
pathDelay = (1/Fs)*([0:L-1]);
% PDP_EXP = Tm*exp((-1*Tm)*([0:L-1]));
PDP_EXP = exp((-1*Tm)*([0:L-1]));
% PDP_EXP_norm = reNormalize(PDP_EXP);
pathGains=10*log10(PDP_EXP);

end