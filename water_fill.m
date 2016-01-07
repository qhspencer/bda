function [rate,power]	= water_fill(gain,rho);
%[rate,power]	= water_fill(gain,rho)
%  Computes the "water filling" solution for a vector communications
%  channel with power (not amplitude) gain values in GAIN and an SNR
%  of RHO (linear, not dB). Returns the achievable RATE (bits/sec/Hz),
%  and the power distribution vector POWER required to get it.

test1	= [1:length(gain)]./gain;
test2	= cumsum(1./gain) + rho;
q	= max(find(test1<test2));
if(isempty(q))
    rate	= 0;
    power	= zeros(size(gain));
else
    X	= test2(q)/q;
    power	= max(X - 1./gain,0);
    rate	= sum(log2(1+power.*gain));
end

if(nargout>1 & q<length(gain))
    power(length(gain))	= 0;
end
