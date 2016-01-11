function [rho,power,rate_alloc]	= water_fill_inv(gain,rate);
%[rho,power,rate_alloc]	= water_fill_inv(gain,rate);
%  Computes the power required to achieve a certain RATE for a vector
%  communications channel using the "water filling" solution, based on
%  the power gain vector GAIN, which must be in row vector
%  form. Returns the required SNR (linear, not dB) in RHO, and the
%  required power distribution in POWER.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for xi and q:
%
%      /        RATE       \ 1/q                2
%      |       2           |                 sig
% xi = |-------------------|            xi > ----
%      |  ---     2      2 |                    2
%      |  | |  lam  / sig  |                 lam
%      \   q      q        /  ,                 q ,
%
% so that q is the largest index for which the inequality on the
% right is true.

xi	= ((2^rate)./(cumprod(gain))).^(1./[1:length(gain)]);
q	= max(find(xi>1./gain));
xi	= xi(q);
power	= xi - 1./gain(1:q);
rho	= sum(power);

if(q<length(gain))
    power(length(gain))	= 0;
end
rate_alloc	= log2(1+power);
