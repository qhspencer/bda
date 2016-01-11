function [P,M_bd,lam_H,lam_Z,D_bd] = block_diag(H,n_ru,ratepoint);
% [P,M_bd,lam_H,lam_Z,D_bd] = block_diag(H,n_ru,ratepoint);
%
% Performs block-diagonalization of H where N_RU is a vector
% containing the number of receivers for each user and RATEPOINT is
% the corresponding desired rate (bits/sec/Hz). N_RU and RATEPOINT
% must be the same size, and sum(n_ru) must be the number of rows in
% H. Returns the power P required to achieve the rate point (in vector
% form), and the corresponding singular vectors (M_BD) and singular
% values (LAM_H). Other optional outputs are the power loading
% coefficients (LAM_Z) and the demodulation vectors D_BD.

n_u	= length(n_ru);
if(length(ratepoint) ~= n_u)
    error('ratepoint and n_ru must be the same size');
end
[n_r,n_t]	= size(H);
if(n_r > n_t)
    error('illegal matrix size -- n_R > n_T');
end
ii	= cumsum(n_ru);
if(ii(n_u) ~= n_r)
    error('n_ru specification does not match H');
end

for m = 1:n_u
    H_m		= H(ii(m)-n_ru(m)+1:ii(m),:);
    H_m_	= H([1:ii(m)-n_ru(m),ii(m)+1:ii(n_u)],:);
    H_m_null	= null(H_m_);
    [U,S,V]	= svd(H_m*H_m_null);
    lam_H_m	= (diag(S(1:n_ru(m),1:n_ru(m))).^2)';
    [P_m(m),lZ]	= water_fill_inv(lam_H_m,ratepoint(m));

    %%% Create matrices of the relevant singular values and vectors
    %%% (optional output).
    M_bd(:,ii(m)-n_ru(m)+1:ii(m))	= H_m_null*V(:,1:n_ru(m));
    lam_H(ii(m)-n_ru(m)+1:ii(m))	= lam_H_m;
    D_bd{m}				= U;
    lam_Z(ii(m)-n_ru(m)+1:ii(m))	= lZ;
end

P	= P_m;
