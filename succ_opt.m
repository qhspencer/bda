function [P,M,lam_H] = succ_opt(H,n_ru,ratepoint,order);
% [P, M, lam_H] = succ_opt(H,n_ru,ratepoint,order);
%
% Performs successive optimization of H where N_RU is a vector
% containing the number of receivers for each user and RATEPOINT is
% the corresponding desired rate (bits/sec/Hz). Optimization is done
% in the order specified by the ORDER vector. ORDER, N_RU and
% RATEPOINT must be the same size, and sum(n_ru) must be the number of
% rows in H. Returns the power required to achieve the rate point (in
% vector form), and the corresponding singular vectors and singular
% values.

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

h	= [];
for m = 1:n_u
    h	= [h; H(ii(order(m))-n_ru(order(m))+1:ii(order(m)),:)];
end
H	= h;
n_ru	= n_ru(order);
ii	= cumsum(n_ru);

h_m_	= [];	h_m	= [];
for m = 1:n_u
    N	= ii(m)-n_ru(m)+1;
    h_m_= [h_m_; h_m];
    h_m	= H(N:ii(m),:);

    if(m==1)
	[u,s,v]	= svd(h_m);
	M	= v(:,1:ii(1));
	lam_H_	= diag(s).';
	lam_H	= lam_H_(1:ii(1));
    else
	[U,S,V]	= svd(h_m_);
	V	= fliplr(V);
	tmp	= V(:,1:n_t-ii(m-1));
	R	= eye(n_ru(m))+h_m*M*diag(lam_Z)*M'*h_m';
	[u,s,v]	= svd(tmp'*h_m'*inv(R)*h_m*tmp);
	M(:,N:ii(m))	= tmp*v(:,1:n_ru(m));
	lam_H_	= diag(s).';
	lam_H(N:ii(m))	= lam_H_(1:n_ru(m));
    end
    [P_m(m),lam_Z_m]	= water_fill_inv(lam_H(N:ii(m)),ratepoint(m));
    lam_Z(N:ii(m))	= lam_Z_m;
end

P	= P_m;

