function [x,y] = bd_capacity(n_t,n_ru,snr,params);
% function [x,y] = bd_capacity(n_t,n_ru,snr,params);
%    returns the capacity for a Joint Transmission (JT) channel,
%    with N_T transmitters, and a vector N_RU specifying the number
%    of receivers per user, at the specified SNR.
%    PARAMS is a cell array containing the following items
%	1 - number of monte-carlo trials	[ 2000 ]
%	2 - 0 for CDF, 1 for CCDF		[  1   ]
%	3 - correlation of H members		[  0   ]
%	4 - single-beam mode			['1beam'], 'all'

if(nargin<1)		n_t	= 4;		end
if(nargin<2)		n_ru	= [2 2];	end
if(nargin<3)		snr	= 20;		end
if(nargin<4)		params	= {2000, 1, 0};	end
if(length(params)<2)	params{2}	= 1;	end
if(length(params)<3)	params{3}	= 0;	end
if(sum(n_ru)>n_t)
    if(length(params)<4)	params{4}	= '1beam';	end
end
n_u	= length(n_ru);
ii	= cumsum(n_ru);
n_r	= ii(n_u);
rho	= 10^(snr/10);
n_trials= params{1};
corr	= params{3};
if(corr)
    %%% R_r stands for either "row" or "receiver" covariance
    R_r	= toeplitz(corr.^[0:n_ru(1)-1]);
    for m = 2:n_u
	block	= zeros(size(R_r,1),n_ru(m));
	R_r	= [R_r, block; block', toeplitz(corr.^[0:n_ru(m)-1])];
    end
    [U,S,V] = svd(R_r);
    S_r	= V*S^(1/2)*V';
else
    S_r	= eye(n_r);
end

for n = 1:n_trials
    H	= S_r*random('c',[n_r,n_t],1/2);

    if(n_r>n_t)
	C(1,n)	= bd_lg(H,n_ru,snr,{params{4},'fast'});
	if(strcmp(params{4},'1beam'))
	    C(2,n)	= bd_lg(H,n_ru,snr,{params{4},'opt'});
	end
    else
	if(n_u>1)
	    for m = 1:n_u
		m_ru	= n_ru(m);
		h_m	= H(ii(m)-m_ru+1:ii(m),:);
		h_m_	= H([1:ii(m)-m_ru,ii(m)+1:ii(n_u)],:);
		[U,S,V]	= svd(h_m_);
		V	= fliplr(V);
		tmp	= V(:,1:m_ru+n_t-n_r);
		[u,s,v]	= svd(h_m*tmp);
		M_wf(:,m_ru*(m-1)+1:m_ru*m) = tmp*v(:,1:m_ru);
		lam_H(m_ru*(m-1)+1:m_ru*m) = (diag(s(1:m_ru,1:m_ru)).^2)';

		%%% blind transmitter stuff
		C_b(m)	= log2(real(det(eye(m_ru)+rho/n_t*h_m*h_m')))/n_u;

		%%% TDMA
		C_tdma(m)	= water_fill((svd(h_m).^2).',rho)/n_u;
	    end
	else
	    [u,s,v]	= svd(H);
	    M_p		= M_wf	= v;
	    lam_H	= diag(s)'.^2;
	    C_b		= log2(real(det(eye(n_r)+rho/n_t*H*H')));
	    C_tdma	= water_fill((svd(H).^2).',rho);
%	    if(enable_experimental_stuff)
%		M_p2	= M_p;
%	    end
	end

	%%% C(1,:) is capacity for channel inversion 
	C(1,n)	= n_r*log2(1+rho/sum(svd(H).^(-2)));
	%%% C(2,:) is capacity for water filling (lower bound)
	C(2,n)	= water_fill(sortinv(lam_H),rho);
	%%% C(3,:) is capacity for equal power transmission
	C(3,n)	= log2(real(det(eye(n_r)+rho/size(M_wf,2)*H*M_wf*M_wf'*H')));
	%%% C(4,:) is the "blind" transmitter
	C(4,n)	= sum(C_b);
	%%% C(5,:) is TDMA with channel knowledge
%	C(5,n)	= sum(C_tdma);
    end

    waitbar(n/n_trials);
end

if(params{2})
    cdfstr	= 'ccdf';
    keypos	= 'bottom';
    rel		= '$>$';
else
    cdfstr	= 'cdf';
    keypos	= 'top';
    rel		= '$<$';
end
for n = 1:size(C,1)
    eval(['[x(n,:),y(n,:)] = ',cdfstr,'(C(n,:));']);
end
