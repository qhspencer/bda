function [x,y] = succ_opt_ord(n_t,n_ru,ratepoint,power,params);
% [x,y] = succ_opt_ord(n_t,n_ru,ratepoint,power,params);
%
% Tests a proposed algorithm for ordering the selection of user
% matrices for successive optimization.

if(nargin<1)	n_t	= 4;		end
if(nargin<2)	n_ru	= [2 2];	end
if(nargin<3)	ratepoint = [3 3];	end
if(nargin<4)	power	= [1 1];	end
if(nargin<5)	params	= {2000, 1};	end
n_u	= length(n_ru);
n_r	= sum(n_ru);
if(ratepoint==0) % random ratepoint
    randrate	= 1;
else
    randrate	= 0;
    if(length(n_ru)~=length(ratepoint))
	error('n_ru and ratepoint must be the same size.');
    end
end

numTrials	= params{1};
perm_vecs	= perms([1:n_u]);
num_perms	= size(perm_vecs,1);
ii		= cumsum(n_ru);

if(power==0)	% random power
    randpower	= 1;
else
    randpower	= 0;
    for m	= 1:n_u	% Create power scaling matrix Pm
	Pm(ii(m)-n_ru(m)+1:ii(m),:)	= ones(n_ru(m),n_t)*power(m);
    end
end

% preallocate memory
P	= zeros(5,numTrials);

for n = 1:numTrials
    if(randrate)	%% random rate, fixed to range: [2,10]
	ratepoint = random('u',[1 n_u])*8+2;
    end
    if(randpower)
	power	= random('u',[1 n_u])*1.5+0.5;
	for m	= 1:n_u	% Power Multiplier
	    Pm(ii(m)-n_ru(m)+1:ii(m),:)	= ones(n_ru(m),n_t)*power(m);
	end
    end
    H	= random('c',[n_r,n_t],1/2).*Pm;

    %% Block-Diagonalization
    P(1,n)	= sum(block_diag(H,n_ru,ratepoint));

    %% Sucessive Optimization -- Find optimal choice
    P_so	= zeros(1,num_perms);
    for m	= 1:num_perms
	P_so(m)	= sum(succ_opt(H,n_ru,ratepoint,perm_vecs(m,:)));
    end
    P(2,n)	= min(P_so);
    order_opt	= perm_vecs(find(P_so==P(2,n)),:);

    %% Sucessive Optimization -- New ordering algorithm
    for m	= 1:n_u
	H_m	= H(ii(m)-n_ru(m)+1:ii(m),:);
	H_m_	= H([1:ii(m)-n_ru,ii(m)+1:ii(n_u)],:);

	%% Angle calculation -- the best so far
	[QA,x]	= qr(H_m',0);
	[QB,x]	= qr(H_m_',0);
	s	= svd(QA'*QB);
	alpha(m)= acos(min(min(s),1))*n_ru(m);

	%% Option 2 -- the frobenius norm (cheap, but not so effective).
	alpha2(m)= norm(H_m,'fro');
    end
    [x,order]	= sort(alpha);
    P(3,n)	= sum(succ_opt(H,n_ru,ratepoint,order));
    [x,order2]	= sort(alpha2);
    P(4,n)	= sum(succ_opt(H,n_ru,ratepoint,order2));

    %% Randomly chosen order -- for comparison
    P(5,n)	= P_so(ceil(rand*size(perm_vecs,1)));

    waitbar(n/numTrials);
end
clear x

if(params{2})
    cdfstr	= 'ccdf';
    keypos	= 'top';
    rel		= '$>$';
else
    cdfstr	= 'cdf';
    keypos	= 'bottom';
    rel		= '$<$';
end
for n = 1:size(P,1)
    eval(['[x(n,:),y(n,:)] = ',cdfstr,'(10*log10(P(n,:)),200);']);
end
if(nargout == 0)
    plot(x',y');
    legend('BD','Best SO','SO Alg','SO Exp','Random',1);
    gset(['key right ',keypos,' Right noreverse spacing 2']);
    title(['Random H, ',int2str(n_u),' users']);
    ylabel(['Probability SNR ',rel,' Abscissa']);
    xlabel('Required Average SNR (dB)');
end
