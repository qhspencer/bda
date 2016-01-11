function [R_1,R_2,max_1,max_2,max_cap] = rateregion2d(n_ru,snr,atten);
% [R_1,R_2,max_1,max_2,max_cap] = rateregion2d(n_ru,snr,atten);
%    generates the rate region for a 2-user MIMO channel with N_RU
%    (default 2) antennas per user, and the specified SNR (default
%    10). The second user's channel is attenuated relative to the
%    first by ATTEN dB.
% [R_1,R_2,max_1,max_2,max_cap] = rateregion2d(H,snr);
%    generates the rate region for an arbitrary H matrix at the
%    specified SNR(s). A 2-user MIMO channel is assumed, so the number
%    of rows must be a multiple of 2.

n_u	= 2;
if(nargin==3)
    n_r	= n_ru*n_u;
    n_t	= n_r;
    rho	= 10^(snr/10);

    %% generate random H
    H	= random('c',[n_r,n_t],1/2);
    H_1	= H(1:n_ru,:);
    H_2	= H(n_ru+1:n_r,:) / 10^(atten/20);
elseif(nargin==2)
    H	= n_ru;
    rho	= 10.^(snr/10);
    [n_r,n_t]	= size(H);
    n_ru	= n_r/n_u;

    H_1	= H(1:n_ru,:);
    H_2	= H(n_ru+1:n_r,:);
end

%%%%%%% initialize arrays
n_points= 500;
n_rows	= 3;
max_cap	= zeros(1,n_rows);
sum_cap	= zeros(1,n_rows);
R_1	= zeros(n_rows,n_points);
R_2	= zeros(n_rows,n_points);
%n_rows	= 4;
rho_vec	= linspace(0,rho,n_points);

%% calculate vectors for H_1 -- no interference
[U,S,V_1]	= svd(H_2);
V_1		= fliplr(V_1);
tmp		= V_1(:,1:n_ru);
[u,s_1,v_1]	= svd(H_1*tmp);
M_1		= tmp*v_1(:,1:n_ru);
lam_H_1		= fliplr(diag(s_1)).'.^2;

%% calculate vectors for H_1 -- w/ interference
[u,s_1,v_1]	= svd(H_1);
M_1_i		= v_1(:,1:n_ru);
lam_H_1_i	= fliplr(diag(s_1)).'.^2;

%% calculate vectors for H_2 -- no interference
[U,S,V_2]	= svd(H_1);
V_2		= fliplr(V_2);
tmp		= V_2(:,1:n_ru);
[u,s_2,v_2]	= svd(H_2*tmp);
M_2		= tmp*v_2(:,1:n_ru);
lam_H_2		= fliplr(diag(s_2)).'.^2;

%% calculate vectors for H_2 -- w/ interference
[u,s_2,v_2]	= svd(H_2);
M_2_i		= v_2(:,1:n_ru);
lam_H_2_i	= fliplr(diag(s_2)).'.^2;

%% calculate rate region
for n = 1:n_points
    R_1(1,n)	= water_fill(lam_H_1,rho_vec(n));
    R_2(1,n)	= water_fill(lam_H_2,rho-rho_vec(n));

    [R_1(2,n),lam_Z_1] = water_fill(lam_H_1_i,rho_vec(n));
    R_i_1	= eye(n_ru) + H_2*M_1_i*diag(lam_Z_1)*M_1_i'*H_2';
    lam_H_2_a	= svd(M_2'*H_2'*inv(R_i_1)*H_2*M_2).';
    R_2(2,n)	= water_fill(lam_H_2_a,rho-rho_vec(n));

    [R_2(3,n),lam_Z_2] = water_fill(lam_H_2_i,rho_vec(n));
    R_i_2	= eye(n_ru) + H_1*M_2_i*diag(lam_Z_2)*M_2_i'*H_1';
    lam_H_1_a	= svd(M_1'*H_1'*inv(R_i_2)*H_1*M_1).';
    R_1(3,n)	= water_fill(lam_H_1_a,rho-rho_vec(n));

    for m = 1:n_rows
	sum_cap(m)	= R_1(m,n)+R_2(m,n);
	if(sum_cap(m)>max_cap(m))
	    max_cap(m)	= sum_cap(m);
	    max_ind(m)	= n;
	end
    end
end

for n	= 1:n_rows
    max_1(n)	= R_1(n,max_ind(n));
    max_2(n)	= R_2(n,max_ind(n));
end

%%% experimental -- I don't think this changes the final result
%%% after all
%
%[lam_H,ind]	= sort([lam_H_1,lam_H_2]);
%lam_H	= fliplr(lam_H);
%ind	= fliplr(ind);
%
%%% below is a slightly modified version of the function water_fill:
%test1	= [1:length(lam_H)]./lam_H;
%test2	= cumsum(1./lam_H) + rho;
%q	= max(find(test1<test2));
%if(isempty(q))
%    rate	= 0;
%    lam_Z	= zeros(size(lam_H));
%else
%    X	= test2(q)/q;
%    lam_Z	= max(X - 1./lam_H,0);
%    rate	= log2(1+lam_Z.*lam_H);
%end
%%% end of water_fill stuff
%
%[x,ind2]	= sort(ind);
%max_1(1)	= sum(rate(ind2([1 2])));
%max_2(1)	= sum(rate(ind2([3 4])));

if(nargout == 0)
    plot(R_1',R_2',max_1,max_2,'6*');
    ylabel('Capacity for user 2');
    xlabel('Capacity for user 1');
    title(['Rate Regions for an SNR of ',int2str(snr),' dB']);
    legend('Rate Region BD','Rate Region U1',...
	   'Rate Region U2','Max. Sum Rates',3);
end
