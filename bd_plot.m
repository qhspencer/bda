function [result] = bd_plot(N,figpath,color);
%bd_plots(N,figpath,color);

if(nargin<3)	color	= 1;			end
if(nargin<2)	figpath = 'block_diag/';	end
if(strcmp(figpath,'test')) test	= 1;	else	test	= 0;	end

% channel strings
str1	= 'Inversion';		str2	= 'Block-Diag.';
str3	= 'PCI';		str4	= 'Blind Tx';
str5	= '1 User';		str6	= 'TDMA';

ii	= [0:90]*pi/45;
a	= sin(ii);
b	= cos(ii);

n_trials= 5000;
result	= 1;

if(length(N)>1)
    for n	= 1:length(N)
	bd_plots(N(n),figpath,color);
    end
end

switch(N)
 case 0,
  % generate all plots
  q	= 0;
  while(result)
      q		= q+1;
      result	= bd_plots(q,figpath,color);
  end

 case 1,
  figname	= 'cap_cdf';
  [x_1,y_1]	= bd_capacity(4,[1 1 1 1],10);
  [x_2,y_2]	= bd_capacity(4,[2 2],10);
  [X_1user,Y_1user]	= bd_capacity(4,4,10);
  circ1		= [a*1.0+6.6;	a*0.5+3.45];
  circ2		= [b*0.02+0.8;	b*0.02+0.8];
  X	= [x_1(2,:); x_1(1,:); x_1(3,:); x_1(4,:); X_1user];
  Y	= [y_1(2,:); y_1(1,:); y_1(3,:); y_1(4,:); Y_1user];
  X1	= [x_2(2,:)];	Y1	= [y_2(2,:)];
  X3	= [x_2(3,:)];	Y3	= [y_2(3,:)];
  X4	= [x_2(4,:)];	Y4	= [y_2(4,:)];
  plot(X',Y',X1',Y1','1',X3',Y3','3',X4',Y4','4',circ1',circ2','1');
  legend(str1,str2,str4,str5,3);
  text();
  text(1.8,0.65,'{1,1,1,1}\times 4','HorizontalAlignment','center');
  text(1.8,0.4,'{2,2}\times 4','HorizontalAlignment','center');
  title('Capacity CCDFS for n_T=4 at 10 dB SNR');
  ylabel('Probability Capacity > Abscissa');
  xlabel('Capacity (bits/use)');
  axis('normal');

 case 2,
  figname	= 'cap_SNR';
  outage = 0.1;
  SNR	 = [0:2:20];
  for n = 1:length(SNR)
      disp(sprintf('SNR = %i\n',SNR(n)));
      [x,y]	= bd_capacity(4,[1,1,1,1],SNR(n));
      for m = 1:size(x,1)
	  C_1(m,n)	= x(m,max(find(y(m,:)>=(1-outage))));
      end
      [x,y]	= bd_capacity(4,[2,2],SNR(n));
      for m = 1:size(x,1)
	  C_2(m,n)	= x(m,max(find(y(m,:)>=(1-outage))));
      end
      [x,y]	= bd_capacity(4,4,SNR(n));
      C_1user(n)= x(2,max(find(y(2,:)>=(1-outage))));
  end
  circ1		= [a*0.25+14;	a*0.25+7];
  circ2		= [b*0.7+4.25;	b*0.8+4.4];
  C	= [C_1(2,:); C_1(1,:); C_1(3,:); C_1(4,:); C_1user];
  C1	= [C_2(2,:)];
  C3	= [C_2(3,:)];
  C4	= [C_2(4,:)];
  plot(SNR,C,SNR,C1','1',SNR,C3','3',SNR,C4','4',circ1',circ2','1');
  legend(str2,str1,str3,str4,str5,3);
  text();
  text(15,1,'{1,1,1,1}\times 4','HorizontalAlignment','center');
  text(7,10,'{2,2}\times 4','HorizontalAlignment','center');
  title(['Capacity as a function of SNR, $n_T=4$, $n_R=4$']);
  ylabel(['Capacity at ',num2str(outage),' outage probability']);
  xlabel('SNR (dB)');
  axis('normal');

 case 3,
  figname	= 'cap_n_T';
  outage = 0.1;
  sizes  = [4:8];
  for n = 1:length(sizes)
      disp(sprintf('size = %i\n',sizes(n)));
      [x,y]	= bd_capacity(sizes(n),[1,1,1,1],10);
      for m = 1:size(x,1)
	  C_1(m,n)	= x(m,max(find(y(m,:)>=(1-outage))));
      end
      [x,y]	= bd_capacity(sizes(n),[2 2],10);
      for m = 1:size(x,1)
	  C_2(m,n)	= x(m,max(find(y(m,:)>=(1-outage))));
      end
      [x,y]	= bd_capacity(sizes(n),4,10);
      C_1user(n)= x(2,max(find(y(2,:)>=(1-outage))));
  end
  circ1		= [a*0.05+4.275;	a*0.05+4.2];
  circ2		= [b*0.5+3.2;	b*1.0+6.2];
  C	= [C_1(2,:); C_1(1,:); C_1(3,:); C_1(4,:); C_1user];
  C1	= [C_2(2,:)];
  C3	= [C_2(3,:)];
  C4	= [C_2(4,:)];
  plot(sizes,C',sizes,C1,'1',sizes,C3,'3',sizes,C4,'4',circ1',circ2','1');
  legend(str2,str1,str3,str4,str5,3);
  text();
  text(4.5,1,'{1,1,1,1}\times 4','HorizontalAlignment','left');
  text(4.9,3.5,'{2,2}\times 4','HorizontalAlignment','left');
  title(['Capacity as a function of n_T at 10 dB SNR']);
  ylabel(['Capacity at ',num2str(outage),' outage probability']);
  xlabel('Number of transmit antennas (n_T)');
  axis('normal');

 case 4,
  figname	= 'cap_corr';
  outage = 0.1;
  corr  = [0, 0.25, 0.5, 0.75, 0.85, 0.9, 0.95, 1.00];
  [x,y]	= bd_capacity(4,[1,1,1,1],10,{n_trials*2,1});
  for m = 1:size(x,1)
      C_1(m)	= x(m,max(find(y(m,:)>=(1-outage))));
  end
  for n = 1:length(corr)
      disp(sprintf('corr = %1.2f\n',corr(n)));
      [x,y]	= bd_capacity(4,[2 2],10,{n_trials,1,corr(n)});
      for m = 1:size(x,1)
	  C_2(m,n)	= x(m,max([1, find(y(m,:)>=(1-outage))]));
      end
      [x,y]	= bd_capacity(4,4,10,{n_trials,1,corr(n)});
      C_1user(n)= x(2,max(find(y(2,:)>=(1-outage))));
  end
  circ1		= [a*0.015+0.1;	a*0.015+0.95;	a*0.01+0.7];
  circ2		= [b*0.9+6;	b*1.0+2.0;	b*0.2+0.56];
  C_1	= C_1([2,1,3,4]);
  C1	= [C_2(2,:)];
  C2	= [C_2(1,:)];
  C3	= [C_2(3,:)];
  C4	= [C_2(4,:)];
  plot([0,1],[C_1; C_1]',corr,C_1user,corr,C1,'1',corr,C2,'2',corr, ...
       C3,'3',corr,C4,'4',circ1',circ2','1');
  legend(str2,str1,str3,str4,str5,3);
  axis([0 1 0 10]);
  text();
  text(0.925,3.25,'{1,1,1,1}\times 4','HorizontalAlignment','right');
  text(0.125,4,'{2,2}\times 4','HorizontalAlignment','left');
  title(['Capacity as a function of channel correlation at 10 dB SNR']);
  ylabel(['Capacity at ',num2str(outage),' outage probability']);
  xlabel('Rx antenna element correlation');
  axis('normal');

 case 5,
  figname	= 'rate_region';
  snr	= [3, 10, 20];
  H	= random('c',[4,4],1/2);
  maxtmp_1	= [];	maxtmp_2	= [];
  for n = 1:length(snr)
      [R_1,R_2,max_1,max_2]	= rateregion2d(H,snr(n));
      ind	= int2str(n);
      eval(['BD_1(',ind,',:)	= [R_1(1,:)];']);
      eval(['BD_2(',ind,',:)	= [R_2(1,:)];']);
      eval(['U1_1(',ind,',:)	= [R_1(2,:)];']);
      eval(['U1_2(',ind,',:)	= [R_2(2,:)];']);
      eval(['U2_1(',ind,',:)	= [R_1(3,:)];']);
      eval(['U2_2(',ind,',:)	= [R_2(3,:)];']);
      maxtmp_1	= [maxtmp_1, max_1];
      maxtmp_2	= [maxtmp_2, max_2];
      rho	= 10^(snr(n)/10);
      c_1	= log2(real(det(eye(2) + rho/4*H(1:2,:)*H(1:2,:)')));
      c_2	= log2(real(det(eye(2) + rho/4*H(3:4,:)*H(3:4,:)')));
      Bl_1(n,:)	= [0,c_1];
      Bl_2(n,:)	= [c_2,0];
  end
  max_1	= maxtmp_1;	max_2	= maxtmp_2;
  n	= length(snr);
  plot(BD_1(1,:)',BD_2(1,:)','1;Rate Region BD;',...
       BD_1(2:n,:)',BD_2(2:n,:)','1;;',...
       U1_1(1,:)',U1_2(1,:)','2;Rate Region U1;',...
       U1_1(2:n,:)',U1_2(2:n,:)','2;;',...
       U2_1(1,:)',U2_2(1,:)','3;Rate Region U2;',...
       U2_1(2:n,:)',U2_2(2:n,:)','3;;',...
       Bl_1(1,:)',Bl_2(1,:)','4;Blind Transmitter;',...
       Bl_1(2:n,:)',Bl_2(2:n,:)','4;;',...
       max_1,max_2,'5*;Max. Sums;');
  ylabel('Capacity for user 2');
  xlabel('Capacity for user 1');
  title(['Rate Regions at SNRs of 3, 10, and 20 dB']);
  axis('square');

 case 6,
  figname	= 'rate_region_atten';
  snr	= [3, 10, 20];
  atten	= 10;
  H	= random('c',[4,4],1/2);
  H(3:4,:)	= H(3:4,:)/sqrt(10^(atten/10));	% attenuate last two rows
  maxtmp_1	= [];	maxtmp_2	= [];
  for n = 1:length(snr)
      [R_1,R_2,max_1,max_2]	= rateregion2d(H,snr(n));
      ind	= int2str(n);
      eval(['BD_1(',ind,',:)	= [R_1(1,:)];']);
      eval(['BD_2(',ind,',:)	= [R_2(1,:)];']);
      eval(['U1_1(',ind,',:)	= [R_1(2,:)];']);
      eval(['U1_2(',ind,',:)	= [R_2(2,:)];']);
      eval(['U2_1(',ind,',:)	= [R_1(3,:)];']);
      eval(['U2_2(',ind,',:)	= [R_2(3,:)];']);
      maxtmp_1	= [maxtmp_1, max_1];
      maxtmp_2	= [maxtmp_2, max_2];
      rho	= 10^(snr(n)/10);
      c_1	= log2(real(det(eye(2) + rho/4*H(1:2,:)*H(1:2,:)')));
      c_2	= log2(real(det(eye(2) + rho/4*H(3:4,:)*H(3:4,:)')));
      Bl_1(n,:)	= [0,c_1];
      Bl_2(n,:)	= [c_2,0];
  end
  max_1	= maxtmp_1;	max_2	= maxtmp_2;
  n	= length(snr);
  plot(BD_1(1,:)',BD_2(1,:)','1;Rate Region BD;',...
       BD_1(2:n,:)',BD_2(2:n,:)','1;;',...
       U1_1(1,:)',U1_2(1,:)','2;Rate Region U1;',...
       U1_1(2:n,:)',U1_2(2:n,:)','2;;',...
       U2_1(1,:)',U2_2(1,:)','3;Rate Region U2;',...
       U2_1(2:n,:)',U2_2(2:n,:)','3;;',...
       Bl_1(1,:)',Bl_2(1,:)','4;Blind Transmitter;',...
       Bl_1(2:n,:)',Bl_2(2:n,:)','4;;',...
       max_1,max_2,'5*;Max. Sums;');
  ylabel('Capacity for user 2');
  xlabel('Capacity for user 1');
  title(['Rate Regions at SNRs of 3, 10, and 20 dB']);
  axis('normal');

 case 7,
  figname	= 'so_222';
  [x,y]	= succ_opt_ord(6,[2 2 2],0,0,{n_trials,1});
  plot(x',y');
  legend('BD','SO-Optimal','SO-Angle Algorithm',...
	 'SO-Frobenius Norm','SO-Random',1);
  title('(a) n_{R_j}=(2,2,2)');
  ylabel(['Probability SNR > Abscissa']);
  xlabel('Required Average SNR (dB)');
  axis([0 50 0 1]);

 case 8,
  figname	= 'so_123';
  [x,y]	= succ_opt_ord(6,[1 2 3],0,0,{n_trials,1});
  plot(x',y');
  legend('BD','SO-Optimal','SO-Angle Algorithm',...
	 'SO-Frobenius Norm','SO-Random',1);
  title('(b) n_{R_j}=(1,2,3)');
  ylabel(['Probability SNR > Abscissa']);
  xlabel('Required Average SNR (dB)');
  axis([0 50 0 1]);

 otherwise,
  result	= 0;
end
