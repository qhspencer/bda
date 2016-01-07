function [x,y]	=	cdf(in,num,style);
%   function [x,y]	=	cdf(in,num,style);
%
%	Plots a sample CDF of the elements in the vector IN. The
%	variable NUM is the number of histogram bins (default is 100),
%	and STYLE is the linestyle for the plot.
%

if nargin<3	style	= '-';		end
if nargin<2	num	= 100;		end

in	= in(:);
[b,a]	= hist(in,num);

x	= a;
y	= cumsum(b)/length(in);
if nargout<2
    plot(x,y,style);
%    axis([min(x) max(x) 0 1]);
end