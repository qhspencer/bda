function [s,i] = sortinv(X);
% [s,i] = sortinv(X);
% Sorts the elements of X in decending order, and returns S,
% equivalent to fliplr(sort(X)). The vector I 

[S,I]	= sort(X);
s	= fliplr(S);
if(nargout>1)
    II	= fliplr(I);
    for n = 1:length(X)
	i(II(n))	= n;
    end
end
