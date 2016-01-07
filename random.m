function out	= random(type,siz,var);
%	function out	= random(type,size,var);
%
%	Generates random variables of the specified type and size
%	var specifies variance for all but uniform and exponential.
%
%	Supported Distributions:
%
%	'uniform'	uniform
%	'exp'		exponential
%	'gauss'		gaussian
%	'compgauss'	complex gaussian
%	'rayleigh'	rayleigh
%	'laplace'	laplacian

%	last modified April 30, 1998
%	by Quentin Spencer

if nargin==1
    var	= 1;
    siz	= [1 1];
elseif nargin==2
    if ~all(size(siz)==[1 2])
	var	= siz;
	siz	= size(var);
    else
	var	= 1;
    end
end
if all(size(var)~=siz)
    var	= var*ones(siz);
end


if strcmp(type,'uniform') | strcmp(type,'u')
    out	= rand(siz);

elseif strcmp(type,'chi')
    if(var<2)
	error('Degrees of freedom must be >= 1');
    end
    out	= sum((sqrt(-2*log(rand(var/2,siz(2)))).*...
		  cos(rand(var/2,siz(2))*2*pi)).^2);

elseif strcmp(type,'exp') | strcmp(type,'e')
    out	= -log(rand(siz)).*sqrt(var);

elseif strcmp(type,'laplace') | strcmp(type,'l')
    out	= round(rand(siz))*2-1;
    out	= -log(rand(siz)).*sqrt(var/2) .* out;

elseif strcmp(type,'rayleigh') | strcmp(type,'r')
    out	= sqrt(-2*var.*log(rand(siz)));
%   var_out will be  (2 - pi/2)*var_in

elseif strcmp(type,'gauss') | strcmp(type,'g')
    out	= sqrt(-2*var.*log(rand(siz))).*cos(rand(siz)*2*pi);

elseif strcmp(type,'compgauss') | strcmp(type,'c')
    out	= sqrt(-2*var.*log(rand(siz))).*exp(j*rand(siz)*2*pi);
end
