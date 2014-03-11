function y=movingvar(X,N)
% y=movingvar(X,N)
% Calculates N-point moving variance of  Vector X
% Highly recommend that N be odd (no error checking)
% Note: first and last N/2 points will be unreliable.
% Output will be a column vector.
% Please don't distribute without header info
% Scott Seidman  1/23/99

d

X=X(:);
XSQR=X.*X;
convsig=ones(1,N);
y=(conv(convsig,XSQR)-(conv(convsig,X).^2)/N)/(N-1);

y=y(ceil(N/2):length(X)+floor(N/2));