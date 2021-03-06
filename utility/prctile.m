function y = prctile(x,p);
%PRCTILE gives the percentiles of the sample in X.
%   Y = PRCTILE(X,P) returns a value that is greater than P percent
%   of the values in X. For example, if P = 50  Y is the median of X.
%
%   P may be either a scalar or a vector. For scalar P, Y is a row
%   vector containing Pth percentile of each column of X. For vector P,
%   the ith row of Y is the P(i) percentile of each column of X.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 2.10 $  $Date: 2002/01/17 21:31:44 $

[prows pcols] = size(p);
if prows ~= 1 & pcols ~= 1
    error('P must be a scalar or a vector.');
end
if any(p > 100) | any(p < 0)
    error('P must take values between 0 and 100');
end

if (~any(isnan(x)))
   y = prctilecol(x,p);
else                    % if there are NaNs, process each column
   if (size(x,1) == 1)
      x = x';
   end
   c = size(x,2);
   np = length(p);
   y = zeros(np,c);
   for j=1:c
      xx = x(:,j);
      xx = xx(~isnan(xx));
      y(:,j) = prctilecol(xx,p)';
   end
   if (min(size(x)) == 1)
      y = y';
   end
end

function y = prctilecol(x,p);
xx = sort(x);
[m,n] = size(x);

if m==1 | n==1
    m = max(m,n);
        if m == 1,
           y = x*ones(length(p),1);
           return;
        end
    n = 1;
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx(:); max(x)];
else
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx; max(x)];
end

q = [0 q 100];
y = interp1(q,xx,p);
