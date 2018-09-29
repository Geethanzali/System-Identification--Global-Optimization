%--------------------------------------------
% golden line search
%
% iLS-ident
% A. B. Smith
% absmith813@gmail.com
%
% Description: This is a standard golden line search script. The basic 
% idea is as follows: We wish to find the minimum of a monotonically
% decreasing function, so we partition into intervals which obey the
% golden ratio, and systematically converge the interval onto the solution
% by throwing away the largest-valued interval. This routine requires that
% the object passed (obj) have an (.fx) method, which produces a function 
% evaluation. 
%--------------------------------------------
function [optx,x,f] = golden(obj,bnd,tol)
% Initialize ------------------------------------
g = (1+sqrt(5))/2;
a = bnd(1);
b = bnd(2);
c = a + (b-a)/(1+g);
d = c + (b-c)/(1+g);

[X,fc] = obj.fx(c);
x = [b-a];
f = [fc];
% Begin iteration for optimal solution ----------
maxitr = 500;
fd=0;

for i=1:maxitr
    % Termination
    if (b-a) < tol
        break;
    end
    % Test point
    if (b-c) > (c-a)
        d = c + (b-c)/(1+g);
    else
        d = c + (a-c)/(1+g);
    end
    [X,fd] = obj.fx(d);
    
    % Bnd update
    h = d;
    if fc > fd
        h = c;
        c = d;
        fc = fd;
    end
    if h<c
        a = h;
    else
        b = h;
    end
    x = [x;b-a];
    f = [f;fd];   
end

optx = (b+a)/2;
end