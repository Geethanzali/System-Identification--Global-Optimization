%--------------------------------------------
% fast SISO 1st-order solution
% with time-delay
%
% iLS-ident
% A. B. Smith
% absmith813@gmail.com
%
% Description: Driver for SISO ident. First a delay is estimated based on a
% heuristic. Then the optimal subproblem is solved iteratively to determine
% an appropriate time constant. 
%--------------------------------------------
function [a,b,xi,d,r,X,v,da,dda] = fast(U,Y,WDELAY)
    t = length(Y);
    % time-constant bounds
    AUBND = 1;
    ALBND = 0;
    % dead-time estimation parameters
    DSTEP = 10;
    % opt tolerance
    tol = 10^-7;
    
    d = 0;
    STOP = ~WDELAY;
    DMIN = 0;
    DMAX = DMIN + DSTEP;
    % Delay estimation, the idea here is as follows: start with a discrete
    % delay DMIN, solve the optimization problem for a, b, xi, and v, then
    % increment DMIN and repeat until DMAX = DMIN+DSTEP is reached. Go back to
    % the solutions computed and locate the lowest value of the objective
    % function/residual, (Y-X-v)^T(Y-X-v). If this minimum is on the interior
    % of the interval DMIN to DMAX, select this minimum as the delay.
    % Otherwise, if the minimum lies on the edge of the interval, set DMIN
    % = DMAX, DMAX = DMIN + DSTEP and repeat the process. 
    while ~STOP
        f = zeros(DMAX-DMIN+1,1);
        a = zeros(DMAX-DMIN+1,1);
        for i=DMIN:DMAX
            siso = fsiso(U(DMAX-i+1:t-i),Y(DMAX+1:t));
            [a(i-DMIN+1),g,r] = golden(siso,[ALBND,AUBND],tol);
            f(i-DMIN+1) = r(length(r));
        end
        [C,I] = min(f);
        
        if (I-1+DMIN) == DMAX
            DMIN = DMAX;
            DMAX = DMIN + DSTEP;
        else
            STOP = 1;
            d = I-1+DMIN;
        end
    end
    % NOTE: this routine could now be replaced with the
    % quality estimate, ie Q = 1 - integral{f_u=U(a)/f_u=C(a) da}
    
    % Setup the ident object
    siso = fsiso(U(1:t-d),Y(d+1:t)); 
    % Solve for the time constant using a golden ratio line search
    [a,g,r] = golden(siso,[ALBND,AUBND],tol); 
    % Calculate the subproblem solution at the optimal a
    [X,r,b,xi,v] = siso.fx(a);
    
    % basic numerical derivative calculation
    xstep = 10^-3;
    [X,rf] = siso.fx(a+xstep);
    [X,rb] = siso.fx(a-xstep);
    
    da = (rf - rb)/(2*xstep);
    dda = (rf + rb - 2*r)/(xstep^2);
end
