%--------------------------------------------
% fast SISO 1st-order solution
% with time-delay
%
% iLS-ident
% A. B. Smith
% absmith813@gmail.com
%
% Description: Matlab class (haha). This object solves the first order
% SISO ident problem optimally for fixed time constant. 
% See comments in code.
%--------------------------------------------
classdef fsiso
    properties
        U = [];     % property for input data
        Y = [];     % property for output data
        
        t = 0;      % number of discrete time steps under consideration
    end
    methods
% Constructor ------------------------------------------------------------
        function obj = fsiso(U,Y) 
        %this method reads in U,Y data
            obj.U = U;
            obj.Y = Y;
% -----------------------------------------
% Parameter definition --------------------
            t = length(Y);
            obj.t = t;    
        end
% Function evaluation ----------------------------------------------------
        function [X,z,b,XI,v] = fx(obj,a) 
        % This method solves the optimization subproblem by directly 
        % computing the solution to the linear system: 
        % df/dxi + sum_ lm*dg/dxi = 0 
        % g = 0,
        % or in street lingo; the Lagrangian ... well maybe just the
        % streets of France where they were passing out Lecons sur le
        % calcul des fonctions.
        % X is the optimal process trajectory, z the residual, b the gain,
        % XI the initial condition for X, and v the offset.
            U = obj.U;
            Y = obj.Y;
            t = obj.t;
            tol = 10^-7;
            
        % From here till where I denote END, the computations are focused,
        % literally, on a detailed, symbolic solution of the nonzero 
        % elements of the linear systems arising from the above 
        % optimization subproblem.
            g = (1+a^2)*ones(t-1,1);
            itr = 2;
            while (a^2/g(itr-1)) > tol && itr < t
                g(itr) = g(1) - a^2/g(itr-1);
                itr=itr+1;
            end
            
            r = -a*Y(1:t-1) + Y(2:t);
            beta = U(1:t-1);
            phii = (-a+1)*ones(t-1,1);
            rrho = r;
            for i=2:t-1
                beta(i) = U(i) + (a/g(i-1))*beta(i-1);
                phii(i) = phii(i) + (a/g(i-1))*phii(i-1);
                rrho(i) = r(i) + (a/g(i-1))*rrho(i-1);
            end
           
            A = zeros(2,2); y = zeros(2,1);
            for i=1:t-1
                A(1,1) = A(1,1) + g(i)^-1*beta(i)^2;
                A(1,2) = A(1,2) + g(i)^-1*beta(i)*phii(i);
                A(2,2) = A(2,2) + g(i)^-1*phii(i)^2;
                y(1,1) = y(1,1) + g(i)^-1*beta(i)*rrho(i);
                y(2,1) = y(2,1) + g(i)^-1*phii(i)*rrho(i);
            end
            A(2,1) = A(1,2);
            
            % Check the conditioning of A, if poorly scaled, just set b
            % equal to zero. This was to avoid various scaling issues
            % matlab was reporting, however, it is simply a quick fix. 
            if   det(A) < tol, b = 0; v = y(2,1)/A(2,2);
            else x = A\y; b = x(1); v = x(2); end
            % END
            
            % Lagrange multipliers, lm, in the above subproblem statement
            lm = 2*rrho - 2*beta*b - 2*phii*v;
            lm(t-1) = -lm(t-1)/g(t-1);
            for i=2:t-1
                lm(t-i) = (a*lm(t-i+1) - lm(t-i))/g(t-i);
            end
            
            % generate X data from transfer function model
            XI = 0.5*(2*Y(1) - a*lm(1) - 2*v);
            X = zeros(t,1);
            xn = XI;
            for i=1:t
                X(i) = xn;
                xn = a*X(i) + b*U(i);
            end

            % residual
            z = (Y - X - v)'*(Y - X - v);
        end
    end
end