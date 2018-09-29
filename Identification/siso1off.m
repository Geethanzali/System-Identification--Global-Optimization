%--------------------------------------------%
% fast SISO 1st-order solution               %
%                                            %
%                                            %
% iLS-ident                                  %
% absmith813@gmail.com                       %
%--------------------------------------------%
classdef siso1off
    properties
        U = [];
        Y = [];
        
        t = 0;
    end
    methods
    % Constructor ---------------------------------------------------
        function obj = siso1off(U,Y)
            obj.U = U;
            obj.Y = Y;

            % Parameter definition 
            t = length(Y);
            obj.t = t;    
        end
    % Function evaluation -------------------------------------------
        function [X,z,a,b,XI,v] = fx(obj,a)
            U = obj.U; Y = obj.Y;
            t = obj.t;
            tol = 10^-7;
            
            f = ones(t-1,1);
            f(1) = f(1) + a^2; itr = 2;
            while (f(itr-1) - 1) > tol && itr < t
                f(itr) = f(1) - a^2/f(itr-1); itr = itr+1;
            end
            
            r = -a*Y(1:t-1) + Y(2:t);
            beta = U(1:t-1);
            phii = (-a+1)*ones(t-1,1);
            rrho = r;
            for i=2:t-1
                beta(i) = U(i) + (a/f(i-1))*beta(i-1);
                phii(i) = phii(i) + (a/f(i-1))*phii(i-1);
                rrho(i) = r(i) + (a/f(i-1))*rrho(i-1);
            end
           
            A = zeros(2,2); y = zeros(2,1);
            for i=1:t-1
                A(1,1) = A(1,1) + f(i)^-1*beta(i)^2;
                A(1,2) = A(1,2) + f(i)^-1*beta(i)*phii(i);
                A(2,2) = A(2,2) + f(i)^-1*phii(i)^2;
                y(1,1) = y(1,1) + f(i)^-1*beta(i)*rrho(i);
                y(2,1) = y(2,1) + f(i)^-1*phii(i)*rrho(i);
            end
            A(2,1) = A(1,2);
            if   det(A) < tol, b = 0; v = y(2,1)/A(2,2);
            else x = A\y; b = x(1); v = x(2); end
            
            % Lagrange multipliers
            lm = 2*rrho - 2*beta*b - 2*phii*v;
            lm(t-1) = -lm(t-1)/f(t-1);
            for i=2:t-1
                lm(t-i) = (a*lm(t-i+1) - lm(t-i))/f(t-i);
            end
            
            % generate X data from transfer function model
            XI = 0.5*(2*Y(1) - a*lm(1) - 2*v);
            X = zeros(t,1);
            
            xn = XI;
            for i=1:t
                X(i) = xn;
                xn = a*X(i) + b*U(i);
            end
            X = X + v*ones(t,1);
            
            % residual
            z = (Y - X)'*(Y - X);
        end
    % H calc --------------------------------------------------------
        function [H] = hessian(obj,a)
            tol = 10^-4;
            
            [X,g]   = obj.fx(a);
            [X,bcf] = obj.fx(a-tol);
            [X,fwf] = obj.fx(a+tol);
            fdf = (fwf - g)/(tol);
            bdf = (g - bcf)/(tol);
            H = (fdf - bdf)/(tol);
        end
    end
end
