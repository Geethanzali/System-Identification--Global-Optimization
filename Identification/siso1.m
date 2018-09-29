%--------------------------------------------%
% fast SISO 1st-order solution               %
%                                            %
% iLS-ident                                  %
% absmith813@gmail.com                       %
%--------------------------------------------%
classdef siso1
    properties
        U = [];
        Y = [];
        
        t = 0;
    end
    methods
    % Constructor ---------------------------------------------------
        function obj = siso1(U,Y)
            obj.U = U;
            obj.Y = Y;

            % Parameter definition 
            t = length(Y);
            obj.t = t;    
        end
     % Function evaluation ------------------------------------------
        function [X,z,a,b,XI] = fx(obj,a)
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
            rrho = r;
            for i=2:t-1
                beta(i) = U(i) + (a/f(i-1))*beta(i-1);
                rrho(i) = r(i) + (a/f(i-1))*rrho(i-1);
            end
           
            A = 0; y = 0;
            for i=1:t-1
                A = A + f(i)^-1*beta(i)^2;
                y = y + f(i)^-1*beta(i)*rrho(i);   
            end
            if   A < tol, b = 0;
            else x = A\y; b = x; end
            
            % Lagrange multipliers
            lm = 2*rrho - 2*beta*b;
            lm(t-1) = -lm(t-1)/f(t-1);
            for i=2:t-1
                lm(t-i) = (a*lm(t-i+1) - lm(t-i))/f(t-i);
            end
            
            % generate X data from transfer function model
            XI = 0.5*(2*Y(1) - a*lm(1));
            X = zeros(t,1);
            
            xn = XI;
            for i=1:t
                X(i) = xn;
                xn = a*X(i) + b*U(i);
            end

            % residual
            z = (Y - X)'*(Y - X);
        end
    end
end
