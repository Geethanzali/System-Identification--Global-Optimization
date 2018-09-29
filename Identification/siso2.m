%-------------------------------------------------------------------%
% fast SISO 2nd-order class def                                     %
%                                                                   %
% iLS Ident                                                         %
% absmith813@gmail.com                                              %
%-------------------------------------------------------------------%
classdef siso2
    properties
        U = []
        Y = [];
        t = 0;
    end
    methods
    % Constructor ---------------------------------------------------
        function obj = siso2(U,Y)
            obj.U = U;
            obj.Y = Y;

            % Parameter definition
            t = length(Y);
            obj.t = t;  
        end
    % Function evaluation -------------------------------------------
        function [X,z,a,b,XI] = fx(obj,a)
            U = obj.U; Y = obj.Y;
            t = obj.t;
            tol = 10^-7;
            
            ap = a;
            a(1) = ap(1) + ap(2); 
            a(2) = -ap(1)*ap(2);
            
            % Information life parameter
            f = (a(1)^2+a(2)^2+1)*ones(t-2,1); 
            g = a(1)*(1-a(2))*ones(t-2,1);
            
            f(2) = f(1) - g(1)^2/f(1); g(2) = g(1) + a(2)*g(1)/f(1); 
            for i=3:t-2
                f(i) = f(1) - g(i-1)^2/f(i-1) - a(2)^2/f(i-2);
                g(i) = g(1) + a(2)*g(i-1)/f(i-1); 
            end

            % adjusted residual and input contributions
            r = -a(2)*Y(1:t-2) - a(1)*Y(2:t-1) + Y(3:t);
            rrho = r; rrho(2) = rrho(2) + rrho(1)*g(1)/f(1);
            beta = U(1:t-2); beta(2) = beta(2) + beta(1)*g(1)/f(1);
            gamm = U(2:t-1); gamm(2) = gamm(2) + gamm(1)*g(1)/f(1);
            for i=3:t-2
                rrho(i) = rrho(i) + rrho(i-1)*(g(i-1)/f(i-1)) + rrho(i-2)*(a(2)/f(i-2));
                beta(i) = beta(i) + beta(i-1)*(g(i-1)/f(i-1)) + beta(i-2)*(a(2)/f(i-2));
                gamm(i) = gamm(i) + gamm(i-1)*(g(i-1)/f(i-1)) + gamm(i-2)*(a(2)/f(i-2));
            end
           
            % gain solution matrix (2 x 2) b = y
            A = zeros(2,2);
            y = zeros(2,1);
            for i=1:t-2
                A(1,1) = A(1,1) + (1/f(i))*gamm(i)^2;
                A(1,2) = A(1,2) + (1/f(i))*beta(i)*gamm(i);
                A(2,2) = A(2,2) + (1/f(i))*beta(i)^2;
                y(1,1) = y(1,1) + rrho(i)*gamm(i)/f(i);
                y(2,1) = y(2,1) + rrho(i)*beta(i)/f(i);
            end
            A(2,1) = A(1,2); 
            if   det(A) < tol, b = [0;0];
            else x = A\y; b = x; end
            
            % Lagrange multipliers
            lm = 2*rrho - 2*[gamm,beta]*b;
            lm(t-2) = -lm(t-2)/f(t-2);
            lm(t-3) = (g(t-3)*lm(t-2) - lm(t-3))/f(t-3);
            for i=4:t-1
                lm(t-i) = (g(t-i)*lm(t-i+1) + a(2)*lm(t-i+2) - lm(t-i))/f(t-i);
            end
            
            % generate X data from transfer function model
            XI = [0.5*(2*Y(2) - a(1)*lm(1) - a(2)*lm(2));0.5*(2*Y(1) - a(2)*lm(1))];
            X = zeros(t,1);
            X(1) = XI(2); xn = XI(1);
            for i=2:t
                X(i) = xn;
                xn = a(1)*X(i) + a(2)*X(i-1) + b(1)*U(i) + b(2)*U(i-1);
            end
            
            
            % residual
            z = Y - X;
            z = z'*z;   
        end
    % H calc --------------------------------------------------------
        function [H] = hessian(obj,a)
            tol = 10^-4;
            
            [X,g]    = obj.fx(a);
            [X,bcs1] = obj.fx([a(1)-tol,a(2)]);
            [X,fws1] = obj.fx([a(1)+tol,a(2)]);
            [X,bcs2] = obj.fx([a(1),a(2)-tol]);
            [X,fws2] = obj.fx([a(1),a(2)+tol]);
            fds1 = (fws1 - g)/(tol); fds2 = (fws2 - g)/(tol);
            bds1 = (g - bcs1)/(tol); bds2 = (g - bcs2)/(tol);
            H = [(fds1-bds1)/(tol),(fds1-bds2)/(tol),(fds2-bds1)/(tol),(fds2-bds2)/(tol)];
            
            maxd = max(H(1,2),H(1,3));
            H(1,2) = maxd; H(1,3) = maxd;
        end
    end
end