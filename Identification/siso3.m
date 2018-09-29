%-------------------------------------------------------------------%
% fast SISO 3nd-order class def                                     %
%                                                                   %
% iLS Ident                                                         %
% absmith813@gmail.com                                              %
%-------------------------------------------------------------------%
classdef siso3
    properties
        U = []
        Y = [];
        t = 0;
    end
    methods
    % Constructor ---------------------------------------------------
        function obj = siso3(U,Y)
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
            a(1) = ap(1) + ap(2) + ap(3); 
            a(2) = -(ap(1)*ap(2) + ap(1)*ap(3) + ap(2)*ap(3)); 
            a(3) = ap(1)*ap(2)*ap(3);
            
            % Information life parameter
            f = (a(1)^2 + a(2)^2 + a(3)^2 + 1)*ones(t-3,1);
            g = (a(1) - a(1)*a(2) - a(2)*a(3))*ones(t-3,1); 
            h = (a(2) - a(1)*a(3))*ones(t-3,1);
            
            f(2) = f(1) - g(1)^2/f(1);               g(2) = g(1) + g(1)*h(1)/f(1);                  h(2) = h(1) + g(1)*a(3)/f(1);  
            f(3) = f(1) - g(2)^2/f(2) - h(1)^2/f(1); g(3) = g(1) + g(2)*h(2)/f(2) + h(1)*a(3)/f(1); h(3) = h(1) + g(2)*a(3)/f(2); 
            for i=4:t-3
                f(i) = f(1) - g(i-1)^2/f(i-1) - h(i-2)^2/f(i-2) - a(3)^2/f(i-3);
                g(i) = g(1) + g(i-1)*h(i-1)/f(i-1) + h(i-2)*a(3)/f(i-2);
                h(i) = h(1) + g(i-1)*a(3)/f(i-1);  
            end
           
            % adjusted residual and input contributions
            r = -a(3)*Y(1:t-3) - a(2)*Y(2:t-2) - a(1)*Y(3:t-1) + Y(4:t);
            rrho = r;        rrho(2) = rrho(2) + rrho(1)*g(1)/f(1); rrho(3) = rrho(3) + rrho(2)*g(2)/f(2) + rrho(1)*h(1)/f(1);
            beta = U(1:t-3); beta(2) = beta(2) + beta(1)*g(1)/f(1); beta(3) = beta(3) + beta(2)*g(2)/f(2) + beta(1)*h(1)/f(1);
            gamm = U(2:t-2); gamm(2) = gamm(2) + gamm(1)*g(1)/f(1); gamm(3) = gamm(3) + gamm(2)*g(2)/f(2) + gamm(1)*h(1)/f(1);
            delt = U(3:t-1); delt(2) = delt(2) + delt(1)*g(1)/f(1); delt(3) = delt(3) + delt(2)*g(2)/f(2) + delt(1)*h(1)/f(1);
            for i=4:t-3
                rrho(i) = rrho(i) + rrho(i-1)*(g(i-1)/f(i-1)) + rrho(i-2)*(h(i-2)/f(i-2)) + rrho(i-3)*(a(3)/f(i-3));
                beta(i) = beta(i) + beta(i-1)*(g(i-1)/f(i-1)) + beta(i-2)*(h(i-2)/f(i-2)) + beta(i-3)*(a(3)/f(i-3));
                gamm(i) = gamm(i) + gamm(i-1)*(g(i-1)/f(i-1)) + gamm(i-2)*(h(i-2)/f(i-2)) + gamm(i-3)*(a(3)/f(i-3));
                delt(i) = delt(i) + delt(i-1)*(g(i-1)/f(i-1)) + delt(i-2)*(h(i-2)/f(i-2)) + delt(i-3)*(a(3)/f(i-3));
            end
           
            % gain solution matrix (2 x 2) b = y
            A = zeros(3,3);
            y = zeros(3,1);
            for i=1:t-3
                A(1,1) = A(1,1) + 0.5*(1/f(i))*delt(i)^2;
                A(2,2) = A(2,2) + 0.5*(1/f(i))*gamm(i)^2;
                A(3,3) = A(3,3) + 0.5*(1/f(i))*beta(i)^2;
                A(1,2) = A(1,2) + (1/f(i))*delt(i)*gamm(i);
                A(1,3) = A(1,3) + (1/f(i))*delt(i)*beta(i);
                A(2,3) = A(2,3) + (1/f(i))*gamm(i)*beta(i);
                y(1,1) = y(1,1) + rrho(i)*delt(i)/f(i);
                y(2,1) = y(2,1) + rrho(i)*gamm(i)/f(i);
                y(3,1) = y(3,1) + rrho(i)*beta(i)/f(i);
            end
            A = A+A';
            
            %zero reduction
            if(a(3) < tol), A = A(1:2,1:2); y = y(1:2,1); end
            if(a(2) < tol && a(3) < tol), A = A(1,1); y = y(1,1); end
            if   det(A) < tol, b = [0;0;0];
            else x = A\y; b = x; end
            while(length(b)<length(a)) b = [b;0]; end
            
            % Lagrange multipliers
            lm = 2*rrho - 2*[delt,gamm,beta]*b;
            lm(t-3) = -lm(t-3)/f(t-3);
            lm(t-4) = (g(t-3)*lm(t-3) - lm(t-4))/f(t-4);
            lm(t-5) = (h(t-3)*lm(t-3) + g(t-4)*lm(t-4) - lm(t-5))/f(t-5);
            for i=6:t-1
                lm(t-i) = (a(3)*lm(t-i+3) + h(t-i+2)*lm(t-i+2) + g(t-i+1)*lm(t-i+1) - lm(t-i))/f(t-i);
            end
            
            % generate X data from transfer function model
            XI = [0.5*(2*Y(3) - a(1)*lm(1) - a(2)*lm(2) - a(3)*lm(3));0.5*(2*Y(2) - a(2)*lm(1) - a(3)*lm(2));0.5*(2*Y(1) - a(3)*lm(1))];
            X = zeros(t,1);
            X(1) = XI(3); X(2) = XI(2); xn = XI(1);
            for i=3:t
                X(i) = xn;
                xn = a(1)*X(i) + a(2)*X(i-1) + a(3)*X(i-2) + b(1)*U(i) + b(2)*U(i-1) + b(3)*U(i-2);
            end
            
            % residual
            z = Y - X;
            z = z'*z;   
        end
    % H calc --------------------------------------------------------
        function [H] = hessian(obj,a)
            H = 1;
        end
    end
end