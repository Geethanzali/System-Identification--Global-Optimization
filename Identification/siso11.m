classdef siso11
    properties
        U = [];
        Y = [];

        t = 0;
    end
    methods
    % Constructor ---------------------------------------------------
        function obj = siso11(U,Y)
            obj.U = U;
            obj.Y = Y;
            %obj.StartStep = StartStep;

            % Parameter definition
            t = length(Y);
            obj.t = t;  
        end
    % Function evaluation -------------------------------------------
        function [X,z,a,b,XI] = fx(obj,a,StartStep)

            UD = obj.U; YD = obj.Y;
            t = obj.t;
            %StartStep = obj.StartStep;
            tol = 10^-7;
            %StartStep = StartStep;
            %ap = a;
            %a = zeros(1,3);
            %StartStep = a(4);
%             a(1) = ap(1) + ap(2) + ap(3); 
%             a(2) = -(ap(1)*ap(2) + ap(1)*ap(3) + ap(2)*ap(3)); 
%             a(3) = ap(1)*ap(2)*ap(3);
            
%             dumdum=zeros(2,1);
%             dumdum=[eye(2) dumdum];
%             A=[a;dumdum];

            A=a;

            B=1;
            X0=zeros(1,1);
            dumt=length(YD);
            qmat=zeros(dumt,1);

            dumA=eye(1);
            %Av=dumA;
            dumB=B;
            %Bv=dumB;
            qmat(StartStep,:)=(dumA*X0).';

            %dumt=7;

            for i=StartStep+1:dumt
                dumA=dumA*A;
                if i==StartStep+1
                    dumB=B;
                else
                    dumB=[A*dumB B];
                end

                %Av=[Av;dumA];
                %Bv=[Bv zeros((i-1)*2,1);dumB];
                X0=[X0;UD(i-1)];
                qmat(i,:)=([dumA dumB]*X0).';
            end

            dumW=StartStep-1;
            Weight=[1*ones(dumW,1);1*ones(dumt-dumW,1)];
            WD=diag(Weight);
            
            %EE=zeros(length(C));EE(1,1)=-1;EE(2,2)=-1;EE(3,3)=1;
            %CC=iLS_qp(qmat.'*WD*qmat,-qmat.'*WD*YD,[],[],EE,zeros(size(C)),C);
            if max(max(abs(qmat)))==0
                b=zeros(1,1);
            else
                b=inv(qmat.'*WD*qmat)*qmat.'*WD*YD;
            end

            X = qmat*b;
            z = YD - X;
            z = z.'*z;
            
            XI = zeros(1,1);
            
            
        end
    end
end