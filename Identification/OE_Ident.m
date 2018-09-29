clear all;
clc;

cont2 = 1;
sys = 2;
 
a=0;
sim('id');

z = iddata(y5,u5);

     
                  % 1 Minute sampling time required
                   % Three columns Y, U and D
                  
MV = u5;   % Set point to aeration
DV = u5;   % Flow 
CV=  y5;        % Nitrate measurement

tst = 0:1:max(size(MV)-1);
tttst=tst;

GpMatrix = [ tf(0,[1 1]) tf(0,[1 1])];  % Initial data



% while cont2 == 1;
% figure(1)
% subplot(3,1,1)
% plot(CV(:,1))
% legend('CV')
% ylabel('Nitrate [units]')
% subplot(3,1,2)
% plot(MV(:,1))
% legend('MV')
% ylabel('Temperature [F]')
% subplot(3,1,3)
% plot(DV(:,1))
% xlabel('Time [s]')
% ylabel('Temperature [%]')


    
    estimind = input('estimate Gp (1) or estimate Gd (2) = ')
    
    CVst = CV;
    if estimind == 1; MVst=MV; end
    if estimind == 2; MVst=DV; end
      
     tst = 0:1:max(size(MV)-1);
     
     MVstsc = (MVst-mean(MVst))/(max(MVst-min(MVst)));
     CVstsc = (CV-mean(CV))/(max(CV-min(CV)));
     
%      figure(2)
% 
%      plot(tst,CVstsc,tst,MVstsc,'LineWidth',1);
%      xlabel('Time (min)','FontSize',20);
%      ylabel('Output','FontSize',20);
     
     %[xst,yst] = ginput ;  % Set crosshair
    
     tinit = 1;%round(xst(1),0);
     tfinal= 4000;%round(xst(2),0);
     
     if tfinal > max(length(MV)); tfinal = max(length(MV));end
     
     initial_time_final_time = [tinit tfinal]

     U=MVst(tinit:tfinal,1);  
     Y=CVst(tinit:tfinal,1); 

%   figure(2)
%   ttt=1:1:max(length(U));
%   yplot=(Y-mean(Y))/(max(Y)-min(Y));
%   uplot=(U-mean(U))/(max(U)-min(U));
%   plot(ttt,yplot,ttt,uplot); 
  %pause

smpt=1;
opt = [0,1,0,1,1,-3,-4,1,-1];

tol = 10^-7; 
t = length(Y);
MIND = 5; 

du = 1; 
%while du<t && abs(U(du+1)-U(du))<tol, du=du+1; end 
%if du > t - MIND || (length(U) ~= t)
    %fh = fopen('err.txt','at');
    %fprintf(fh,'<p><b>DATA ERROR:</b> Check to make sure dim(U) = dim(Y) and input not constant.</p>\n');
    %fclose(fh);
    %exit
%end

%if t < MIND
    %fh = fopen('err.txt','at');
    %fprintf(fh,'<p><b>DATA ERROR:</b> Check to make sure no. points > %i.</p>\n',MIND);
    %fclose(fh);
    %exit
%end

%Normalize data
umean = mean(U); ymean = mean(Y);
%U = U - umean; rngu = max(U) - min(U);
%Y = Y - ymean; rngy = max(Y) - min(Y);
%U = (1/rngu)*U; Y = (1/rngy)*Y;
FFMD = 0.50;
siso = siso1off(U,Y);
[X,z,a,b,XI,w] = siso.fx(FFMD);

%Dead-time estimation  (this does not work too great)

DSTEP = 10; DLENG = 40; STOP = 0; td = 0; 
tn = max(1,du-floor(0.5*DLENG)); %1
UN = U(tn:t); 
YN = Y(tn:t) - w*ones(t-tn+1,1); 
tt = length(YN); 

if tt < DSTEP+MIND, STOP = 1; end
if opt(9) ~=-1, td = opt(9); STOP = 1; end

DMIN = 0; DMAX = DMIN + DSTEP;
tmin = DMAX+1; tmax = min(tt,tmin+DLENG);
while ~STOP
    f = zeros(DMAX-DMIN+1,1);
    for i=DMIN:DMAX
        siso = siso1(UN(tmin-i:tmax-i),YN(tmin:tmax));
        [a,g,r] = golden(siso,[opt(1),opt(2)],tol);
        f(i-DMIN+1) = r(length(r));
    end
    [C,I] = min(f);
    
    if (I-1+DMIN) == DMAX, DMIN = DMAX; DMAX = DMIN + DSTEP; 
    else STOP = 1; td = I-1+DMIN; end
    
    tmin = DMAX+1; tmax = min(tt,tmin+DLENG);
    if((tmax-tmin)<MIND), STOP = 1; end
end

% td = Delayestimate;

dataid = iddata(Y,U,1);
GpMat = procest(dataid,'P1D'); % Matlab ID

UD = U(1:t-td); YD = Y(td+1:t);

%Order Estimation
OSST = 75; OSDST = 50; MORD = 2; RAVE = 2; CRITN = 1.4; %MORD= Maximum Order
odr = 1; siso = {}; 

siso{1}=siso1off(UD,YD);  
siso{2}=siso2(UD,YD);
siso{3}=siso3(UD,YD);                           %SISO1off with golden
opta1 = golden(siso{1},[opt(1),opt(2)],tol);    %opta1= optimum point of x with lowest residual
X = siso{odr}.fx(opta1); 
rf = (YD - X);

opta=opta1

% if opt(5) > 0
%    if opt(5) == 1, opta = opta1;
%    else
%        lbnd = zeros(1,opt(5)); ubnd = ones(1,opt(5));
%        lbnd(1,1) = opt(1); lbnd(1,2) = opt(3);
%        ubnd(1,1) = opt(2); ubnd(1,2) = opt(4);
%        opta = simplex(siso{opt(5)},[lbnd;ubnd],tol);
%    end
% end

% while opt(5) == -1 && odr<=MORD
%     lbnd = zeros(1,odr+1); ubnd = ones(1,odr+1);
%     lbnd(1,1) = opt(1); lbnd(1,2) = opt(3);
%     ubnd(1,1) = opt(2); ubnd(1,2) = opt(4);
%     
%     [opta2,g,r] = simplex(siso{odr+1},[lbnd;ubnd],tol);
%     for i=1:RAVE
%         [optav,g,rv] = simplex(siso{odr+1},[lbnd;ubnd],tol);
%         if rv(length(rv))<r(length(r)), r = rv; opta2 = optav; end
%     end
%     X = siso{odr+1}.fx(opta2);
%     rs = (YD - X);    
% 
%     k = OSST; ostst = 1;
%     while k < length(rf)
%         rfm = (rf(1:k)'*rf(1:k))/(YD'*YD);
%         rsm = (rs(1:k)'*rs(1:k))/(YD'*YD);
%         if rfm > CRITN*rsm, ostst = 0; end
%         k = k+OSDST;
%     end
%     
%     if(ostst == 1), opt(5) = odr; opta = opta1;
%     else odr = odr+1; rf = rs; opta1 = opta2; end
% end
if odr > MORD, opt(5) = MORD+1; opta = opta2; end

%type = input('Order: ');
switch opt(5)  %opt(5)
    case 1
        siso = siso1off(UD,YD);
        [X,g,a,b,xi,w] = siso.fx(opta);
        
        %H = siso.hessian(opta);
    case 2
        siso = siso2(UD,YD); w = 0;
        [X,g,a,b,xi] = siso.fx(opta); 
        %H = siso.hessian(opta);
    case 3
        siso = siso3(UD,YD); w = 0;
        [X,g,a,b,xi] = siso.fx(opta); 
        %H = siso.hessian(opta);
end

Gd = tf(b,[1 -a],1); %(num, denom, sampletime in discrete)
[num,den] = tfdata(d2c(Gd)); %d2c- produces cont time model from discrete time
Gp=tf(num,den,'inputdelay',td);

[nDV,dDV]=tfdata(Gp,'v'); %returns numerator & denominator

GpMatrix(estimind)=Gp

GpMatlab(estimind) = GpMat;

z1 = dtrend(dataid);
opt = procestOptions;
opt.Focus = 'prediction';
m= procest(z1,'p1d','Tolerance',0.000000001,opt);

%using ARMA model
opt = procestOptions;
opt.Focus = 'prediction';
opt.DisturbanceModel='ARMA1';
opt.Regularization.Lambda = 1; 
m1 = procest(dataid,'p1d','Tolerance',0.000000001,opt);

compare(dataid,m,m1,GpMatrix(estimind),GpMatlab(estimind))
legend('Actual','Detrend','ARMA','GpMatrix','GpMatlab')


