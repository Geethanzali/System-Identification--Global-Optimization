%------------------------------------------------------%
% Simplex Method for                                   %
% non-negative objective function                      %
%                                                      %
% A. B. Smith                                          %
% absmith813@gmail.com                                 %
%------------------------------------------------------%
function [optx,x,f] = simplex(obj,bnd,tol)
search_time = cputime;
% Initialize -------------------------------------------
[n,d] = size(bnd);
maxitr = 2000;
stepx = 10;
f = [];
x = [];

spx = zeros(1,d);
[X,rmin] = obj.fx(spx); smin = spx;
while spx(1,1) < 1
    [X,r] = obj.fx(spx);
    if r < rmin, rmin = r; smin = spx; end
    
    k = d; while spx(1,k) > 1-stepx^-1 && k>1, k=k-1; end
    for i=k+1:d, spx(1,i) = 0; end
    spx(1,k) = spx(1,k)+stepx^-1;
end
spx = smin;

bndn = [(1-stepx^-1)*spx;(1+stepx^-1)*spx];
for i=1:d
    if bnd(1,i) < bndn(1,i), bnd(1,i) = bndn(1,i); end
    if bnd(2,i) > bndn(2,i), bnd(2,i) = bndn(2,i); end
end

spx = zeros(d+1,d);
fsx = zeros(d+1,1);
for i=1:d+1
    for j=1:d
        spx(i,j) = rand*(bnd(2,j)-bnd(1,j)) + bnd(1,j);
    end
    [X,fsx(i)] = obj.fx(spx(i,:));
    
    x = [x;spx(i,:)];
    f = [f;fsx(i)];
end

% Begin iteration for optimal solution -----------------
for i=1:maxitr 
    [fsx,P] = sort(fsx);
    spxc = spx;
    for j=1:d+1
        spx(j,:) = spxc(P(j),:);
    end
    
    % Center of mass -----------------------------------
    g = zeros(1,d);
    for j=1:d
        g = g + spx(j,:);
    end
    g = (1/d)*g;
    
    % Termination --------------------------------------
    sdx = g - spx(d+1,:);
    if sqrt(sdx*sdx') < tol
        break;
    end
    
    % Reflection ---------------------------------------
    r = g + (g - spx(d+1,:));
    frx = 0;
    c = spx(d+1,:) + 0.5*(g - spx(d+1,:));
    fcx = fsx(d+1)+1;
    for j=1:d % Bound check ---------------------------%
        if r(j) < bnd(1,j) || r(j) > bnd(2,j)          %
            frx = -1;                                  %
        end                                            %
    end       % ---------------------------------------%                           
    if frx ~= -1
        [X,frx] = obj.fx(r);
        if frx < fsx(1)
            % Expansion <><><><><><><><><><><><><><><><>
            e = g + 2*(g - spx(d+1,:));                %
            fex = 0;                                   %
            for j=1:d % Bound check -------------------%
                if e(j) < bnd(1,j) || e(j) > bnd(2,j)  %
                    fex = -1;                          %
                end                                    %
            end       % -------------------------------%
            if fex ~= -1                               %
                [X,fex] = obj.fx(e);                   %
                if fex < frx                           %
                    r = e;                             %
                    frx = fex;                         %
                end                                    %
            end                                        %
            % <><><><><><><><><><><><><><><><><><><><><>
        end    
    end
    
    if frx > fsx(d) || frx == -1
        % Contraction ----------------------------------
        frx = fsx(d+1)+1;
        [X,fcx] = obj.fx(c);
    end
    
    % Update set ---------------------------------------
    if frx > fsx(d) && fcx > fsx(d+1)
        % Shrink ---------------------------------------
        for j=2:d+1
            spx(j,:) = spx(1,:) + 0.5*(spx(j,:) - spx(1,:));
            [X,fsx(j)] = obj.fx(spx(j,:));
            
            x = [x;spx(j,:)];
            f = [f;fsx(j)];
        end
    else
        [C,I] = min([frx;fcx]);
        z = [r;c];
        spx(d+1,:) = z(I,:);
        fsx(d+1) = C;
        
        x = [x;spx(d+1,:)];
        f = [f;fsx(d+1)];
    end
    % End of Iteration ---------------------------------
    %if mod(i,50) == 0
    %fh = fopen('simplex.html','wt');
    %fprintf(fh,'Simplex method iteration %i complete.<br>\n',i);
    %fclose(fh);
    %end
end

optx = g;
search_time = cputime - search_time;
end