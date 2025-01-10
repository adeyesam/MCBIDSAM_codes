%% The file loads and simulates hierarchically ranked models (based on AICc+2Cve values)
% from the branch and bound algorithm. 'solnrank' is the rank of model we desire to simulate.
% The rows of 'inputs' are the input variables while the columns are the time instants


function y = simulate_model_dynamic2(inputs, y_initial, solnrank, constr)

load('solution_workspaceC')
u = inputs;
[lengthu, lengthu2] = size(u);
N2 = N^2;
xqs(:,1) = y_initial;

F = 4;
nci = [1 0 1 0];                   % No. of C in each input variable
ncio = [1 0 1 0];                  % No. of C in each output variable
nhi = [0 2 0 2];                   % Net No. of H in each species i on LHS 
nhio = [0 2 0 2];                  % Net No. of H in each species i on RHS

% Scaling inputs and initial outputs
for i = 1:lengthu
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end
for i=1:N
    y_initial(i,:)=1 + 2*(y_initial(i,:)-mnd(i))/(mxd(i)-mnd(i));
end

XXQ = y_initial;
XQ(:,1) = XXQ;

nddr = 0;
for jj=1:lengthu2
    
    % Compute input transformations
    if lengthu>1
        i_UiUj = nchoosek(1:lengthu,2);                     % Indices for Mixed second order inputs
        UiUj = zeros(size(i_UiUj,1),1);
        for i = 1:size(i_UiUj,1)
                UiUj(i,:) = u(i_UiUj(i,1),jj) .* u(i_UiUj(i,2),jj);
        end
    end
    n_uy = max(N,lengthu);
    
    mii=[];
    for i=1:min(N,lengthu)
        mii=[mii;[i i]];
    end
    iUY = nchoosek(1:n_uy,2);                            % Combination of indices if n(u)=n(measured outputs)=n_uy
    iUY=[iUY;flip(iUY,2)];
    [~,iUY1]=sort(iUY(:,1));
    iUY=[mii;iUY(iUY1,:)];
    if n_uy==N                                          % Eliminate indices aloted to u greater than n(u)
        elim = find(iUY(:,1)>lengthu);
        iUY(elim,:) = [];
    else                                                % Eliminate indices alotted to y greater than n(measured outputs)
        elim = find(iUY(:,2)>N);
        iUY(elim,:) = [];
    end
    UY = zeros(size(iUY,1),1);
    for i=1:size(iUY,1)
        UY(i,:) = u(iUY(i,1),jj).*XXQ(iUY(i,2),:);     % Compute U*Y
    end
    U2 = u(:,jj) .* u(:,jj);
    LogU = log(u(:,jj));                                % Natural log transformation
    ExpU = exp(-u(:,jj));                               % Exponential transformation
    InvU = 1./u(:,jj);                                  % Inverse of inputs
    InvUsq = (u(:,jj).*u(:,jj)).\1;                     % Inverse of squared inputs
    sqrtU = sqrt(u(:,jj));                              % Square root of inputs
    sqrtInvU = sqrtU.\1;                                % 1/sqrt(u)
    sigmoidU = (1+exp(-u(:,jj))).\1;                    % Sigmoid

    if lengthu>1
        Utrans = [ones(1,1); u(:,jj); UY; UiUj; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
    else
        Utrans = [ones(1,1); u(:,jj); UY; U2; LogU; ExpU; InvU; InvUsq; sqrtU; sqrtInvU; sigmoidU];
    end
    Ntrans = size(Utrans,1);                            % NUMBER OF TRANSFORMATIONs
    
    % Arranging Utrans according to rank
    Utrans_r = Utrans(idOb,:);                          

    if Ntrans>85
        Ntrans = 85;
        Utrans_r = Utrans_r(1:Ntrans,:);
    end

    % Model parameters for particular rank
    sln = solnrank;
    indx = find(Obtrnd(sln,:));

    Cm1=Cmat(sln,:);                 Am1=Amat(sln,:);  
    Cm_1=reshape(Cm1,Ntrans,N)';      A1=reshape(Am1,N,N)';
    C1=Cm_1(:,indx);      
    
    Cm2=Cmatc(sln,:);                 Am2=Amatc(sln,:);  
    Cm_2=reshape(Cm2,Ntrans,N)';      A2=reshape(Am2,N,N)';
    C2=Cm_2(:,indx);      
    
    w = wmat(:,1,sln);
    
    % Model simulation
    Utrans_ = Utrans_r(indx,:);                         % Selecting relevant input transformations
    if constr
        if jj>2                                         % define steady state
            flag1 = all(abs(1 - u(3:4,jj)./u(3:4,jj-1))<0.002);
            flag2 = all(abs(1 - xqs(3:4,jj)./xqs(3:4,jj-1))<0.02);
        end
        
        if jj==1 || jj==2                               % initial
            XQ(:,jj+1) = A1*XQ(:,jj) + C1*Utrans_;
        elseif flag1 && flag2                           % will use constrained model for steady state
            XQ(:,jj+1) = A1*XQ(:,jj) + C1*Utrans_;
        else
            XQ(:,jj+1) = A1*XQ(:,jj) + C1*Utrans_;
            nddr = 0;
        end

        for j=1:size(XQ,1)
            xqs(j,jj+1) = mnd(j) + 0.5*(XQ(j,jj+1)-1)*(mxd(j)-mnd(j));
        end
        
        if jj>2 && flag1 && flag2 && nddr==0
            nddr = nddr+1;
            objfun = @(x)vddr_obj(x, XQ(:,jj+1), N, w);
            nlcon = @(x)vddr_nlcon(x, inputs(:,jj), nci, ncio, N, F, mnd, mxd);
        
            nlrhs = 0; 
            nle = 0;  
            
            lb = 0.8*XQ(:,jj+1); 
            ub = 1.2*XQ(:,jj+1); 
        
            yddr0 = XQ(:,jj+1);
    
            optns = optiset('solver', 'ipopt', 'display', 'iter', 'maxiter', 1e2, 'maxtime', 1000);
            Optim = opti('fun', objfun, 'ineq', [], [], 'nlmix', nlcon, nlrhs, nle, 'bounds', lb, ub, 'options', optns)

            for i=1:N 
                xxs(i,jj+1) = mnd(i) + 0.5*(yddr_(i,jj+1)-1)*(mxd(i)-mnd(i));
            end

            cin = nci*inputs(:,jj);
            cout = ncio*xxs(:,jj+1);

            XQ(:,jj+1) = yddr_(:,jj+1);
    
            for j=1:size(XQ,1)
                xqs(j,jj+1) =  mnd(j) + 0.5*(XQ(j,jj+1)-1)*(mxd(j)-mnd(j));
            end
            
        elseif jj>2 && flag1 && nddr>0
            XQ(:,jj+1) = yddr_(:,jj+1);
            for j=1:size(XQ,1)
                xqs(j,jj+1) =  mnd(j) + 0.5*(XQ(j,jj+1)-1)*(mxd(j)-mnd(j));
            end
        end
    else
        if jj==1
            XQ(:,jj+1) = A1*XXQ + C1*Utrans_;
        else
            XQ(:,jj+1) = A1*XQ(:,jj) + C1*Utrans_;
        end
        for j=1:size(XQ,1)
            xqs(j,jj+1) = mnd(j) + 0.5*(XQ(j,jj+1)-1)*(mxd(j)-mnd(j));
        end
    end
    XXQ = XQ(:,jj+1);

end
if constr
    figure, plot(nlcon_)
end
y = xqs;
% Model interpretation
fprintf(['\nThe identified model is of the form \n \n y(k+1) = A*y(k) + C*Utrans(k) \n \n'])
fprintf(['Selected basis functions in Utrans for model ', num2str(sln), ': \n \n']);
BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,sln,true);

fprintf(['\n Model parameters: \n']);
if constr
    A2
    C2
else
    A1
    C1
end


return

