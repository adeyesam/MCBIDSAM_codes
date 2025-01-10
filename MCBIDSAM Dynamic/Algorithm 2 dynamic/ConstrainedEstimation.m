function [Obj_fun, A_mat, C1_mat, C2_mat, D_mat, l_mat, fval_trend] = ConstrainedEstimation(indx, iUY,i_UiUj,UYindex,Ncons, Ain, C1in, C2in, Din, Lin, nci, ncio, nhi, nhio)

load('data4est.mat','data_est','u_est','utrans_est','idOb','nut','mnd','mxd', 'umin', 'umax')
data = data_est;
u = u_est;
Utrans = utrans_est;
clear data_est u_est utrans_est

f_indx = find(indx);

%% Estimation and model predictions
sgma = 0.06;   

tic
    %% define the parameters
    lengthu = size(u,1);
    Nwmcon = size(data,1);
    N = Nwmcon - Ncons;                                     % NUMBER OF OUTPUT
    N2 = N^2;                                               % NUMBER OF ELEMENTS IN A
    Ntrans = size(Utrans,1);                                % NUMBER OF TRANSFORMATIONs
    ntrans = numel(f_indx);
    lengthu = size(u,1);

    AtomIN = data(N+1:end,:);
    data = data(1:N,:);

    %% Indices for parameters in matrices
    setA= 1:N2;
    setC1 = 1:N*lengthu;
    setC2 = [];
    for stc = 1:N
        setC2 = [setC2 f_indx+(stc-1)*Ntrans];
    end
    
    numsetA = numel(setA);         
    ar_a = zeros(numsetA,1);
    ac_a = ar_a;

    numsetC1 = numel(setC1);         
    numsetC2 = numel(setC2);    
    
    ar_c1 = zeros(numsetC1,1);
    ac_c1 = ar_c1;
    ar_c2 = zeros(numsetC2,1);
    ac_c2 = ar_c2;
    
    setD = 1:N;
    numsetD = numel(setD);
    ar_d = (1:N)';
    ac_d = ones(N,1);
    
    lenreading = size(data,2)-1;  
    lenreadpone = lenreading+1;     
    Lenvec = numsetA + numsetC1 + numsetC2 + numsetD;               
    sz = lenreadpone*N;                        
    szvec = sz+Lenvec;                         
    vecparams = zeros(Lenvec,1); 
    szlam = N^2 + (N-N^2)/2;                      % Number of covariance
    meanprior = vecparams;
    
    %% define C and D and obtain vecparams for reduced system
    A = reshape(Ain,[N,N])';
    C1 = reshape(C1in,[lengthu,N])';
    C2 = reshape(C2in,[Ntrans,N])';
    D = Din';
    lam = Lin';

    for i = 1:numsetA
        ar_a(i) = ceil(setA(i)/N);
        ac_a(i) = N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N);
        vecparams(i) = A(ar_a(i),ac_a(i));
    end

    for i = 1:numsetC1
        ar_c1(i) = ceil((setC1(i))/lengthu);
        ac_c1(i) = setC1(i) - lengthu*fix((setC1(i))/lengthu); 
        if ac_c1(i) == 0
            ac_c1(i) = lengthu;
        end
        vecparams(i+numsetA) = C1(ar_c1(i),ac_c1(i));
    end
    
    for i = 1:numsetC2
        ar_c2(i) = ceil((setC2(i))/Ntrans);
        ac_c2(i) = setC2(i) - Ntrans*fix((setC2(i))/Ntrans); 
        if ac_c2(i) == 0
            ac_c2(i) = Ntrans;
        end
        vecparams(i+numsetA+numsetC1) = C2(ar_c2(i),ac_c2(i));
    end

    for i = 1:numsetD
        vecparams(i+numsetA+numsetC1+numsetC2) = D(ar_d(i),ac_d(i));
    end
    
    %% Initialize xq
    ybml = data(1:N,:);
    for i = 2:lenreadpone
        ybml(:,i) = A*ybml(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
    end
    ybmls = ybml;
    datas = data;
    for i=1:N 
        ybmls(i,:) = mnd(i) + 0.5*(ybml(i,:)-1)*(mxd(i)-mnd(i));
        datas(i,:) = mnd(i) + 0.5*(data(i,:)-1)*(mxd(i)-mnd(i));
    end
    
    %% Jacobian remains unchanged with model parameters since model is linear in parameters
    dAA = zeros(sz,N);
    dCC1 = zeros(sz,lengthu);   dCC2 = zeros(sz,ntrans);  dDD = zeros(sz,N);
    for ij = 1:N
        dAA((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*N+1:ij*N) = data(:,1:end-1)';
        dCC1((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*lengthu+1:ij*lengthu) = u(:,1:end-1)';
        dCC2((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*ntrans+1:ij*ntrans) = Utrans(f_indx,1:end-1)';
        dDD((ij-1)*lenreading + 1:ij*lenreading,ij) = ones(lenreading,1);
    end

    hwhole = [dAA dCC1 dCC2 dDD];
    clear dAA dCC1 dCC2 dDD

    %% Constranied parameter estimation    
    objfun = @(x)obj3(x, data, Utrans, u, N, sz, lenreading, numsetA, numsetC1, numsetC2, numsetD, ar_a, ac_a, ar_c1, ac_c1, ar_c2, ac_c2, ar_d, ac_d);
    nlcon = @(x)mass_bal_constr3(x, indx, N, lengthu, Ncons,  nci, ncio, nhi, nhio, mnd, mxd, umin, umax);

    nlrhs = zeros((N+lengthu+ntrans)*Ncons+Ncons,1); 
    nle = nlrhs; 

    lb = vecparams - 0.1e3*abs(vecparams); 
    ub = vecparams + 0.1e3*abs(vecparams); 

    if numsetD==0
        x0 = [reshape(A',N2,1); reshape(C1',numsetC1,1); reshape(C2(:,f_indx)',Lenvec,1)]; %oldparams
    else
        x0 = [reshape(A',N2,1); reshape(C1',numsetC1,1); reshape(C2(:,f_indx)',numsetC2,1); D]; %oldparams
    end

    optns = optiset('solver', 'ipopt', 'display', 'iter', 'maxiter', 1.22e2, 'maxtime', 2800);
    Optim = opti('fun', objfun, 'ineq', [], [], 'nlmix', nlcon, nlrhs, nle, 'bounds', lb, ub, 'options', optns)
    
    [x_new,fval,exitflag,info] = solve(Optim,x0);

    vecparams = x_new(1:Lenvec);
        
    fval_trend = fval;
    %% Update A,C,D
    for j = 1:numsetA
        ar = ar_a(j);
        ac = ac_a(j);
        A(ar,ac) = vecparams(j);
    end
    for j = 1:numsetC1
        ar = ar_c1(j);
        ac = ac_c1(j);
        C1(ar,ac) = vecparams(j+numsetA);
    end

    for j = 1:numsetC2
        ar = ar_c2(j);
        ac = ac_c2(j);
        C2(ar,ac) = vecparams(j+numsetA+numsetC1);
    end

    for k = 1:numsetD
        dr = ar_d(k);
        dc = ac_d(k);
        D(dr,dc) = vecparams(k+numsetA+numsetC1+numsetC2);
    end
        
        % Check if any UY is among selected basis function and calc basis
        % recursively
        UYpresent = find(ismember(f_indx,UYindex),1);
        ybml(:,1:2) = data(:,1:2);
        if ~isempty(UYpresent)
            % disp('UY present')
            %  Utransi=zeros(Ntrans,lenreading);
            for ntr=2:lenreadpone
                [Utransii,~,~] = InputTransformation(ybml(:,ntr-1),u(:,ntr-1),iUY,i_UiUj,ntr);
                Utransii = Utransii(idOb,:);
                Utransi(:,ntr-1) = Utransii(1:nut,:);
                ybml(:,ntr) =  A*ybml(:,ntr-1) + C1*u(:,ntr-1) + C2*Utransi(:,ntr-1) + D;
            end
        else
            for i = 2:lenreadpone
                ybml(:,i) =  A*ybml(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
            end
        end
  
%% Objective function calculation
J = hwhole;
clear hwhole
size(J)
rJ=rank(J);

if rJ>=Lenvec
  epss = data-ybml;
    epssTepss = epss.*epss;
    E_R_E = (1/sgma)*sum(epssTepss,"all");                  % R^-1 = 1/sigma * eye(N) so we cann simply bring sigma out
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                     % Determinant of diagonal matrix=product of the diagonal  elements

    L_theta = -0.5*N*lenreadpone*log(2*pi) - 0.5*lenreadpone*log(det_R) - 0.5*E_R_E;

    JtJu = zeros(Lenvec);
    for ir=1:Lenvec
        j=ir:Lenvec;
        for ic=j
            JtJu(ir,ic) = sum(J(:,ir).*J(:,ic));
        end
    end
    JtJ = JtJu + triu(JtJu,1)';
    Co_theta = sgma*inv(JtJ);
    C_ve = 0.5*sum(log(diag(Co_theta))) - 0.5*log(det(Co_theta));

    % AICc + 2Cve
    Obj_fun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreadpone-Lenvec-1) + 2*C_ve
else
    Obj_fun=inf                                                  % This is for a situation we cannot obtain full rank Jacobian for selected sunset
end

A_mat = reshape(A',1,N2);
C1_mat = reshape(C1',1,numel(C1));
C2_mat = reshape(C2',1,numel(C2));
D_mat = D';
l_mat = lam';

end
