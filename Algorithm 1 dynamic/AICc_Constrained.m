function [Obj_fun, A_mat, C_mat, l_mat, ybml_trend, yddr_trend, w_trend, fval_trend] = AICc_Constrained(indx, iUY, i_UiUj, UYindex, Ncons, Ain, Cin, Lin,nci,nhio)

load('data4est.mat','data_est','u_est', 'sel_indx', 'utrans_est','idOb','nut','mnd','mxd')
data = data_est;
u = u_est;
Utrans = utrans_est;
clear data_est u_est utrans_est
f_indx = find(indx);

%% Estimation and model predictions
sgma = 0.06;   
tic
    %% define the parameters
    Nwmcon = size(data,1);
    N = Nwmcon - Ncons;                                     
    N2 = N^2;                                               
    Ntrans = size(Utrans,1);                                
    ntrans = numel(f_indx);
    
    AtomIN = data(N+1:end,:);
    data = data(1:N,:);

    %% Indices for parameters in matrices
    setA= 1:N2;       
    setC=[];
    for stc = 1:N
        setC = [setC f_indx+(stc-1)*Ntrans];
    end
    
    numsetA = numel(setA);         
    ar_a = zeros(numsetA,1);
    ac_a = ar_a;
    numsetC = numel(setC);         
    ar_c = zeros(numsetC,1);
    ac_c = ar_c;

    count1 = 0;
    lenreading = size(data,2)-1;               
    lenreadpone = lenreading+1;     
    Lenvec = numsetA + numsetC;                
    sz = lenreading*N;                       
    sz1 = lenreading*Nwmcon;
    szvec = sz+Lenvec;                         
    vecparams = zeros(Lenvec,1); 
    szlam = N^2 + (N-N^2)/2;                    
    meanprior = vecparams;
    
    %% define A and C and obtain vecparams for reduced system
        A = reshape(Ain,N,N)';
        C = reshape(Cin,[Ntrans,N])';
        lam = Lin';

    for i = 1:numsetA
        ar_a(i) = ceil(setA(i)/N);
        ac_a(i) = N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N);
        vecparams(i) = A(ar_a(i),ac_a(i));
    end
    for i = 1:numsetC
        ar_c(i) = ceil((setC(i))/Ntrans);
        ac_c(i) = setC(i) - Ntrans*fix((setC(i))/Ntrans); 
        if ac_c(i) == 0
            ac_c(i) = Ntrans;
        end
        vecparams(i+numsetA) = C(ar_c(i),ac_c(i));
    end

    
    %% Monitors        
    Trend = zeros(Lenvec,30); 
    firstorder = zeros(10,1); 
    
    %% Initialize xq
    ybml = data(1:N,:);
    for i = 1:lenreading
        ybml(:,i+1) = A*ybml(:,i) + C*Utrans(:,i) ; 
    end
    ybmls = ybml;
    datas = data;
    for i=1:N 
        ybmls(i,:) = mnd(i) + 0.5*(ybml(i,:)-1)*(mxd(i)-mnd(i));
        datas(i,:) = mnd(i) + 0.5*(data(i,:)-1)*(mxd(i)-mnd(i));
    end


    %% M-step variables
    der = lam;                 
    secder = zeros(szlam);   
    cth = 20;
    
    %% Initialize variables for J matrix
    hwhole = zeros(sz,Lenvec); %h_whole=hwhole;
    DF = hwhole(:,1);
    %% Loop utilities
    tol = 1e-4; 
    conv =1;
    n = 0;
    %% Loop Flags
    diff = ones(30,1)*Inf;

    % Exploiting the linearity of model in parameters to form Jacobian from 
    % (transformed) data to minimise computational burden. (time difference is
    % in the order of 10e3!!!)
    
    dAA = zeros(sz,N);
    dCC = zeros(sz,ntrans);
    for ij = 1:N
        dAA((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*N+1:ij*N) = data(:,1:end-1)';
        dCC((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*ntrans+1:ij*ntrans) = Utrans(f_indx,1:end-1)';
    end

    hwhole = [dAA dCC];

    % Checking for estimability: If Jacobian is rank defficient, we do not ...
    % need to carry out estimation...skip the current candidate model
    % Jcbn = hwhole;
    if rank(hwhole)<Lenvec                           
        A_mat = zeros(1,N2);
        C_mat = zeros(1,N*Ntrans);
        l_mat = zeros(1,szlam);
        Obj_fun = inf;
        ybml_trend = zeros(sz,10);
        yddr_trend = zeros(sz,10);
        w_trend = zeros(N,10);
        return
    end
    tryy = 1;
    
    %% 
        len = size(sel_indx,2);         % No. of steady state points
        objfun = @(x)ddr_obj(x, datas(:,2:end), ybmls(:,2:end), N, lenreading);
        nlcon = @(x)ddr_nlcon(x, AtomIN, N, Ncons, sel_indx, len, mnd, mxd,nci,nhio);
    
        nlrhs = [zeros(Ncons*len,1)]; 
        nle = zeros(Ncons*len,1); 
        
        lb = [0.9*ones(N,1); -1e5*ones(sz,1)]; 
        ub = [1.1*ones(N,1); 1e5*ones(sz,1)]; 

        w0 = ones(N,1);
        yddr0 = reshape(ybml(:,2:end),sz,1);
        yddr_w0 = [w0; yddr0];

        optns = optiset('solver', 'ipopt', 'display', 'iter', 'maxiter', 1e2, 'maxtime', 1000);
        Optim = opti('fun', objfun, 'ineq', [], [], 'nlmix', nlcon, nlrhs, nle, 'bounds', lb, ub, 'options', optns)
        
        [yddr_w,fval,exitflag,info] = solve(Optim,yddr_w0);
    
    %% Loops

    ybml_trend = zeros(sz,10);
    yddr_trend = ybml_trend;
    fval_trend = zeros(10,1);
    w_trend = zeros(N,10);

    ybml_trend(:,1) = reshape(ybmls(:,2:end),sz,1);
    yddr_trend(:,1) = yddr_w(N+1:N+sz);
    fval_trend(1) = fval;
    w_trend(:,1) = yddr_w(1:N);

    %% Main Estimation Loop
    conv_outer = 1;
    nddr = 0;
    flag_obj = true;
    while conv_outer>1e-6 && flag_obj==true && nddr<8  
        nddr = nddr+1;

        n = 0;
        conv2 =1;

        yddrs = reshape(yddr_trend(:,nddr),N,lenreading);
        yxs = ybmls./w_trend(:,nddr);

        yddr = yddrs;
        yx = yxs;
        for i=1:N
            yddr(i,:)=1 + 2*(yddrs(i,:)-mnd(i))/(mxd(i)-mnd(i));
            yx(i,:)=1 + 2*(yxs(i,:)-mnd(i))/(mxd(i)-mnd(i));
            ybml(i,:)=1 + 2*(ybmls(i,:)-mnd(i))/(mxd(i)-mnd(i));
        end

        while conv2 > tol && n<80 && flag_obj
            n=n+1;
            n1 = n+1;
        
            for i=2:lenreadpone
                ybml(:,i) = A*data(:,i-1) + C*Utrans(:,i-1);
                df = yx(:,i)- ybml(:,i); 
                DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(1:end); % em for missing data could come here
            end
            m_v = meanprior - vecparams;
            ybar = [DF; m_v];
            Jbar = [hwhole; eye(Lenvec)];
    
            %% M step
            ceph_u = zeros(N); Ce_invu = zeros(sz);
            lamn = 0;
            for ii=1:N
                jj=ii:N;
                for kk=jj
                    lamn=lamn+1;
                    ceph_u(ii,kk) = lam(lamn);
                end
            end
            ceph = ceph_u + triu(ceph_u,1)';
    
            ceph_inv = ceph\eye(N);
            % Expand inverse of unit covariance matrix for lenreading points
            lamn = 0;
            for ii=1:N
                jj=ii:N;
                for kk=jj
                    lamn=lamn+1;
                    Ce_invu((ii-1)*lenreading+1:ii*lenreading,(kk-1)*lenreading+1:kk*lenreading) = ceph_inv(ii,kk)*eye(lenreading);
                end
            end
            Ce_inv = sparse(Ce_invu + triu(Ce_invu,1)');  
    
            clear Ce_invu
        
            Cth_inv = (1/cth)*eye(Lenvec);
            Ce_y1 = Ce_inv*DF;
            Ce_y2 = Cth_inv*m_v;
            Ce_y = [Ce_y1; Ce_y2];
            JCey = Jbar'*Ce_y;
            clear Ce_y
    
            Ce_J = [Ce_inv*hwhole; Cth_inv];
            JCeJ = Jbar'*Ce_J;
            Cth_y = JCeJ\eye(Lenvec);
            clear JCeJ
    
            Cebar_inv = sparse([Ce_inv zeros(sz,Lenvec); zeros(Lenvec,sz)  Cth_inv]);
            JCe = Jbar'*Cebar_inv;
            CeJCth = Ce_J*Cth_y;
            CeJCthJCe = CeJCth*JCe;         
            clear CeJCth Ce_J Ce_inv
    
            P = Cebar_inv - CeJCthJCe;
            clear Cebar_inv
    
            %% 
            lamn=0;
            for ii=1:N
                jj=ii:N;
                for kk=jj
                    vpy=zeros(szvec,1);
                    vpy((ii-1)*lenreading+1:ii*lenreading) = P((kk-1)*lenreading+1:kk*lenreading,:)*ybar;
                    % Trace of pv is the trace of nonzero square matrix that falls on the diagonal
                    tr_pv = sum(diag(P((kk-1)*lenreading+1:kk*lenreading,(ii-1)*lenreading+1:ii*lenreading)));
                    
                    lamn=lamn+1;
                    der(lamn) = -0.5*tr_pv + 0.5*ybar'*P'*vpy;
            
                    rc(lamn,:) = [ii kk];                   % row-col indices for lam          
            
                end
            end
    
            for sdr=1:szlam                             % Secder row
                sdci=sdr:szlam;                         % Secder column indices
                for sdc=sdci
                    rwi = rc(sdr,:);                    % row-col indices of vi
                    cli = rc(sdc,:);                    % row-col indices of vj
    
                    % For pv/pvi/pvj, only selecting elements of P that will
                    % eventually matter for obtaining trace(PViPVj) as P is
                    % multiplied with each v/vi/vj
    
                    if sdr==sdc                         % pvipvj for principal diagonal
                        pv = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading,(rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                        pvpv = pv.*pv';
                        clear pv
                        secder(sdr,sdc) = -0.5*sum(pvpv,'all');
                    elseif rwi(1)==cli(1)               % If vi & vj are in same row
                        % pvi calculates p * v that corresponds to present ROW OF SECDER
                        pvi = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading,(rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                        % pvj calculates p * v that corresponds to present COLUMN OF SECDER
                        pvj = P((cli(2)-1)*lenreading+1:cli(2)*lenreading,(cli(1)-1)*lenreading+1:cli(1)*lenreading);
                        pvipvj = pvi.*pvj';
                        clear pvi pvj
                        secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
                        clear pvipvj
                        secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
                    else
                        % pvi calculates p * v that corresponds to present ROW OF SECDER
                        pvi = P((cli(2)-1)*lenreading+1:cli(2)*lenreading,(rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                        % pvj calculates p * v that corresponds to present COLUMN OF SECDER
                        pvj = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading,(cli(1)-1)*lenreading+1:cli(1)*lenreading);
                        pvipvj = pvi.*pvj';
                        clear pvi pvj
                        secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
                        clear pvipvj
                        secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
                        
                    end
                end
            end
    
        %%
            secderr = rcond(secder);       
            if secderr <1e-13             
                updatelam = 0;
            else
                updatelam = secder\der;
            end
            conv2 = rms(updatelam);
            lam = lam - 0.6*updatelam;
    
            %% Update parameters
            updateparams = Cth_y*JCey;
            oldparams = vecparams; %?
            vecparams = oldparams + updateparams;

            %% Update A,C
            for j = 1:numsetA
                ar = ar_a(j);
                ac = ac_a(j);
                A(ar,ac) = vecparams(j);
            end
            
            for j = 1:numsetC
                ar = ar_c(j);
                ac = ac_c(j);
                C(ar,ac) = vecparams(j+numsetA);
            end
    
            % Check if any UY is among selected basis function and calc basis
            % recursively
            UYpresent = find(ismember(f_indx,UYindex),1);
            if ~isempty(UYpresent)
                % disp('UY present')
                ybml(:,1) = data(:,1); 
                for ntr=1:lenreading
                    [Utransii,~,~] = InputTransformation(ybml(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                    Utransii = Utransii(idOb,:);
                    Utransi(:,ntr) = Utransii(1:nut,:);
                    ybml(:,ntr+1) = A*ybml(:,ntr) + C*Utransi(:,ntr); 
                end
            else
                for i = 1:lenreading
                    ybml(:,i+1) = A*ybml(:,i) + C*Utrans(:,i);              % Model Predictions
                end
            end
    
            if Ncons>0
                for i=1:N 
                    ybmls(i,:) = mnd(i) + 0.5*(ybml(i,:)-1)*(mxd(i)-mnd(i));
                end
                for i = 1:lenreadpone
                    cbalhats(i) = nci*ybmls(:,i);        % Predicted carbon in outlet
                    hbalhats(i) = nhio*ybmls(:,i);       % Predicted net H on RHS
                end
                ybml_wmcon = [ybml;cbalhats;hbalhats];   % Prediction with mass constraint
            else
                ybml_wmcon = ybml;
            end
            
            diff(n1) = norm(yx - ybml); 
            %% Search for best parameters \convergence
            halv=0;
            while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)          
                halv = halv + 1;
                updateparams = 0.6*updateparams;
                vecparams = oldparams + vecparams;
                
                for j = 1:numsetA
                    ar = ar_a(j);
                    ac = ac_a(j);
                    A(ar,ac) = vecparams(j);
                end
                       
                for j = 1:numsetC
                    ar = ar_c(j);
                    ac = ac_c(j);
                    C(ar,ac) = vecparams(j+numsetA);
                end
    
                % Check if any UY is among selected basis function and calc basis
                % recursively
                UYpresent = find(ismember(f_indx,UYindex),1);
                if ~isempty(UYpresent)
                    ybml(:,1) = data(:,1);
                    for ntr=1:lenreading
                        [Utransii,~,~] = InputTransformation(ybml(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                        Utransii = Utransii(idOb,:);
                        Utransi(:,ntr) = Utransii(1:nut,:);
                        ybml(:,ntr+1) = A*ybml(:,ntr) + C*Utransi(:,ntr); 
                     end
                else
                    for i = 1:lenreading
                        ybml(:,i+1) = A*ybml(:,i) + C*Utrans(:,i);              
                    end
                end

                if Ncons>0
                    for i=1:N 
                        ybmls(i,:) = mnd(i) + 0.5*(ybml(i,:)-1)*(mxd(i)-mnd(i));
                    end
                    for i = 1:lenreadpone
                        cbalhats(i) = nci*ybmls(:,i);        % Predicted carbon in outlet
                        hbalhats(i) = nhio*ybmls(:,i);       % Predicted net H on RHS
                    end
                    ybml_wmcon = [ybml;cbalhats;hbalhats];       % Prediction with mass constraint
                else
                    ybml_wmcon = ybml;
                end
                diff(n1) = norm(yx - ybml);
            end
    
            if diff(n1)>diff(n)  && tryy<3               % Trying another starting point and see how it behaves   
                tryy = tryy +1
                vecparams = rand(Lenvec,1)*tryy^2;
                for j = 1:numsetA
                    ar = ar_a(j);
                    ac = ac_a(j);
                    A(ar,ac) = vecparams(j);
                end
                for j = 1:numsetC
                    ar = ar_c(j);
                    ac = ac_c(j);
                    C(ar,ac) = vecparams(j+numsetA);
                end
                lam=rand(szlam,1); n=0;
                diff = ones(30,1)*Inf; conv2=1;
                
            elseif diff(n1)>diff(n)
                A_mat = zeros(1,N2);
                C_mat = zeros(1,N*Ntrans);
                l_mat = zeros(1,szlam);
                Obj_fun = inf;
                ybml_trend = zeros(sz,10);
                yddr_trend = zeros(sz,10);
                % fval_trend = inf;
                w_trend = zeros(N,10);
                return
            end
        
            updateparams = vecparams - oldparams;
            if n>0
                Trend(:,n) = vecparams;
                firstorder(n) = norm(hwhole,inf);  
                %% Check conditions for proper exit
                conv2 = rms(updateparams);
            end
            
        end

        objfun = @(x)ddr_obj(x, datas(:,2:end), ybmls(:,2:end), N, lenreading);
        nlcon = @(x)ddr_nlcon(x, AtomIN, N, Ncons, sel_indx, len, mnd, mxd,nci,nhio);
    
        nlrhs = [zeros(Ncons*len,1)]; 
        nle = zeros(Ncons*len,1); 
        
        lb = [0.9*ones(N,1); -1e5*ones(sz,1)]; 
        ub = [1.1*ones(N,1); 1e5*ones(sz,1)]; 

        w0 = w_trend(:,nddr);
        yddr0 = yddr_trend(:,nddr);
        yddr_w0 = [w0; yddr0];

        optns = optiset('solver', 'ipopt', 'display', 'iter', 'maxiter', 1e2, 'maxtime', 1000);
        Optim = opti('fun', objfun, 'ineq', [], [], 'nlmix', nlcon, nlrhs, nle, 'bounds', lb, ub, 'options', optns)
        
        [yddr_w1,fval1,exitflag,info] = solve(Optim,yddr_w0);
    
        flag_obj = fval1<fval_trend(nddr);

        if flag_obj
            ybml_trend(:,nddr+1) = reshape(ybmls(:,2:end),sz,1);
            yddr_trend(:,nddr+1) = yddr_w1(N+1:N+sz);
            fval_trend(nddr+1) = fval1;
            w_trend(:,nddr+1) = yddr_w1(1:N);
        end
        conv_outer = rms(w_trend(:,nddr) - w_trend(:,nddr+1));
    end
    t=1:size(ybml,2);
    subplot(N,1,1), plot(t,yx(1,:),'k-',t,ybml(1,:),'b-.')

%% Jacobian and Likelihood fcn
J = hwhole;
clear hwhole
size(J)
rJ=rank(J);

if rJ>=Lenvec
    epss = data-ybml;
    epssTepss = epss.*epss;
    E_R_E = (1/sgma)*sum(epssTepss,"all");                
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                    

    L_theta = -0.5*N*lenreading*log(2*pi) - 0.5*lenreading*log(det_R) - 0.5*E_R_E;
   
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
    Obj_fun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreading-Lenvec-1) + 2*C_ve
else
    Obj_fun=inf                                                  % This is for a situation we cannot obtain full rank Jacobian for selected sunset
end

A_mat = reshape(A',1,numel(A));
C_mat = reshape(C',1,numel(C));
l_mat = lam';

% fprintf('tryy                       %d \n',tryy)
% fprintf('n                          %d \n',n)
% fprintf('halv                       %d \n',halv)
% % fprintf('UY bases                   %d \n',UYpresent)
% fprintf('Cve                        %d \n',C_ve)
% fprintf('Loglikelihood              %d \n',L_theta)
% fprintf('Objfun                     %d \n',Obj_fun)
% fprintf('estimation time            %2f \n',toc)

end
