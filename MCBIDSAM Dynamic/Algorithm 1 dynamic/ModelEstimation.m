function [Objfun, Amat, Cmat, Lmat] = ModelEstimation(indx, iUY, i_UiUj, UYindex, Ncons, Ain, Cin, Lin)

load('data4est.mat','data_est','u_est','utrans_est','idOb','nut')
data = data_est;
u = u_est;
Utrans = utrans_est;
clear data_est u_est utrans_est
f_indx = find(indx);

%% Estimate model and its predictions
sgma = 0.06; 

tic
    %% define the parameters
    Nwmcon = size(data,1);
    N = Nwmcon - Ncons;                                     % NUMBER OF OUTPUT
    N2 = N^2;                                               % NUMBER OF ELEMENTS IN A
    Ntrans = size(Utrans,1);                                % NUMBER OF TRANSFORMATIONs
    ntrans=numel(f_indx);
    
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
    
    lenreading = size(data,2)-1;                            % NUMBER OF DATA POINTS/SAMPLES EXCLUDE INTIAL
    lenreadpone = lenreading+1;                             % NUMBER OF DATA POINTS
    Lenvec = numsetA+numsetC;                               % NUMBER OF PARAMETERS
    sz = lenreading*N;                                      % NUMBER OF DATA READINGS
    szvec = sz+Lenvec;                                      % DIMENSION OF LENVEC
    vecparams = zeros(Lenvec,1); 
    szlam = N2 + (N-N2)/2;                                  % Number of covariance
    meanprior = vecparams;
    
    %% define A, B and lam and obtain vecparams for reduced system
    nulc = isequal(Ain,zeros(1,N2));
    if ~nulc
        A = reshape(Ain,[N,N])';
        C = reshape(Cin,[Ntrans,N])';
        lam = Lin';
    else
        A = zeros(N);
        C = zeros(N,Ntrans);
        lam = rand (szlam,1);
    end
    for i = 1:numsetA
        ar_a(i) = ceil(setA(i)/N);
        ac_a(i) = N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N);
        if nulc
            A(ar_a(i),ac_a(i)) = rand;
        end
        vecparams(i) = A(ar_a(i),ac_a(i));
    end

    for i = 1:numsetC
        ar_c(i) = ceil((setC(i))/Ntrans);
        ac_c(i) = setC(i) - Ntrans*fix((setC(i))/Ntrans); 
        if ac_c(i) == 0
            ac_c(i) = Ntrans;
        end
        if nulc
            C(ar_c(i),ac_c(i)) = rand;
        end
        vecparams(i+numsetA) = C(ar_c(i),ac_c(i));
    end
    
    Trend = zeros(Lenvec,30);
    %% Initialize xq
    y = data;
    
    for i = 1:lenreading
        y(:,i+1) = A*y(:,i) + C*Utrans(:,i) ; 
    end
    
    %% M-step variables
    der = lam;                                              % DERIVATIVE OF F
    secder = zeros(szlam);                                  % SECOND DERIVATIVE OF F
    cth = 20;
    
    %% Initialize variables for J matrix
    hwhole = zeros(sz,Lenvec); %h_whole=hwhole;
    DF = hwhole(:,1);
    
    %% Loop utilities
    tol = 2e-3;
    conv =1;
    n = 0;
    
    %% Loop Flags
    diff = ones(30,1)*Inf;

    %% Exploiting the linearity of model in parameters to form Jacobian from 
    dAA = zeros(sz,N);
    dCC = zeros(sz,ntrans);
    for ij = 1:N
        dAA((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*N+1:ij*N) = data(:,1:end-1)';
        dCC((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*ntrans+1:ij*ntrans) = Utrans(f_indx,1:end-1)';
    end
    
    hwhole = [dAA dCC];
    
    % Checking for estimability: If Jacobian is rank defficient, we do not ...
    % need to carry out estimation...skip the current candidate model
    if rank(hwhole)<Lenvec                           
        Amat = zeros(1,N2);
        Cmat = zeros(1,N*Ntrans);
        Lmat = zeros(1,szlam);
        return
    end
    tryy = 1;

    %% Loop
    while conv > tol && n<50
        n=n+1;
        n1 = n+1;
    
        %% ybar,Jbar
        for i=2:lenreadpone
            df = data(:,i)- (A*data(:,i-1) + C*Utrans(:,i-1)); % data row is y, col is time
            DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(1:end); 
        end
        m_v = meanprior - vecparams;
        ybar = [DF; m_v];
        Jbar = [hwhole; eye(Lenvec)];

        
        %% M step
        ceph_u = zeros(N); Ce_invu = zeros(sz);
        % Construct unit covariance matrix
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
        % Expand inverse of unit covarianc matrix for lenreading points
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
        
        secderr = rcond(secder);        
        if secderr <1e-13               
            updatelam = 0
            %             break
        else
            updatelam = secder\der;
        end
        lam = lam - 0.6*updatelam;
    
        %% Update parameters
        updateparams = Cth_y*JCey;
        oldparams = vecparams; %?
        vecparams = oldparams + updateparams;

        %% Update A,C. 
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
            y(:,1) = data(:,1); 
            for ntr=1:lenreading
                [Utransii,~,~] = InputTransformation(y(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                Utransii = Utransii(idOb,:);
                Utransi(:,ntr) = Utransii(1:nut,:);
                y(:,ntr+1) = A*y(:,ntr) + C*Utransi(:,ntr); 
            end
        else
            for i = 1:lenreading
                y(:,i+1) = A*y(:,i) + C*Utrans(:,i);              % Model Predictions
            end
        end

        diff(n1) = norm(data(:,2:end) - y(:,2:end)); 
        %% Search for best parameters \convergence
        halv = 0;
        while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)        
            halv = halv + 1;
            updateparams = 0.6*updateparams;
            vecparams = oldparams + updateparams;

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
                y(:,1) = data(:,1);
                % Utransi=zeros(Ntrans,lenreading);
                for ntr=1:lenreading
                    [Utransii,~,~] = InputTransformation(y(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                    Utransii = Utransii(idOb,:);
                    Utransi(:,ntr) = Utransii(1:nut,:);
                    y(:,ntr+1) = A*y(:,ntr) + C*Utransi(:,ntr); 
                 end
            else
                for i = 1:lenreading
                    y(:,i+1) = A*y(:,i) + C*Utrans(:,i);              
                end
            end
        %     
            diff(n1) = norm(data(:,2:end) - y(:,2:end));

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
            diff = ones(30,1)*Inf; conv=1;

        elseif diff(n1)>diff(n)
            Amat = zeros(1,N2);
            Cmat = zeros(1,N*Ntrans);
            Lmat = zeros(1,szlam);
            Objfun = inf;
            return
        end
        if n>0
            Trend(:,n) = vecparams;
            % Check conditions for proper exit
            conv = rms(updateparams);
        end
    end
    
    toc   
%% Jacobian and Likelihood fcn
J = hwhole;
clear hwhole
size(J)

if rank(J)>=Lenvec
    epss = data-y;
    epssTepss = epss.*epss;
    E_R_E = sgma\sum(epssTepss,"all");                  % R^-1 = 1/sigma * eye(N) so we cann simply bring sigma out
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                     % Determinant of diagonal matrix=product of the diagonal  elements
    L_theta = -0.5*N*lenreading*log(2*pi) - 0.5*lenreading*log(det_R) - 0.5*E_R_E;
   ttoc(16) = toc;

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
    Objfun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreading-Lenvec-1) + 2*C_ve
    if isnan(Objfun)
        Objfun=inf
    end
else
    Objfun=inf                       
end
Amat = reshape(A',1,numel(A));
Cmat = reshape(C',1,numel(C));
Lmat = lam';
end