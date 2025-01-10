function [Objfun, Cmat, Dmat] = RankBasis(indx)

load('data4rank.mat','data_rank','utrans_rank')
data = data_rank;
Utrans = utrans_rank;
clear data_rank utrans_rank
f_indx = find(indx);

%% Estimation and model predictions
sgma = 0.06/5;   
tic
    %% define the parameters
    N = size(data,1);
    N2 = N^2;                                              
    Ntrans = size(Utrans,1);                                
    ntrans = numel(f_indx);                                 

    %% Indices for parameters in matrices
    setC=[];
    for stc = 1:N
        setC = [setC f_indx+(stc-1)*Ntrans];
    end

    numsetC = numel(setC);         
    ar_c = zeros(numsetC,1);
    ac_c = ar_c;

    setD = [];                                          % Algorithm I does not require constant basis function to be present
    numsetD = numel(setD);
    ar_d = (1:N)';
    ac_d = ones(N,1);

    lenreading = size(data,2);                
    Lenvec = numsetC + numsetD;                          % No of model parameters            
    sz = lenreading*N;             
    szvec = sz+Lenvec;              
    vecparams = zeros(Lenvec,1); 
    szlam = N2 + (N-N2)/2;                               % Number of covariance
    meanprior = vecparams;
    
    %% define C and D and obtain vecparams for reduced system
    C = zeros(N,Ntrans);   
    for i = 1:numsetC
        ar_c(i) = ceil((setC(i))/Ntrans);
        ac_c(i) = setC(i) - Ntrans*fix((setC(i))/Ntrans); 
        if ac_c(i) == 0
            ac_c(i) = Ntrans;
        end
        C(ar_c(i),ac_c(i)) = rand;
        vecparams(i) = C(ar_c(i),ac_c(i));
    end

    D = zeros(N,1);
    for j = 1:numsetD
        D(ar_d(j),ac_d(j)) = rand;
        vecparams(j+numsetC) = D(ar_d(j),ac_d(j));
    end

    %% Monitors        
    Trend = zeros(Lenvec,30); 
    firstorder = zeros(10,1); 
    
    %% Initialize xq
    y = data(1:N,:);
    for i = 1:lenreading
        y(:,i) = C*Utrans(:,i) + D; 
    end
    
    %% M-step variables
    lam =  rand (szlam,1); 
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
   
    %% Jacobian remains unchanged with model parameters since model is linear in parameters
    dCC = zeros(sz,ntrans);   
    dDD = zeros(sz,N);
    for ij = 1:N
        dCC((ij-1)*lenreading + 1:ij*lenreading,(ij-1)*ntrans+1:ij*ntrans) = Utrans(f_indx,1:end)';
        dDD((ij-1)*lenreading + 1:ij*lenreading,ij) = ones(lenreading,1);
    end

    if numsetD==0 
        hwhole = dCC;
    else
        hwhole = [dCC dDD];                 
    end
    clear dCC dDD
    
    % Checking for estimability
    if rank(hwhole)<Lenvec                           
        Cmat = zeros(1,N*Ntrans);
        Dmat = zeros(1,N);
        Objfun = inf;
        return
    end
    tryy = 1;

    %% Loop
    while conv > tol
        n=n+1;
        n1 = n+1;
    
        %% ybar,Jbar
        for i=1:lenreading
            yp = C*Utrans(:,i) + D;
            df = data(:,i)- yp; 
            DF(i:lenreading:(N-1)*lenreading +i) = df(1:end); 
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
            updatelam = 0
        else
            updatelam = secder\der;
        end
        lam = lam - 0.6*updatelam;
    
        %% Update parameters
        % updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
        updateparams = Cth_y*JCey;
        oldparams = vecparams; 
        vecparams = oldparams + updateparams;
        
        %% Update C, D
        for j = 1:numsetC
            ar = ar_c(j);
            ac = ac_c(j);
            C(ar,ac) = vecparams(j);
        end

        for k = 1:numsetD
            dr = ar_d(k);
            dc = ac_d(k);
            D(dr,dc) = vecparams(k+numsetC);
        end

        for i = 1:lenreading
            y(:,i) = C*Utrans(:,i) + D; 
        end

        diff(n1) = norm(data - y); 
        %% Search for best parameters \convergence
        halv=0;
        while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)           
    %         disp('halve')
            halv = halv + 1;
            updateparams = 0.6*updateparams;
            vecparams = oldparams + updateparams;
            
            for j = 1:numsetC
                ar = ar_c(j);
                ac = ac_c(j);
                C(ar,ac) = vecparams(j);
            end
    
            for k = 1:numsetD
                dr = ar_d(k);
                dc = ac_d(k);
                D(dr,dc) = vecparams(k+numsetC);
            end

            for i = 1:lenreading
                y(:,i) = C*Utrans(:,i) + D; 
            end
            diff(n1) = norm(data - y);
        end
    
        if diff(n1)>diff(n)  && tryy<3               % Trying another starting point and see how it behaves   
            tryy = tryy +1
            vecparams = rand(Lenvec,1)*tryy^2;
            for j = 1:numsetC
                ar = ar_c(j);
                ac = ac_c(j);
                C(ar,ac) = vecparams(j);
            end
            for k = 1:numsetD
                dr = ar_d(k);
                dc = ac_d(k);
                D(dr,dc) = vecparams(k+numsetC);
            end
            lam=rand(szlam,1); n=0;
            diff = ones(30,1)*Inf; conv=1;
            
        elseif diff(n1)>diff(n)
            Cmat = zeros(1,N*Ntrans);
            Dmat = zeros(1,N);
            Objfun = inf;
            return
        end
    
        if n>0
            Trend(:,n) = vecparams;
            firstorder(n) = norm(hwhole,inf);  
            % Check conditions for proper exit
            conv = rms(updateparams);
        end
        
    end
    
    % t=1:size(data,2);
    % subplot(N,1,1), plot(t,data(1,:),'k-',t,y(1,:),'b:')

    toc

%% Objective function calculation
J = hwhole;
clear hwhole
rJ=rank(J);

if rJ>=Lenvec
    epss = data-y;
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
    Objfun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreading-Lenvec-1) + 2*C_ve
else
    Objfun=inf                                                 
end

Cmat = reshape(C',1,numel(C));
Dmat = D';

end
