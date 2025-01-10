function [Objfun, Cmat, Dmat, lmat] = ModelEstimation(indx, Ncons, Cin, Din, Lin)

load('data4est.mat','data_est','u_est','utrans_est')
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
    
    data = data(1:N,:);

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

    lenreading = size(data,2)-1;    
    lenreadpone = lenreading+1;     
    Lenvec = numsetC + numsetD;                         % No. of parameters        
    sz = lenreadpone*N;              
    szvec = sz+Lenvec;              
    vecparams = zeros(Lenvec,1); 
    szlam = N2 + (N-N2)/2;          
    meanprior = vecparams;
    
    %% define C and D and obtain vecparams for reduced system
    nulc =1; isequal(Cin,zeros(1,N*Ntrans));
    if ~nulc
        C = reshape(Cin,[Ntrans,N])';
        D = Din';
        lam = Lin';
    else
        C = zeros(N,Ntrans);
        D = zeros (N,1);
        lam = rand(szlam,1);
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
        vecparams(i) = C(ar_c(i),ac_c(i));
    end

    for i = 1:numsetD
        if nulc
            D(ar_d(i),ac_d(i)) = rand;
        end
        vecparams(i+numsetC) = D(ar_d(i),ac_d(i));
    end
    
    %% Monitors        
    Trend = zeros(Lenvec,30);
    firstorder = zeros(10,1); 

    %% Initialize xq
    y = data(1:N,:);
    for i = 1:lenreadpone
        y(:,i) = C*Utrans(:,i) + D; 
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

    %% Jacobian remains unchanged with model parameters since model is linear in parameters
    dCC = zeros(sz,ntrans);   dDD = zeros(lenreadpone,N);
    for ij = 1:N
        dCC((ij-1)*lenreadpone + 1:ij*lenreadpone,(ij-1)*ntrans+1:ij*ntrans) = Utrans(f_indx,1:end)';
        dDD((ij-1)*lenreadpone + 1:ij*lenreadpone,ij) = ones(lenreadpone,1);
    end

    if numsetD==0 
        hwhole = dCC;
    else
        hwhole = [dCC dDD];            
    end
    
    % Checking for estimability
    if rank(hwhole)<Lenvec                           
        Cmat = zeros(1,N*Ntrans);
        Dmat = zeros(1,N);
        lmat = zeros(1,szlam);
        Objfun = inf;
        return
    end
    tryy = 1;
    
    %% Loops
    while conv > tol && n<80
        n=n+1;
        n1 = n+1;
    
        %% ybar,Jbar
        for i=1:lenreadpone
            yp = C*Utrans(:,i) + D;
            df = data(:,i)- yp; 
            DF(i:lenreadpone:(N-1)*lenreadpone +i) = df(1:end); 
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
                Ce_invu((ii-1)*lenreadpone+1:ii*lenreadpone,(kk-1)*lenreadpone+1:kk*lenreadpone) = ceph_inv(ii,kk)*eye(lenreadpone);
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
                vpy((ii-1)*lenreadpone+1:ii*lenreadpone) = P((kk-1)*lenreadpone+1:kk*lenreadpone,:)*ybar;
                % Trace of pv is the trace of nonzero square matrix that falls on the diagonal
                tr_pv = sum(diag(P((kk-1)*lenreadpone+1:kk*lenreadpone,(ii-1)*lenreadpone+1:ii*lenreadpone)));
                
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
                    pv = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone, (rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
                    pvpv = pv.*pv';
                    clear pv
                    secder(sdr,sdc) = -0.5 * sum(pvpv, 'all');
                elseif rwi(1)==cli(1)               % If vi & vj are in same row
                    % pvi calculates p * v that corresponds to present ROW OF SECDER
                    pvi = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone, (rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
                    % pvj calculates p * v that corresponds to present COLUMN OF SECDER
                    pvj = P((cli(2)-1)*lenreadpone+1:cli(2)*lenreadpone, (cli(1)-1)*lenreadpone+1:cli(1)*lenreadpone);
                    pvipvj = pvi.*pvj';
                    clear pvi pvj
                    secder(sdr,sdc) = -0.5*sum(pvipvj, "all");
                    clear pvipvj
                    secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
                else
                    % pvi calculates p * v that corresponds to present ROW OF SECDER
                    pvi = P((cli(2)-1)*lenreadpone+1:cli(2)*lenreadpone,(rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
                    % pvj calculates p * v that corresponds to present COLUMN OF SECDER
                    pvj = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone,(cli(1)-1)*lenreadpone+1:cli(1)*lenreadpone);
                    pvipvj = pvi.*pvj';
                    clear pvi pvj
                    secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
                    clear pvipvj
                    secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
                    
                end
            end
        end

    %%
        secderr = rcond(secder);        %CHECK CONDITIONALITY OF SECOND DERIVATIVE BEFORE INVERSION
        if secderr <1e-13               %IF POORLY CONDITIONED, DON'T UPDATE LAMDA IN PRESENT ITERATION
            updatelam = 0;
            %             break
        else
            updatelam = secder\der;
        end
        lam = lam - 0.6*updatelam;

        %% Update parameters
        % updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
        updateparams = Cth_y*JCey;
        oldparams = vecparams; %?
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

        
        for i = 1:lenreadpone
            ybml(:,i) = C*Utrans(:,i) + D; 
        end

        diff(n1) = norm(data - ybml); 
        %% Search for best parameters \convergence
        halv=0;
        while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)            %LOOP IS TO ENSURE NORM IS DECREASING. PARAMETER UPDATE IS SCALED FOR STEPS IN WHICH CURRENT DEVIATION IS GREATER THAN PREVIOUS
        % disp('halve')
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

            for i = 1:lenreadpone
                ybml(:,i) = C*Utrans(:,i) + D; 
            end
            diff(n1) = norm(data - ybml);
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
            lmat = zeros(1,szlam);
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
    
    %% 
    y = ybml;
    yddr = data;
    % t=1:size(yddr,2);
    % subplot(N,1,1), plot(t,yddr(1,:),'k+',t,y(1,:),'bo')

%% Objective function calculation
J = hwhole;
clear hwhole
rJ=rank(J);

if rJ>=Lenvec
    epss = yddr-y;
    epssTepss = epss.*epss;
    E_R_E = (1/sgma)*sum(epssTepss,"all");                
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                     

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
    Objfun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreadpone-Lenvec-1) + 2*C_ve
else
    Objfun=inf                                                  % This is for a situation we cannot obtain full rank Jacobian for selected sunset
end

Cmat = reshape(C',1,numel(C));
Dmat = D';
lmat = lam';

end
