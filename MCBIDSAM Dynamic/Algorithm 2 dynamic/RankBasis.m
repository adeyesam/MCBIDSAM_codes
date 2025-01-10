function [Objfun, Amat, C1mat, C2mat, Dmat] = RankBasis(indx,iUY,i_UiUj,UYindex)
    
    load('data4rank.mat','data_rank','u_rank','utrans_rank')
    data = data_rank;
    u = u_rank;
    Utrans = utrans_rank;
    clear data_rank u_rank utrans_rank
    f_indx = find(indx);
    
    %% Estimation and model predictions
    sgma = 0.06;   
tic
    %% define the parameters
    N = size(data,1);
    N2 = N^2;                                              
    Ntrans = size(Utrans,1);                                
    ntrans = numel(f_indx);                                 % No of basis fcns currently considered
    lengthu = size(u,1);

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
    Lenvec = numsetA + numsetC1 + numsetC2 +numsetD;               
    sz = lenreading*N;             
    szvec = sz+Lenvec;              
    vecparams = zeros(Lenvec,1); 
    szlam = N2 + (N-N2)/2;    
    meanprior = vecparams;
    
    %% define C and D and obtain vecparams for reduced system
    A = zeros(N);                          
    for i = 1:numsetA
        ar_a(i) = ceil(setA(i)/N);
        ac_a(i) = N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N);
        A(ar_a(i),ac_a(i)) = rand;
        vecparams(i) = A(ar_a(i),ac_a(i));
    end
    
    C1 = zeros(N,lengthu);   
    C2 = zeros(N,Ntrans);   

    for i = 1:numsetC1
        ar_c1(i) = ceil((setC1(i))/lengthu);
        ac_c1(i) = setC1(i) - lengthu*fix((setC1(i))/lengthu); 
        if ac_c1(i) == 0
            ac_c1(i) = lengthu;
        end
        C1(ar_c1(i),ac_c1(i)) = rand;
        vecparams(i+numsetA) = C1(ar_c1(i),ac_c1(i));
    end
    
    for i = 1:numsetC2
        ar_c2(i) = ceil((setC2(i))/Ntrans);
        ac_c2(i) = setC2(i) - Ntrans*fix((setC2(i))/Ntrans); 
        if ac_c2(i) == 0
            ac_c2(i) = Ntrans;
        end
        C2(ar_c2(i),ac_c2(i)) = rand;
        vecparams(i++numsetA+numsetC1) = C2(ar_c2(i),ac_c2(i));
    end

    for i = 1:numsetD
        D(ar_d(i),ac_d(i)) = rand;
        vecparams(i+numsetA+numsetC1+numsetC2) = D(ar_d(i),ac_d(i));
    end
    
    %% Monitors        
    Trend = zeros(Lenvec,30); 
    firstorder = zeros(10,1); % new line to monitor first order optimality.......
    
    %% Initialize xq
    y = data(1:N,:);
    for i = 2:lenreadpone
        y(:,i) = A*y(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
    end
    

    %% M-step variables
    lam =  rand (szlam,1); 
    der = lam;                  
    secder = zeros(szlam);    
    cth = 20;
    
    %% Initialize variables for J matrix
    hwhole = zeros(sz,Lenvec); 
    DF = hwhole(:,1);

    %% Loop utilities
    tol = 1e-4; 
    conv =1;
    n = 0;
    diff = ones(30,1)*Inf;

    %% Exploiting the linearity of model in parameters to form Jacobian from 
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
    
    % Checking for estimability
    if rank(hwhole)<Lenvec   
        Amat = zeros(1,N2);
        C1mat = zeros(1,N*lengthu);
        C2mat = zeros(1,N*Ntrans);
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
        for i=2:lenreadpone
            yp = A*data(:,i-1) +  C1*u(:,i-1) + C2*Utrans(:,i-1) + D;
            df = data(:,i)- yp; 
            DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(1:end); 
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
                    pvi = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading,(rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                    pvj = P((cli(2)-1)*lenreading+1:cli(2)*lenreading,(cli(1)-1)*lenreading+1:cli(1)*lenreading);
                    pvipvj = pvi.*pvj';
                    clear pvi pvj
                    secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
                    clear pvipvj
                    secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
                else
                    pvi = P((cli(2)-1)*lenreading+1:cli(2)*lenreading,(rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
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
        updateparams = Cth_y*JCey;
        oldparams = vecparams; %?
        vecparams = oldparams + updateparams;
        
        %% Update A, C, D
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
        if ~isempty(UYpresent)
            % disp('UY present')
            y(:,1) = data(:,1); 
            Utransi=zeros(Ntrans,lenreading);
            for ntr=1:lenreading
                [Utransi(:,ntr),~,~] = InputTransformation(y(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                y(:,ntr+1) = A*y(:,ntr) + C1*u(:,ntr) + C2*Utransi(:,ntr) + D;
            end
        else
            for i = 2:lenreadpone
                y(:,i) =  A*y(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
            end
        end

        diff(n1) = norm(data - y); 
        %% Search for best parameters \convergence
        halv=0;
        while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)            
            halv = halv + 1;
            updateparams = 0.6*updateparams;
            vecparams = oldparams + updateparams;
            
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
        if ~isempty(UYpresent)
            % disp('UY present')
            y(:,1) = data(:,1); 
            Utransi=zeros(Ntrans,lenreading);
            for ntr=1:lenreading
                [Utransi(:,ntr),~,~] = InputTransformation(y(:,ntr),u(:,ntr),iUY,i_UiUj,ntr);
                y(:,ntr+1) = A*y(:,ntr) + C1*u(:,ntr) + C2*Utransi(:,ntr) + D;
            end
        else
            for i = 2:lenreadpone
                y(:,i) =  A*y(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
            end
        end
            diff(n1) = norm(data - y);
        end
    
        if diff(n1)>diff(n)  && tryy<3               
            tryy = tryy +1
            vecparams = rand(Lenvec,1)*tryy^2;
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

            lam=rand(szlam,1); n=0;
            diff = ones(30,1)*Inf; conv=1;
            
        elseif diff(n1)>diff(n)
            Amat = zeros(1,N2);
            C1mat = zeros(1,N*lengthu);
            C2mat = zeros(1,N*Ntrans);
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
size(J)
rJ=rank(J);

if rJ>=Lenvec
    epss = data-y;
    epssTepss = epss.*epss;
    E_R_E = (1/sgma)*sum(epssTepss,"all");               
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                     % Determinant of diagonal matrix=product of the diagonal  elements

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
    Objfun = -2*L_theta + 2*Lenvec + (2*Lenvec*(Lenvec+1))/(lenreading-Lenvec-1) + 2*C_ve
else
    Objfun=inf                                                  % This is for a situation we cannot obtain full rank Jacobian for selected sunset
end
Amat = reshape(A',1,numel(A));
C1mat = reshape(C1',1,numel(C1));
C2mat = reshape(C2',1,numel(C2));
Dmat = D';

end
