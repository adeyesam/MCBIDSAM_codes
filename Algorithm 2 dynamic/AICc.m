function [Objfun, Amat, C1mat, C2mat, Dmat, lmat] = AICc(indx,iUY,i_UiUj,UYindex, Ncons, Ain, C1in, C2in, Din, Lin)
% data = datawmcon(:,1:120); Utrans = Utrans(:,1:120); indx=zeros(1,Ntrans); indx([1:10])=1;
% N=4; Ncons = 2; Cin = zeros(4,Ntrans); Din =zeros(1,4); Lin=zeros(1,N^2 +(N-N^2)/2);

load('data4est.mat','data_est','u_est','utrans_est','idOb','nut','mnd','mxd')
data = data_est;
u = u_est;
Utrans = utrans_est;
clear data_est u_est utrans_est

f_indx = find(indx);

%% Estimation and model predictions
sgma = 0.06;   

    % nci = [0 5 5 5 10];               % No. of C in each input variable
    % nhi = [0 6 6 6 12];               % Net No. of H in each species i on LHS 
    % ncio = [0 5 5 5 10];               % No. of C in each species i
    % nhio = [6 6 6 12];              % Net No. of H in each species i on RHS


% [C_act,XQ,Jcbn] = SecOdrInpTransDCM(f_indx,data,u,Utrans,t);
tic
    %% define the parameters
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
    if sum(indx)>0                 % This helps to estimate model with only linear and bias
        for stc = 1:N
            setC2 = [setC2 f_indx+(stc-1)*Ntrans];
        end
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

    count1 = 0;
    lenreading = size(data,2)-1;    %NUMBER OF DATA POINTS/SAMPLES EXCLUDE INTIAL
    lenreadpone = lenreading+1;     
    Lenvec = numsetA + numsetC1 + numsetC2 + numsetD;               %NUMBER OF PARAMETERS
    sz = lenreading*N;              %NUMBER OF DATA READINGS
    % sz1 = lenreadpone*Nwmcon;
    szvec = sz+Lenvec;              %DIMENSION OF LENVEC
    vecparams = zeros(Lenvec,1); 
    szlam = N2 + (N-N2)/2;          % Number of covariance
    meanprior = vecparams;
    
    %% define C and D and obtain vecparams for reduced system
    nulc = isequal(C2in,zeros(1,N*Ntrans));
    if ~nulc
        A = reshape(Ain,[N,N])';
        C1 = reshape(C1in,[lengthu,N])';
        C2 = reshape(C2in,[Ntrans,N])';
        D = Din';
        lam = Lin';
    else
        A = zeros(N);
        C1 = zeros(N,lengthu);
        C2 = zeros(N,Ntrans);
        D = zeros (N,1);
        lam = rand(szlam,1);
    end
    % C = zeros(N,Ntrans);                          

    for i = 1:numsetA
        ar_a(i) = ceil(setA(i)/N);
        ac_a(i) = N*(1 - ceil(mod(setA(i),N)/N)) + mod(setA(i),N);
        if nulc
            A(ar_a(i),ac_a(i)) = rand;
        end
        vecparams(i) = A(ar_a(i),ac_a(i));
    end
    for i = 1:numsetC1
        ar_c1(i) = ceil((setC1(i))/lengthu);
        ac_c1(i) = setC1(i) - lengthu*fix((setC1(i))/lengthu); 
        if ac_c1(i) == 0
            ac_c1(i) = lengthu;
        end
        if nulc
            C1(ar_c1(i),ac_c1(i)) = rand;
        end
        vecparams(i+numsetA) = C1(ar_c1(i),ac_c1(i));
    end
    
    for i = 1:numsetC2
        ar_c2(i) = ceil((setC2(i))/Ntrans);
        ac_c2(i) = setC2(i) - Ntrans*fix((setC2(i))/Ntrans); 
        if ac_c2(i) == 0
            ac_c2(i) = Ntrans;
        end
        if nulc
            C2(ar_c2(i),ac_c2(i)) = rand;
        end
        vecparams(i++numsetA+numsetC1) = C2(ar_c2(i),ac_c2(i));
    end
    
    for i = 1:numsetD
        if nulc
            D(ar_d(i),ac_d(i)) = rand;
        end
        vecparams(i+numsetA+numsetC1+numsetC2) = D(ar_d(i),ac_d(i));
    end
    
    %% Monitors        
    Trend = zeros(Lenvec,30); %
    firstorder = zeros(10,1); % new line to monitor first order optimality.......
    secderr = 1;
    
    %% Initialize xq
    % data = [ones(1,1+ lenreading);data];
    y = data(1:N,:);
    for i = 2:lenreadpone
        y(:,i) = A*y(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D; 
    end
    
    %% M-step variables
%     lam = rand (szlam,1);       % PARAMETERS FOR ESTIMATED SIGNAL COVARIANCE
    der = lam;                  % DERIVATIVE OF F
    secder = zeros(szlam);    %secder = zeros(sz,sz);      %SECOND DERIVATIVE OF F
    cth = 20;
    % Cthprior = 20*eye(Lenvec); % larger premultipliers imply less confidence in priors  %COVARIANCE OF PARAMETERS TO BE ESTIMATED
    % Ceph = zeros(sz);
    
    %% Initialize variables for J matrix
    hwhole = zeros(sz,Lenvec); %h_whole=hwhole;

    % h = zeros(N,lenreadpone);
    % h1 = h; h2 = h;
    DF = hwhole(:,1);
    % D_F = DF;
    
    % span = 2;               %SPAN OF VARRYING PARAMETER TO BE SPINED (HOW MANY STEPS FORWARD, AND THEN BACKWARD)
    % ind = -span:span;       %INDICES TO BE USED TO MULTIPLY PARAMETER FOR THE VARIATION
    % dex = ind +span+1;      %SERIAL NUMBERING OF VARIATIONS
    % tbsplinez = zeros(N,2*span+1); %CREATES MATRIX INTO WHICH OUTPUTS CORRESPONDING TO VARRYING PARAMETERS WILL BE COLLECTED FOR A DATA POINT
    % xrange = zeros(1,2*span+1);     %MATRIX WHERE VARRYING PARAMETER IS STORED
    %% Loop utilities
    %inner
    delta = 0.001;
    tol = 1e-4; % set tolerance
    conv =1;
    n = 0;
    %% Loop Flags
    diff = ones(30,1)*Inf;

    %% Jacobian remains unchanged with model parameters since model is linear in parameters
    
%         for l = 1:numsetC
%             ar = ar_c(l);
%             ac = ac_c(l);
%             count1 = count1 +1;
%             Crac = C(ar,ac);
%             for i = 1:lenreadpone
%                 Cdel=C;
%                 tbspline = tbsplinez;
%                 xrange = Crac + ind*delta;
%                 for g = dex
%                     Cdel(ar,ac) = xrange(g);
%                     tbspline(:,g) = Cdel*Utrans(:,i);
%                 end
%                 h2(:,i) = ppval(fnder(spline(xrange, tbspline)), Crac);
%             end
%             hwhole(:,count1) = reshape(h2(1:end,:)',[sz,1]); h_whole(:,count1) = reshape(h2(1:end,:),[sz,1]);
%         end
%         count1 = 0;
    
   % % Exploiting the linearity of model in parameters to form Jacobian from 
    % (transformed) data to minimise computational burden. (time difference is
    % in the order of 10e3!!!)
    
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
    
    % Checking for estimability: If Jacobian is rank defficient, we do not ...
    % need to carry out estimation...skip the current candidate model
    % Jcbn = hwhole;
    if rank(hwhole)<Lenvec                           
        Amat = zeros(1,N2);
        C1mat = zeros(1,N*lengthu);
        C2mat = zeros(1,N*Ntrans);
        Dmat = zeros(1,N);
        lmat = zeros(1,szlam);
        Objfun = inf;
        return
    end
    tryy = 1;
    
    %% Loops

    %% 'Pre-stage': BML MAP based on data
    while conv > tol && n<80
        n=n+1;
        n1 = n+1;
    
        %% ybar,Jbar
        %         tic
        for i=2:lenreadpone
            yp = A*data(:,i-1) +  C1*u(:,i-1) + C2*Utrans(:,i-1) + D;
            df = data(:,i)- yp; 
            DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(1:end); % em for missing data could come here
    %         D_F((i-2)*N+1:(i-1)*N) = df(1:end);
        end
        
        m_v = meanprior - vecparams;
        ybar = [DF; m_v];
        Jbar = [hwhole; eye(Lenvec)];

        %% M step
        %%  New way as @ 09/07/22
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
    
    %     Cebar = [Ceph zeros(sz,Lenvec); zeros(Lenvec,sz)  Cthprior]; 
    %     P = inv(Cebar) - (Cebar\Jbar)*((Jbar'*(Cebar\Jbar))\Jbar')/Cebar; 
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

        %% new trial 3 : Further Exploiting structure that arise from multiplication with sparse matrices to reduce memory requirement and computational expense
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
                    pv = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading, (rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                    pvpv = pv.*pv';
                    clear pv
                    secder(sdr,sdc) = -0.5 * sum(pvpv, 'all');
                elseif rwi(1)==cli(1)               % If vi & vj are in same row
                    % pvi calculates p * v that corresponds to present ROW OF SECDER
                    pvi = P((rwi(2)-1)*lenreading+1:rwi(2)*lenreading, (rwi(1)-1)*lenreading+1:rwi(1)*lenreading);
                    % pvj calculates p * v that corresponds to present COLUMN OF SECDER
                    pvj = P((cli(2)-1)*lenreading+1:cli(2)*lenreading, (cli(1)-1)*lenreading+1:cli(1)*lenreading);
                    pvipvj = pvi.*pvj';
                    clear pvi pvj
                    secder(sdr,sdc) = -0.5*sum(pvipvj, "all");
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
        secderr = rcond(secder);        %CHECK CONDITIONALITY OF SECOND DERIVATIVE BEFORE INVERSION
        if secderr <1e-13               %IF POORLY CONDITIONED, DON'T UPDATE LAMDA IN PRESENT ITERATION
            updatelam = 0;
            %             break
        else
            updatelam = secder\der;
        end
        conv2 = rms(updatelam);
        lam = lam - 0.6*updatelam;
        Trend_updatelam(:,n) = updatelam;

        %% stability monitors
        
        %% Update parameters
        % updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
        updateparams = Cth_y*JCey;
        oldparams = vecparams; %?
        vecparams = oldparams + updateparams;

        %% Update C,D
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


        diff(n1) = norm(data - ybml); %OBTAINS NORM OF DIFFERENCE BETWEEN ORIGINAL DATA AND WHAT IS GENERATED WITH PARAMETERS OBTAINED
        %% Search for best parameters \convergence

        halv=0;
        while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)            %LOOP IS TO ENSURE NORM IS DECREASING. PARAMETER UPDATE IS SCALED FOR STEPS IN WHICH CURRENT DEVIATION IS GREATER THAN PREVIOUS
    %         disp('halve')
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
            ybml(:,1:2) = data(:,1:2); 
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
            diff(n1) = norm(data - ybml);
        end

        if diff(n1)>diff(n)  && tryy<3               % Trying another starting point and see how it behaves   
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
            % Fmat = zeros(1,N);
            lmat = zeros(1,szlam);
            Objfun = inf;
            return
        end
    
        % Cdel = C;
    
        if n>0
            Trend(:,n) = vecparams;
            firstorder(n) = norm(hwhole,inf);   %TO KEEP TRACK OF FIRST ORDER OPTIMALITY AS ITERATION PROGRESSES
            %% Check conditions for proper exit
            conv = rms(updateparams);
        end
        
    end

    yddr_trend = zeros(sz,10);
    conv_outer = 1;
    nddr = 0;
    
    %% Main Estimation Loops
%     while conv_outer>1e-6 && nddr<5
%         %% Stage 1: DDR based on data and yddr(updated in every iteration)
%         nddr = nddr+1
%         obj = @(x)ddr_obj(x, data, ybml, N, lenreadpone);
%         % [obj, grad, hess, Hstr] = ddr_obj(x, data, ybml, N, lenreadpone);
%         grad = @(x)ddr_objgrad(x, data, ybml, N, lenreadpone);
%         hess = 4*eye(sz);
%         Hstr = @() sparse(~hess==0);
%         
%         nlcon = @(x)ddr_nlcon(x, AtomIN, N, Ncons, lenreadpone, mnd, mxd);
%         % [nlcon, nljac, nljacstr, nlrhs, nle] = ddr_nlcon(x, AtomIN, N, lenreadpone);
%     
%         Ac = zeros(lenreadpone,sz);                 % for C bal
%         Ah = Ac;                                    % for H bal
%         
%         for i = 1:lenreadpone
%             Ac(i,(i-1)*N+1:i*N) = nci;
%             Ah(i,(i-1)*N+1:i*N) = nhio;
%         end
%         AA = [Ac;Ah];
%         clear Ac Ah
%         nljac = @(x) sparse(AA);
%         nljacstr = @() sparse(~AA==0);
%         nlrhs = reshape(AtomIN',Ncons*lenreadpone,1);
%         nle = zeros(Ncons*lenreadpone,1);
%     
%         % clear AA
%         lb = .0*reshape(data,sz,1);
%         ub = 1e5*reshape(data,sz,1);
%         x0 = reshape(data,sz,1);
%     
%         % size(obj)
%         % size(grad)
%         % size(hess)
%         Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'ndec', sz); %,'nljac',nljac,'ndec',sz);%,'nljacstr',nljacstr);
%         % Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr,'bounds',lb,ub);
%         [yddr,fval,exitflag,info] = solve(Opt,x0);
% 
%         %yddr = fmincon(obj, x0, [], [], AA, nlrhs, lb, ub);
% 
%         yddr_trend(:,nddr) = yddr;
%         yddr = reshape(yddr,N,lenreadpone);
%         % 
%         % if nddr==1
%         %      conv_outer =1;
%         % else
%         %      conv_outer = rms(yddr_trend(:,nddr) - yddr_trend(:,nddr-1))
%         % end
%     
%         %% Stage 2: Bayesian MAP based on reconciled data, yddr
%         n = 0;
%         conv2 =1;
%         while conv2 > tol && n<80
%             n=n+1;
%             n1 = n+1;
%         
%             %% ybar,Jbar
%             %         tic
%             for i=1:lenreadpone
%                 yp = C*Utrans(:,i) + D;
%                 df = yddr(:,i)- yp; 
%                 DF(i:lenreadpone:(N-1)*lenreadpone +i) = df(1:end); % em for missing data could come here
%         %         D_F((i-2)*N+1:(i-1)*N) = df(1:end);
%             end
%             
%             m_v = meanprior - vecparams;
%             ybar = [DF; m_v];
%             Jbar = [hwhole; eye(Lenvec)];
%     
%             %% M step
%             %%  New way as @ 09/07/22
%             ceph_u = zeros(N); Ce_invu = zeros(sz);
%             % Construct unit covariance matrix
%             lamn = 0;
%             for ii=1:N
%                 jj=ii:N;
%                 for kk=jj
%                     lamn=lamn+1;
%                     ceph_u(ii,kk) = lam(lamn);
%                 end
%             end
%             ceph = ceph_u + triu(ceph_u,1)';
%     
%             ceph_inv = ceph\eye(N);
%             % Expand inverse of unit covariance matrix for lenreading points
%             lamn = 0;
%             for ii=1:N
%                 jj=ii:N;
%                 for kk=jj
%                     lamn=lamn+1;
%                     Ce_invu((ii-1)*lenreadpone+1:ii*lenreadpone,(kk-1)*lenreadpone+1:kk*lenreadpone) = ceph_inv(ii,kk)*eye(lenreadpone);
%                 end
%             end
%             Ce_inv = sparse(Ce_invu + triu(Ce_invu,1)');
%     
%             clear Ce_invu
%         
%         %     Cebar = [Ceph zeros(sz,Lenvec); zeros(Lenvec,sz)  Cthprior]; 
%         %     P = inv(Cebar) - (Cebar\Jbar)*((Jbar'*(Cebar\Jbar))\Jbar')/Cebar; 
%             Cth_inv = (1/cth)*eye(Lenvec);
%             Ce_y1 = Ce_inv*DF;
%             Ce_y2 = Cth_inv*m_v;
%             Ce_y = [Ce_y1; Ce_y2];
%             JCey = Jbar'*Ce_y;
%             clear Ce_y
%     
%             Ce_J = [Ce_inv*hwhole; Cth_inv];
%             JCeJ = Jbar'*Ce_J;
%             Cth_y = JCeJ\eye(Lenvec);
%             clear JCeJ
%     
%             Cebar_inv = sparse([Ce_inv zeros(sz,Lenvec); zeros(Lenvec,sz)  Cth_inv]);
%             JCe = Jbar'*Cebar_inv;
%             CeJCth = Ce_J*Cth_y;
%             CeJCthJCe = CeJCth*JCe;         
%             clear CeJCth Ce_J Ce_inv
%     
%             P = Cebar_inv - CeJCthJCe;
%             clear Cebar_inv
%     
%             %% new trial 3 : Further Exploiting structure that arise from multiplication with sparse matrices to reduce memory requirement and computational expense
%             lamn=0;
%             for ii=1:N
%                 jj=ii:N;
%                 for kk=jj
%                     vpy=zeros(szvec,1);
%                     vpy((ii-1)*lenreadpone+1:ii*lenreadpone) = P((kk-1)*lenreadpone+1:kk*lenreadpone,:)*ybar;
%                     % Trace of pv is the trace of nonzero square matrix that falls on the diagonal
%                     tr_pv = sum(diag(P((kk-1)*lenreadpone+1:kk*lenreadpone,(ii-1)*lenreadpone+1:ii*lenreadpone)));
%                     
%                     lamn=lamn+1;
%                     der(lamn) = -0.5*tr_pv + 0.5*ybar'*P'*vpy;
%             
%                     rc(lamn,:) = [ii kk];                   % row-col indices for lam          
%             
%                 end
%             end
%     
%             for sdr=1:szlam                             % Secder row
%                 sdci=sdr:szlam;                         % Secder column indices
%                 for sdc=sdci
%                     rwi = rc(sdr,:);                    % row-col indices of vi
%                     cli = rc(sdc,:);                    % row-col indices of vj
%     
%                     % For pv/pvi/pvj, only selecting elements of P that will
%                     % eventually matter for obtaining trace(PViPVj) as P is
%                     % multiplied with each v/vi/vj
%     
%                     if sdr==sdc                         % pvipvj for principal diagonal
%                         pv = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone,(rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
%                         pvpv = pv.*pv';
%                         clear pv
%                         secder(sdr,sdc) = -0.5*sum(pvpv,'all');
%                     elseif rwi(1)==cli(1)               % If vi & vj are in same row
%                         % pvi calculates p * v that corresponds to present ROW OF SECDER
%                         pvi = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone,(rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
%                         % pvj calculates p * v that corresponds to present COLUMN OF SECDER
%                         pvj = P((cli(2)-1)*lenreadpone+1:cli(2)*lenreadpone,(cli(1)-1)*lenreadpone+1:cli(1)*lenreadpone);
%                         pvipvj = pvi.*pvj';
%                         clear pvi pvj
%                         secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
%                         clear pvipvj
%                         secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
%                     else
%                         % pvi calculates p * v that corresponds to present ROW OF SECDER
%                         pvi = P((cli(2)-1)*lenreadpone+1:cli(2)*lenreadpone,(rwi(1)-1)*lenreadpone+1:rwi(1)*lenreadpone);
%                         % pvj calculates p * v that corresponds to present COLUMN OF SECDER
%                         pvj = P((rwi(2)-1)*lenreadpone+1:rwi(2)*lenreadpone,(cli(1)-1)*lenreadpone+1:cli(1)*lenreadpone);
%                         pvipvj = pvi.*pvj';
%                         clear pvi pvj
%                         secder(sdr,sdc) = -0.5*sum(pvipvj,"all");
%                         clear pvipvj
%                         secder(sdc,sdr) = secder(sdr,sdc);          % trace is unchanged by transposing a matrix
%                         
%                     end
%                 end
%             end
%     
%         %%
%             secderr = rcond(secder);        %CHECK CONDITIONALITY OF SECOND DERIVATIVE BEFORE INVERSION
%             if secderr <1e-13               %IF POORLY CONDITIONED, DON'T UPDATE LAMDA IN PRESENT ITERATION
%                 updatelam = 0;
%                 %             break
%             else
%                 updatelam = secder\der;
%             end
%             conv2 = rms(updatelam);
%             lam = lam - 0.6*updatelam;
%             Trend_updatelam(:,n) = updatelam;
%     
%             %% stability monitors
%             
%             %% Update parameters
%             % updateparams = (Jbar'*(Cebar\Jbar))\(Jbar'*(Cebar\ybar));
%             updateparams = Cth_y*JCey;
%             oldparams = vecparams; %?
%             vecparams = oldparams + updateparams;
%     
%             %% Update C,D
%             for j = 1:numsetC
%                 ar = ar_c(j);
%                 ac = ac_c(j);
%                 C(ar,ac) = vecparams(j);
%             end
%             
%             for k = 1:numsetD
%                 dr = ar_d(k);
%                 dc = ac_d(k);
%                 D(dr,dc) = vecparams(k+numsetC);
%             end
%     
%             
%             for i = 1:lenreadpone
%                 y(:,i) = C*Utrans(:,i) + D; 
%             end
%     
%             diff(n1) = norm(yddr - y); %OBTAINS NORM OF DIFFERENCE BETWEEN ORIGINAL DATA AND WHAT IS GENERATED WITH PARAMETERS OBTAINED
%             %% Search for best parameters \convergence
%     
%             halv=0;
%             while (diff(n1)>diff(n))  && (halv<100) &&  (tryy<3)            %LOOP IS TO ENSURE NORM IS DECREASING. PARAMETER UPDATE IS SCALED FOR STEPS IN WHICH CURRENT DEVIATION IS GREATER THAN PREVIOUS
%         %         disp('halve')
%                 halv = halv + 1;
%                 updateparams = 0.6*updateparams;
%                 vecparams = oldparams + updateparams;
%                 
%                 for j = 1:numsetC
%                     ar = ar_c(j);
%                     ac = ac_c(j);
%                     C(ar,ac) = vecparams(j);
%                 end
%         
%                 for k = 1:numsetD
%                     dr = ar_d(k);
%                     dc = ac_d(k);
%                     D(dr,dc) = vecparams(k+numsetC);
%                 end
%     
%                 for i = 1:lenreadpone
%                     y(:,i) = C*Utrans(:,i) + D; 
%                 end
%                 diff(n1) = norm(yddr - y);
%             end
%     
%             if diff(n1)>diff(n)  && tryy<3               % Trying another starting point and see how it behaves   
%                 tryy = tryy +1
%                 vecparams = rand(Lenvec,1)*tryy^2;
%                 for j = 1:numsetC
%                     ar = ar_c(j);
%                     ac = ac_c(j);
%                     C(ar,ac) = vecparams(j);
%                 end
%                 for k = 1:numsetD
%                     dr = ar_d(k);
%                     dc = ac_d(k);
%                     D(dr,dc) = vecparams(k+numsetC);
%                 end
%                 lam=rand(szlam,1); n=0;
%                 diff = ones(30,1)*Inf; conv2=1;
%                 
%             elseif diff(n1)>diff(n)
%                 Cmat = zeros(1,N*Ntrans);
%                 Dmat = zeros(1,N);
%                 % Fmat = zeros(1,N);
%                 lmat = zeros(1,szlam);
%                 Objfun = inf;
%                 return
%             end
%         
%             % Cdel = C;
%         
%             if n>0
%                 Trend(:,n) = vecparams;
%                 firstorder(n) = norm(hwhole,inf);   %TO KEEP TRACK OF FIRST ORDER OPTIMALITY AS ITERATION PROGRESSES
%                 %% Check conditions for proper exit
%                 conv2 = rms(updateparams);
%             end
%             
%         end
%         ybml = y;
% 
%         ybml_trend(:,nddr) = reshape(ybml,1,sz);
%         if nddr==1
%              conv_outer =1;
%         else
%              conv_outer = rms(ybml_trend(:,nddr) - ybml_trend(:,nddr-1))
%         end
%      end
    
y = ybml;
yddr = data;
    t=1:size(data,2);
    subplot(N,1,1), plot(t,data(1,:),'k+',t,y(1,:),'bo')

%% Jacobian and Likelihood fcn
J = hwhole;%(1:Lenvec,:);
clear hwhole
size(J)
rJ=rank(J);
% inv_r = (sigma*eye(sz))\eye(sz);

if rJ>=Lenvec
    % E_R_E = 0;
    % epss = zeros(N,lenreading);
    % parfor i = 1 : lenreading
    %     epss(:,i) = (data(:,i) - XQ(:,i));
    %     E_R_E = E_R_E + epss(:,i)'*inv_R*epss(:,i);
    % end
    epss = data-y;
    epssTepss = epss.*epss;
    E_R_E = (1/sgma)*sum(epssTepss,"all");                  % R^-1 = 1/sigma * eye(N) so we cann simply bring sigma out
    clear epss epssTepss data u Utrans

    det_R = sgma^N;                                     % Determinant of diagonal matrix=product of the diagonal  elements

    L_theta = -0.5*N*lenreading*log(2*pi) - 0.5*lenreading*log(det_R) - 0.5*E_R_E;
   
    % JRJ=J'*inv_r*J; rJRJ=rank(JRJ);
    % Co_theta = inv(J'*inv_r*J);

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
    Objfun=inf                                                  % This is for a situation we cannot obtain full rank Jacobian for selected sunset
end

Amat = reshape(A',1,N2);
C1mat = reshape(C1',1,numel(C1));
C2mat = reshape(C2',1,numel(C2));
Dmat = D';
lmat = lam';

fprintf('tryy                       %d \n',tryy)
fprintf('n                          %d \n',n)
fprintf('halv                       %d \n',halv)
% fprintf('UY bases                   %d \n',UYpresent)
fprintf('Cve                        %d \n',C_ve)
fprintf('Loglikelihood              %d \n',L_theta)
fprintf('Objfun                     %d \n',Objfun)
fprintf('estimation time            %2f \n',toc)

end
