function obj_p = obj3(x, data, Utrans, u, N, sz, lenreadpone, numsetC1, numsetC2, numsetD, ar_c1, ac_c1, ar_c2, ac_c2, ar_d, ac_d)
    
    % Distributing paramaters from vectors into matrices
    Ntrans = size(Utrans, 1);
    lengthu = size(u,1);

    C1 = zeros(N,lengthu);
    C2 = zeros(N, Ntrans);
    D = zeros(N,1);

    for j = 1:numsetC1
        ar = ar_c1(j);
        ac = ac_c1(j);
        C1(ar,ac) = x(j);
    end

    for j = 1:numsetC2
        ar = ar_c2(j);
        ac = ac_c2(j);
        C2(ar,ac) = x(j+numsetC1);
    end
    
    for k = 1:numsetD
        dr = ar_d(k);
        dc = ac_d(k);
        D(dr,dc) = x(k+numsetC1+numsetC2);
    end
    %Calculateing model predictions
    ybml = zeros(N,lenreadpone);
    DF = zeros(sz,1);
    for i = 1:lenreadpone
        ybml(:,i) = C1*u(:,i) + C2*Utrans(:,i) + D; 
        df = data(:,i)- ybml(:,i); 
        DF(i:lenreadpone:(N-1)*lenreadpone +i) = df(1:end); % em for missing data could come here
        sse(i) = sum(df.^2, 'all');
    end

    obj_p = sum(sse,'all');

end
