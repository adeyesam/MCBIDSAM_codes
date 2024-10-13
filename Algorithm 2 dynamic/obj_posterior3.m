function obj_p = obj_posterior3(x, data, Utrans, u, Lenvec, N, sz, lenreading, numsetA, numsetC1, numsetC2, numsetD, ar_a, ac_a, ar_c1, ac_c1, ar_c2, ac_c2, ar_d, ac_d, lam, cth, meanprior)
    
    Ntrans = size(Utrans, 1);
    lengthu = size(u,1);
    lenreadpone = lenreading+1;

    A = zeros(N);
    C1 = zeros(N,lengthu);
    C2 = zeros(N, Ntrans);
    D = zeros(N,1);

    for j = 1:numsetA
        ar = ar_a(j);
        ac = ac_a(j);
        A(ar,ac) = x(j);
    end
    
    for j = 1:numsetC1
        ar = ar_c1(j);
        ac = ac_c1(j);
        C1(ar,ac) = x(j+numsetA);
    end

    for j = 1:numsetC2
        ar = ar_c2(j);
        ac = ac_c2(j);
        C2(ar,ac) = x(j+numsetA+numsetC1);
    end
    
    for k = 1:numsetD
        dr = ar_d(k);
        dc = ac_d(k);
        D(dr,dc) = x(k+numsetA+numsetC1+numsetC2);
    end
    %Calculateing model predictions
    ybml = zeros(N,lenreading);
    DF = zeros(sz,1);
    for i = 2:lenreadpone
        ybml(:,i) =  A*ybml(:,i-1) + C1*u(:,i-1) + C2*Utrans(:,i-1) + D;
        df = data(:,i)- ybml(:,i); 
        DF(i-1:lenreading:(N-1)*lenreading +i-1) = df(1:end); 
        sse(i-1) = sum(df.^2, 'all');
    end
    obj_p = sum(sse,'all');

end
