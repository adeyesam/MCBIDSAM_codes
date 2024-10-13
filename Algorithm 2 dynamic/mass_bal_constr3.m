function nlcon = mass_bal_constr3(x, indx, N, lengthu, Ncons, Lenvec, nci, ncio, nhi, nhio, mnd, mxd, umin, umax)
    if Ncons>0
        f_indx = find(indx);                    % Positions of current basis fcns
        ntrans = numel(f_indx);                 % No. of basis fcns in current model
        numsetA = N^2;
        numsetC1 = N*lengthu;
        numsetC2 = ntrans*N;
        numsetD = N;
        A = reshape(x(1:numsetA),N,N)';
        C1 = reshape(x(numsetA+1:numsetA+numsetC1), lengthu, N)';
        C2 = reshape(x(numsetA+numsetC1+1:numsetA+numsetC1+numsetC2), ntrans, N)';
        D = x(numsetA+numsetC1+numsetC2+1:end);
        xvec1 = reshape(A,numsetA,1);
        xvec2 = reshape(C1, numsetC1,1);
        xvec3 = reshape(C2, numsetC2,1);        % Rearranging vecparams such ...
                                                % that it now counts columnwise of A/C1/C2 instead ...
                                                % of rowise as it is allthrough the algorithms
                                                % and also in 'obj_posterior3'

        AA1c = zeros(N,N^2);
        AA1h = AA1c;
        AA2c = zeros(lengthu,N*lengthu);
        AA2h = AA2c;
        AA3c = zeros(ntrans,N*ntrans);
        AA3h = AA3c;

        dy = mxd-mnd;
        du = umax-umin; 

        for g=1:lengthu
            nci_(g) = 0.5*(du(g))*nci(g);
            nhi_(g) = 0.5*(du(g))*nhi(g);
        end

        for g=1:N
            ncio_(g) = 0.5*(dy(g))*ncio(g)/(1-A(g,g));
            nhio_(g) = 0.5*(dy(g))*nhio(g)/(1-A(g,g));
        end
            
        for k=1:N
            ncio_p = ncio_;
            ncio_p(k) = 0;
            nhio_p = nhio_;
            nhio_p(k) = 0;
            AA1c(k, (k-1)*N+1:k*N) = ncio_p;
            AA1h(k, (k-1)*N+1:k*N) = nhio_p;
        end

        for k=1:lengthu
            AA2c(k, (k-1)*N+1:k*N) = ncio_;
            AA2h(k, (k-1)*N+1:k*N) = nhio_;
        end
        
        for k=1:ntrans
            AA3c(k, (k-1)*N+1:k*N) = ncio_;
            AA3h(k, (k-1)*N+1:k*N) = nhio_;
        end

        nlrhs1 = zeros(Ncons*(N+lengthu+ntrans),1);
        nlrhs1(N+1:N+lengthu) = nci_;                            % for linear input variables
        nlrhs1(2*N+lengthu+ntrans+1:2*N+2*lengthu+ntrans) = nhi_;         % for linear input variables

        nlcon1c = [AA1c*xvec1; AA2c*xvec2; AA3c*xvec3] ;         % Coefficients of basis functions individually set to zero for C bal
        nlcon1h = [ AA1h*xvec1; AA2h*xvec2; AA3h*xvec3];         % Coefficients of basis functions individually set to zero for H bal
        nlcon1 = [nlcon1c; nlcon1h] - nlrhs1;
        clear AA1c AA1h AA2c AA2h AA3c AA3h xvec1 xvec2 xvec3

        nlcon2a = 0;                         % constant term from C-bal set to zero
        nlcon3a = 0;                         % constant term from H-bal set to zero
        for l=1:N
            nlcon2a = nlcon2a + ncio(l)*(mnd(l) + 0.5*dy(l)*D(l)/(1-A(l,l)) - 0.5*dy(l));
            nlcon3a = nlcon3a + nhio(l)*(mnd(l) + 0.5*dy(l)*D(l)/(1-A(l,l)) - 0.5*dy(l));
        end

        nlcon2b = 0;
        nlcon3b = 0;
        for j=1:lengthu
            nlcon2b = nlcon2b + nci(j)*(umin(j) - 0.5*du(j));
            nlcon3b = nlcon3b + nhi(j)*(umin(j) - 0.5*du(j));
        end

        nlcon2 = nlcon2a - nlcon2b;              % Connstrains matrix D elts based on C bal
        nlcon3 = nlcon3a - nlcon3b;              % Connstrains matrix D elts based on H bal

        nlcon = [nlcon1; nlcon2; nlcon3];
    else
        nlcon = [];
    end
end
