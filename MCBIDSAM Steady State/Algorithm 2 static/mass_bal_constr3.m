function nlcon = mass_bal_constr3(x, indx, N, lengthu, Ncons, nci, ncio, nhi, nhio, mnd, mxd, umin, umax)
    if Ncons>0
        f_indx = find(indx);                    % Positions of current basis fcns
        ntrans = numel(f_indx);                 % No. of basis fcns in current model
        numsetC1 = N*lengthu;
        numsetC2 = ntrans*N;
        C1 = reshape(x(1:numsetC1), lengthu, N)';
        C2 = reshape(x(numsetC1+1:numsetC1+numsetC2), ntrans, N)';
        xvec1 = reshape(C1, numsetC1,1);
        xvec2 = reshape(C2, numsetC2,1);           % Rearranging vecparams such ...
                                                % that it now counts columnwise of C instead ...
                                                % of rowise as it is allthrough the algorithms
                                                % and also in 'obj_posterior3'

        AA1c = zeros(lengthu,N*lengthu);
        AA1h = AA1c;
        AA2c = zeros(ntrans,N*ntrans);
        AA2h = AA2c;

        dy = mxd-mnd;
        du = umax-umin; 

        for g=1:lengthu
            nci_(g) = 0.5*(du(g))*nci(g);
            nhi_(g) = 0.5*(du(g))*nhi(g);
        end

        for g=1:N
            ncio_(g) = 0.5*(dy(g))*ncio(g);
            nhio_(g) = 0.5*(dy(g))*nhio(g);
        end
            
        for k=1:lengthu
            AA1c(k, (k-1)*N+1:k*N) = ncio_;
            AA1h(k, (k-1)*N+1:k*N) = nhio_;
        end
        
        for k=1:ntrans
            AA2c(k, (k-1)*N+1:k*N) = ncio_;
            AA2h(k, (k-1)*N+1:k*N) = nhio_;
        end
        AA1 = [AA1c; AA1h];
        AA2 = [AA2c; AA2h];
        clear AA1c AA1h AA2c AA2h

        % Check which current basis fcn is which u
        % nlrhs1 = zeros(ntrans*Ncons,1);
        % fu = [];
        % for j=1:ntrans
        %     for k = 1:numel(Uindx)
        %         if f_indx(j) == Uindx(k)
        %             fu = [fu [j;k]];               % [basis fcn; what u it is]
        %         end
        %     end
        % end

        % Only basis functions that are linear forms of u's rep.ing concentration have RHS in atom balance
        % for jj=1:size(fu,2)
        %     fu1 = fu(1,jj);
        %     fu2 = fu(2,jj);
        % 
        %     nlrhs1(fu1) = nci_(fu2);
        %     % AA(fu1,(fu1-1)*N+1:fu1*N) = nci;
        % 
        %     nlrhs1(fu1+ntrans) = nhi_(fu2);
        %     % AA(fu1+ntrans,(fu1-1)*N+1:fu1*N) = nhi;
        % end
        % xvec;
        % x;
        % AA*xvec;

        nlrhs1 = zeros(Ncons*(lengthu+ntrans),1);
        nlrhs1(1:2*lengthu) = [nci_ nhi_];         % for linear input variables

        nlcon1 = [AA1*xvec1; AA2*xvec2] - nlrhs1 ;         % Coefficients of basis functions individually set to zero
        
        nlcon2a = 0;                         % constant term from C-bal set to zero
        nlcon3a = 0;                         % constant term from H-bal set to zero
        for l=1:N
            nlcon2a = nlcon2a + ncio(l)*(mnd(l) + 0.5*dy(l)*x(numsetC1+numsetC2+l) - 0.5*dy(l));
            nlcon3a = nlcon3a + nhio(l)*(mnd(l) + 0.5*dy(l)*x(numsetC1+numsetC2+l) - 0.5*dy(l));
        end

        nlcon2b = 0;
        nlcon3b = 0;
        for j=1:lengthu
            nlcon2b = nlcon2b + nci(j)*(umin(j) - 0.5*du(j));
            nlcon3b = nlcon3b + nhi(j)*(umin(j) - 0.5*du(j));
        end

        nlcon2 = nlcon2a - nlcon2b;
        nlcon3 = nlcon3a - nlcon3b;

        nlcon = [nlcon1; nlcon2; nlcon3];
        % nlcon = nlcon.^2;
    else
        nlcon = [];
    end
end
