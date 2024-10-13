function nlcon =  ddr_nlcon(x, AtomIN, N, Ncons, sel_indx, len, mnd, mxd,nci,nhio)
sz_ss = N*len;
lenreading = size(AtomIN,2)-1;

xx = reshape(x(N+1:end),N,lenreading);
xvec = reshape(xx(:,sel_indx), len*N,1);  % select and vectorise identified steady states

%% Atom bal constraints
Ac = zeros(len,sz_ss);                    % for C bal
Ah = Ac;                                  % for H bal

for i = 1:len
    Ac(i,(i-1)*N+1:i*N) = nci;
    Ah(i,(i-1)*N+1:i*N) = nhio;
end
AA = [Ac;Ah];

nlcon1 = Ac*xvec; % x(N+1:N+sz);          % constraints (they are linear in this case)
nlcon2 = Ah*xvec; % x(N+1:N+sz); 
clear Ac Ah

nlrhs = reshape(AtomIN(:,sel_indx)',Ncons*len,1);
nlcon = [nlcon1; nlcon2] - nlrhs;

clear AA nlcon1 nlcon2 
return