function nlcon =  ddr_nlcon(x, AtomIN, N, Ncons, lenreadpone, ncio, nhio)

sz = N*lenreadpone;

%% Atom bal constraints
Ac = zeros(lenreadpone,sz);                 % for C bal
Ah = Ac;                                    % for H bal

for i = 1:lenreadpone
    Ac(i,(i-1)*N+1:i*N) = ncio;
    Ah(i,(i-1)*N+1:i*N) = nhio;
end
AA = [Ac;Ah];

nlcon1 = Ac*x(N+1:N+sz);                    % constraints (they are linear in this case)
nlcon2 = Ah*x(N+1:N+sz); 

clear Ac Ah
nlrhs = reshape(AtomIN',Ncons*lenreadpone,1);
nlcon = [nlcon1; nlcon2] - nlrhs;

clear AA nlcon1 nlcon2 

return