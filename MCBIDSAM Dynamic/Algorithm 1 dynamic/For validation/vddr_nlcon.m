function nlcon =  vddr_nlcon(x, us, nci, ncio, N, F, mnd, mxd)

xx = reshape(x,N,1);
for i=1:N 
    xxs(i,:) = mnd(i) + 0.5*(xx(i,:)-1)*(mxd(i)-mnd(i));
end

c_in = nci*us;
c_out = ncio*xxs;
nlcon = c_in - c_out;


return