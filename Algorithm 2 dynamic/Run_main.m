%% Load data 
% load('data.mat');
% load('u.mat');
%% 
N = size(data,1);               % NUMBER OF OUTPUT
N2 = N^2;                       % NUMBER OF ELEMENTS IN A
lengthu = size(u,1);            % NUMBER OF INPUTS
lengthu2 = size(u,2);           % Number of input samples
lenreading = size(data,2)-1;
lenreadpone = lenreading+1;

%% Adding element bal to data
nci = [1 0 1 0];                % No. of C in each input variable (F/V Caf Cbf Ccf Cdf)
nhi = [0 2 0 2];                % Net No. of H in each species i on LHS 
ncio = [1 0 1 0];               % No. of C in each input variable (F/V Caf Cbf Ccf Cdf)
nhio = [0 2 0 2];               % Net No. of H in each species i on RHS

cbal = zeros(1,lenreadpone);
hbal = cbal;
for i=1:lenreadpone
    cbal(i) = nci*u(:,i);
    hbal(i) = nhi*u(:,i);
end
% datawmcon = data;
datawmcon = [data;cbal;hbal];    % Data with mass constr
Nwcons = size(datawmcon,1);
Ncons = Nwcons - N;              % No. of mass constr imposed

%% Analyses for steady state
sel_indx_all = Identify_SteadyState(data,u);

%%
umax=zeros(1,lengthu); 
umin=umax;
for i=1:lengthu
    umax(i) = max(u(i,:));
    umin(i) = min(u(i,:));
end
for i=1:Nwcons
    mxd(i)=max(datawmcon(i,:));
    mnd(i)=min(datawmcon(i,:));
end
% Scale to range 1:3 
for i=1:N
    datawmcon(i,:)=1 + 2*(datawmcon(i,:)-mnd(i))/(mxd(i)-mnd(i));
end
for i = 1:size(u,1)
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end

%% Compute input transformations
[Utrans, iUY, i_UiUj] = InputTransformation(data,u,[],[],true);
Ntrans = size(Utrans,1);         

[Uindx, UYindex1] = BasisInterpretation(lengthu,iUY,i_UiUj,[],[],[],false);

%% Selecting specific range(s) of data for ranking and estimation
rank_sel =1:500;  
est_sel =1:800;
fine_sel = est_sel;

data_rank=datawmcon(:,rank_sel);
u_rank=u(:,rank_sel);
utrans_rank=Utrans(:,rank_sel);

data_est=datawmcon(:,est_sel);
u_est=u(:,est_sel);
utrans_est=Utrans(:,est_sel);

data_fine = datawmcon(:,fine_sel);
u_fine = u(:,fine_sel);
utrans_fine = Utrans(:,fine_sel);

%% Analyses of training data for steady state
sel_indx = Identify_SteadyState(data(:,est_sel),u);

%% Ranking basis functions
save('data4rank','data_rank','u_rank','utrans_rank','mnd','mxd')
clear data_rank u_rank utrans_rank

gh=zeros(1,Ntrans,Ntrans);
for i=1:Ntrans
    gh(:,i,i)=1;
end

ind=zeros(1,Ntrans);
parfor i=1:Ntrans
    ind=gh(:,:,i); i
    [Ob(i,:), Am(i,:), C1m(i,:), C2m(i,:), Dmat(i,:)] = AICc_rank(ind,iUY,i_UiUj,UYindex1,Ncons);
end
[Ob, idOb]=sort(Ob);
Utrans = Utrans(idOb,:); 

%%
figure(101)
plot(Ob)

utrans_est = utrans_est(idOb,:);
utrans_fine = utrans_fine(idOb,:);
if Ntrans>80
    nut = 80;                         % No. of top rank Utrans to be considered for BnB
    Ntrans = nut;
    utrans_est = utrans_est(1:nut,:);
    utrans_fine = utrans_fine(1:nut,:);
else
    nut = Ntrans;
end

UYindex2 = [];
for i=1:numel(UYindex1)
    uyi = find(idOb==UYindex1(i));
    if uyi <= nut
        UYindex2 = [UYindex2 uyi];
    end
end

[Uindx2, ~] = find(idOb==Uindx);      % New positions of uis in the ranked basis set

save('data4est','data_est','u_est','utrans_est','Uindx2','UYindex2','sel_indx_all','sel_indx','idOb','nut','mnd','mxd','umin','umax')
clear data_est u_est utrans_est 

save('data4fine','data_fine','u_fine','utrans_fine','idOb','nut')
clear data_fine u_fine utrans_fine gh

%% Branch and Bound for model selection
ndes = 12;                            % Number of desired solution
Nworker = 8;
[AICcve, Amat, C1mat, C2mat, Dmat, Lmat, Obtrnd]=bnbwcsingle(N,ndes,lengthu,Ntrans,Ncons,iUY,i_UiUj,UYindex2,Nworker)


figure
plot(nonzeros(Obtrnd(1,:)))
ylabel('Objective fcn (AICc)')
xlabel('Number of Basis fcns')

%% Post processing of indices of transformed inputs for interpretation
slnrank = 1;
[~, ~] = BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Mass constraints with DDR
for i=1:13
    indx = ~(Obtrnd(i,:)==0);
    [Objfunc(i,:),Amatc(i,:), C1matc(i,:), C2matc(i,:), Dmatc(i,:), lmatc(i,:), ybmlmat(:,:,i), fvalmat(:,i)] = AICc_Constrained(indx,iUY,i_UiUj,UYindex2, Ncons, Amat(i,:), C1mat(i,:), C2mat(i,:), Dmat(i,:), Lmat(i,:));
end
save('solution_workspaceC','N','Nhd','AICcve','Amat','C1mat','C2mat','Dmat','Amatc','C1matc','C2matc','Dmatc','Objfunc','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')

