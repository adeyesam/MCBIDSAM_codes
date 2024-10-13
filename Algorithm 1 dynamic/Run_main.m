%% Load data 
% load('data.mat');
% load('t.mat');
% load('u.mat');

%% 
N = size(data,1);             
N2 = N^2;                     
lengthu = size(u,1);        
lengthu2 = size(u,2);         
lenreading = size(data,2);

%% Adding element bal to data
nci = [1 0 1 0];                    % No. of C in each species i
nhi = [0 2 0 2];                    % Net No. of H in each species i on LHS 
Cf = zeros(numel(nci),lenreading);  % Conc of each spp in feed [Caf Cbf Ccf Cdf]
Cf(:,:) = u;                        % Caf is u2

cbal = zeros(1,lenreading);
hbal = cbal;
for i=1:lenreading
    cbal(i) = nci*Cf(:,i);
    hbal(i) = nhi*Cf(:,i);
end
datawmcon = [data;cbal;hbal];       % Data with mass constr
Nwcons = size(datawmcon,1);
Ncons = Nwcons - N;                 % No. of mass constr imposed

%% Analyses for steady state
sel_indx_all = Identify_SteadyState(data,u);

%%
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
Ntrans = size(Utrans,1)         
UYindex1 = BasisInterpretation(lengthu,iUY,i_UiUj,[],[],[],false);

%% Selecting specific range(s) of data for ranking, model estimation and finetuning
rank_sel =1:200;
est_sel =1:800; 
fine_sel = est_sel;

data_rank = datawmcon(:,rank_sel);
u_rank = u(:,rank_sel);
utrans_rank = Utrans(:,rank_sel);

data_est = datawmcon(:,est_sel);
u_est = u(:,est_sel);
utrans_est = Utrans(:,est_sel);

data_fine = datawmcon(:,fine_sel);
u_fine = u(:,fine_sel);
utrans_fine = Utrans(:,fine_sel);

%% Analyses of training data for steady state
sel_indx = Identify_SteadyState(data(:,est_sel),u);
clear data u Utrans datawmcon

%% Ranking basis functions
save('data4rank','data_rank','u_rank','utrans_rank','UYindex1')
clear data_rank u_rank utrans_rank

gh=zeros(1,Ntrans,Ntrans);
for i=1:Ntrans
    gh(:,i,i)=1;
end

ind=zeros(1,Ntrans);
Nworker = 8;
parfor i=1:Ntrans
    ind=gh(:,:,i); i
    [Ob(i,:), Am(i,:), Cm(i,:)] = AICc_rank(ind,iUY,i_UiUj,UYindex1);
end
[Ob, idOb]=sort(Ob);
figure(101)
plot(Ob)

%% 
utrans_est = utrans_est(idOb,:);
utrans_fine = utrans_fine(idOb,:);
if Ntrans>80
    nut = 80;                          % No. of top rank Utrnas to be considered for BnB
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

save('data4est','data_est','u_est','utrans_est','sel_indx_all','sel_indx','idOb','nut','UYindex2','mnd','mxd')
clear data_est u_est utrans_est u1 u2 u3 uu1 uu2 uu3 ut1 ut2 ut3

save('data4fine','data_fine','u_fine','utrans_fine','idOb','nut','UYindex2')
clear data_fine u_fine utrans_fine

%% Branch and Bound for model selection
ndes = 10;                   % Number of desired solution
Nworker = 4;
[AICcve, Amat, Cmat, Lmat, Obtrnd]=bnbwcsingle(N,ndes,Ntrans,Ncons,iUY,i_UiUj,UYindex2,Nworker)

%% Save solution workspace
save('solution_workspace','N','AICcve','Amat','Cmat','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')

%% Post processing of indices of transformed inputs for interpretation
slnrank = 1;
UYindex = BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Mass constraints with DDR
for i=1:ndes
    indx = ~(Obtrnd(i,:)==0);
    [Objfunc(i,:),Amatc(i,:), Cmatc(i,:), lmatc(i,:), ybmlmat(:,:,i), yddrmat(:,:,i), wmat(:,:,i), fvalmat(:,:,i)] = AICc_Constrained(indx, iUY, i_UiUj, UYindex, Ncons, Amat(i,:), Cmat(i,:), Lmat(i,:),nci,nhi);
end

save('solution_workspaceC','N','Nhd','AICcve','Amat','Cmat','Amatc','Cmatc','Objfunc','est_sel','idOb','Lmat','wmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')

