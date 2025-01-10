%% If you use this algorithm, please cite
% Adeyemo, S., Bhattacharyya, D. Sparse Mass-Constrained Nonlinear Dynamic Model Building from
% Noisy Data Using a Bayesian Approach. Ind. Eng. Chem. Res. 2025 https://doi.org/10.1021/acs.iecr.4c02481

%% Load data 
% load('data.mat');
% load('u.mat');

% The rows of 'data' and 'u' are the measured output and input variables
% respectively while the columns are the time instants. 

%% Dimension of data
N = size(data,1);                                       % No. of output variablees
N2 = N^2;                       
[lengthu, lengthu2] = size(u);                          % Dimension of input data
lenreading = size(data,2);

%% Number of atoms of elements of interest per molecule of species in input/output variables
nci = [0 5 5 5 10];                                     % No. of C in each input variable [F/V Caf, Cbf, Ccf, Cdf]
ncio = [5 5 5 10];                                      % No. of C in each output variable [Ca, Cb, Cc, Cd]
nhi = [0 6 6 6 12];                                     % Net No. of H in each species i on LHS 
nhio = [6 6 6 12];                                      % Net No. of H in each species i on RHS
Cf = zeros(numel(nci),lenreading);                      % Conc of each spp in feed [Caf Cbf Ccf Cdf]
Cf(:,:) = u;                            

% Compute number of moles of C and H atoms at the inlet
cbal = zeros(1,lenreading);
hbal = cbal;
for i=1:lenreading
    cbal(i) = nci*Cf(:,i);
    hbal(i) = nhi*Cf(:,i);
end

% Augment output variables with computed moles of atom for mass balance
datawmcon = [data; cbal; hbal];       
Nwcons = size(datawmcon,1);
Ncons = Nwcons - N;                                     % No. of mass constraints imposed

%% Range of input/output variables
umax=zeros(1,lengthu); umin=umax;
for i=1:lengthu
    umax(i) = max(u(i,:));
    umin(i) = min(u(i,:));
end
for i=1:Nwcons
    mxd(i)=max(datawmcon(i,:));
    mnd(i)=min(datawmcon(i,:));
end

%% Normalization of input/output data
for i=1:N
    datawmcon(i,:)=1 + 2*(datawmcon(i,:)-mnd(i))/(mxd(i)-mnd(i));
end
for i = 1:size(u,1)
    u(i,:) = 1+2*(u(i,:)-umin(i))/(umax(i)-umin(i));
end

%% Compute input transformations
[Utrans, i_UiUj] = InputTransformation(u,[],true);
Ntrans = size(Utrans,1)                                  % No. of transformations
BasisInterpretation(lengthu,i_UiUj,[],[],[],false);

%% Selecting specific range(s) of data for ranking and estimation
rank_sel = randsample(lenreading,0.3*lenreading);
est_sel = randsample(lenreading,0.3*lenreading);

data_rank=datawmcon(1:N,rank_sel);
u_rank=u(:,rank_sel);
utrans_rank=Utrans(:,rank_sel);

data_est=datawmcon(:,est_sel);
u_est=u(:,est_sel);
utrans_est=Utrans(:,est_sel);

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
    [Ob(i,:), Cm(i,:), Dm(i,:)] = RankBasis(ind);
end
[Ob, idOb]=sort(Ob);
Utrans = Utrans(idOb,:); 

%%
utrans_est = utrans_est(idOb,:);
if Ntrans>85
    nut = 85;                                               % No. of top rank Utrans to be considered for BnB
    Ntrans = nut;
    utrans_est = utrans_est(1:nut,:);
else
    nut = Ntrans;
end

save('data4est','data_est','u_est','utrans_est','idOb','nut','mnd','mxd')
clear data_est u_est utrans_est 

%% Branch and Bound for model selection
ndes = 10;                                                   % Number of desired hierarchically ranked models
Nworker = 8;                                                % Number of workers
[AICcve, Cmat, Dmat, Lmat, Obtrnd] = BranchandBound(N,ndes,Ntrans,Ncons,Nworker)

%% Save solution workspace
save('solution_workspace','N','AICcve','Cmat','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj')

%% Post processing of indices of transformed inputs for interpretation
% slnrank = 1;
% BasisInterpretation(lengthu,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Mass constraints with DDR
for i=1:ndes
    indx = ~(Obtrnd(i,:)==0);
    [Objfunc(i,:),Cmatc(i,:), Dmatc(i,:), lmatc(i,:), ybmlmat(:,:,i), yddrmat(:,:,i), wmat(:,:,i), fvalmat(:,:,i)] = ConstrainedEstimation(indx, Ncons, Cmat(i,:), Dmat(i,:), Lmat(i,:), ncio, nhio);
end

%% Save solution workspace
save('solution_workspaceC','N','AICcve','Amat','Cmat','Amatc','Cmatc','Objfunc','est_sel','idOb','Lmat','wmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY')
