%% If you use this algorithm, please cite
% Adeyemo, S., Bhattacharyya, D. Sparse Mass-Constrained Nonlinear Dynamic Model Building from
% Noisy Data Using a Bayesian Approach. Ind. Eng. Chem. Res. 2025 https://doi.org/10.1021/acs.iecr.4c02481

%% Load data 
% load('data.mat');
% load('u.mat');

%% Dimension of data
N = size(data,1);                                       % No. of output variablees
N2 = N^2;                       
[lengthu, lengthu2] = size(u);                          % Dimension of input data
lenreading = size(data,2)-1;
lenreadpone = lenreading+1;

%% Number of atoms of elements of interest per molecule of species in input/output variables
nci = [1 0 1 0];                % No. of C in each input variable 
nhi = [0 2 0 2];                % Net No. of H in each species i on LHS 
ncio = [1 0 1 0];               % No. of C in each input variable
nhio = [0 2 0 2];               % Net No. of H in each species i on RHS

% Compute number of moles of C and H atoms at the inlet
cbal = zeros(1,lenreadpone);
hbal = cbal;
for i=1:lenreadpone
    cbal(i) = nci*u(:,i);
    hbal(i) = nhi*u(:,i);
end

% Augment output variables with computed moles of atom for mass balance
datawmcon = [data;cbal;hbal];    % Data with mass constr
Nwcons = size(datawmcon,1);
Ncons = Nwcons - N;              % No. of mass constr imposed

%% Analyses for steady state
steady_indx_all = Identify_SteadyState(data,u);

%% Range of input/output variables
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

%% Normalization of input/output data
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

data_rank=datawmcon(1:N,rank_sel);
u_rank=u(:,rank_sel);
utrans_rank=Utrans(:,rank_sel);

data_est=datawmcon(:,est_sel);
u_est=u(:,est_sel);
utrans_est=Utrans(:,est_sel);

%% Analyses of training data for steady state
steady_indx = Identify_SteadyState(data(:,est_sel),u);

%% Ranking basis functions
save('data4rank','data_rank','u_rank','utrans_rank','mnd','mxd')
clear data_rank u_rank utrans_rank

gh=zeros(1,Ntrans,Ntrans);
for i=1:Ntrans
    gh(:,i,i)=1;
end

ind=zeros(1,Ntrans);
for i=1:Ntrans
    ind=gh(:,:,i);
    [Ob(i,:), Am(i,:), C1m(i,:), C2m(i,:), Dmat(i,:)] = RankBasis(ind,iUY,i_UiUj,UYindex1);
end
[Ob, idOb]=sort(Ob);
Utrans = Utrans(idOb,:); 

%%
utrans_est = utrans_est(idOb,:);
if Ntrans>85
    nut = 85;                         % No. of top rank Utrans to be considered for BnB
    Ntrans = nut;
    utrans_est = utrans_est(1:nut,:);
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

save('data4est','data_est','u_est','utrans_est','Uindx2','UYindex2','steady_indx_all','steady_indx','idOb','nut','mnd','mxd','umin','umax')
clear data_est u_est utrans_est 

%% Branch and Bound for model selection
ndes = 10;                            % Number of desired solution
Nworker = 8;
[AICcve, Amat, C1mat, C2mat, Dmat, Lmat, Obtrnd] = BranchandBound(N,ndes,lengthu,Ntrans,Ncons,iUY,i_UiUj,UYindex2,Nworker)

%% Save solution workspace
save('solution_workspaceC','N','AICcve','Amat','C1mat','C2mat','Dmat','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')

%% Post processing of indices of transformed inputs for interpretation
% slnrank = 1;
% [~, ~] = BasisInterpretation(lengthu,iUY,i_UiUj,Obtrnd,idOb,slnrank,true);

%% Mass constraints with DDR
for i=1:ndes
    indx = ~(Obtrnd(i,:)==0);
    [Objfunc(i,:),Amatc(i,:), C1matc(i,:), C2matc(i,:), Dmatc(i,:), lmatc(i,:), fvalmat(:,i)] = ConstrainedEstimation(indx,iUY,i_UiUj,UYindex2, Ncons, Amat(i,:), C1mat(i,:), C2mat(i,:), Dmat(i,:), Lmat(i,:), nci, ncio, nhi, nhio);
end

%% Save solution workspace
save('solution_workspaceC','N','AICcve','Amat','C1mat','C2mat','Dmat','Amatc','C1matc','C2matc','Dmatc','Objfunc','est_sel','idOb','Lmat','mnd','mxd','umax','umin','Ntrans','Obtrnd','i_UiUj','iUY','UYindex2')
