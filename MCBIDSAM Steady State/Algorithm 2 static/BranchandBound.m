function [Aic,C1mat,C2mat,Dmat,Lmat,Obtrnd] =BranchandBound(N,ndes,lengthu,Ntrans,Ncons,Nworker)

% ndes is the number of solutions you desire, ie. 10 => returns the global top 10 solutions
% Nworker is the number of workers, ie processors, available to you to solve the problem


%% Base model (with linear basis and bias only)
szlam = N^2 + (N-N^2)/2;
[Obj_base, C1mat_b, C2mat_b, Dmat_b, lmat_b] = ModelEstimation(0, Ncons, zeros(1,N*lengthu),zeros(1,N*Ntrans),zeros(1,N),zeros(1,szlam));
%%

h1j=zeros(1,Ntrans);
timer=0;
szlam = N^2 + (N-N^2)/2;
pct=1e-6;

% h1j is a vector which you use to for your branching strategy
h1j=1:Ntrans;
flag1=true;
flag2=true;

% F is the fixed set of a node, C is the candidate set of a node, both vectors of logical 1s and 0s
F=false(1,Ntrans);
C=true(1,Ntrans);

% ASTB is your initial incumbent, ie perfection
ASTB=inf; C1MATB=zeros(1,N*lengthu); C2MATB=zeros(1,N*Ntrans); DMATB=zeros(1,N); LMATB=zeros(1,szlam); OBTRNDB = zeros(1,Ntrans);
% ASTB=Obj_base; C1MATB=C1mat_b; C2MATB=C2mat_b; DMATB=Dmat_b; LMATB=lmat_b; OBTRNDB = zeros(1,Ntrans);


ita=0;
Aic=Inf(ndes,1);
% Aic=Obj_base*ones(ndes,1);
C1mat=zeros(ndes,N*lengthu);
C2mat=zeros(ndes,N*Ntrans);
Dmat=zeros(ndes,N);
Lmat=zeros(ndes,szlam);
Obtrnd=zeros(ndes,Ntrans);
dumbsolve();

Aic=[Aic ; Obj_base];
C1mat=[C1mat ; C1mat_b];
C2mat=[C2mat ; C2mat_b];
Dmat=[Dmat ; Dmat_b];
Lmat=[Lmat ; lmat_b];
Obtrnd=[Obtrnd ; zeros(1,Ntrans)];

    function dumbsolve()
        while ~isempty(F)
            ita=ita+1;
            Abound=max(Aic);
            %Pruning of supernodes
            if ita>2
                [bumb,~]=size(F);
                discard=false(bumb,1);
                L1=ASTB>Abound;
                discard=find(L1);
                if any(discard)
                    F(discard,:)=[];
                    C(discard,:)=[];
                    ASTB(discard,:)=[];
                    C1MATB(discard,:)=[];
                    C2MATB(discard,:)=[];
                    DMATB(discard,:)=[];
                    LMATB(discard,:)=[];
                    OBTRNDB(discard,:)=[];
                end
            end
            if ~isempty(F)    
                %Tracks Progress
                FC=[];
                CC=[];
                ASTBT=[];
                C1MATBT=[];
                C2MATBT=[];
                DMATBT=[];
                LMATBT=[];
                OBTRNDBT=[];
                [NumberofNodes,~]=size(F)
                t=0;
                if NumberofNodes>28000
                    flagconst=true;
                else
                    flagconst=false;
                end
            end   
            while NumberofNodes<Nworker && ~isempty(F)
                [~,idc]=max(ASTB);F               
                f1=F(idc,:)
                c1=C(idc,:)
                f2=f1;
                h1k=h1j;
                h1k(~c1)=Inf
                [~,id]=min(h1k);
                f2(id)=true
                c1(id)=false
                F(idc,:)=[];
                C(idc,:)=[];
                [a1,c1mat1,c2mat1,dmat1,lmat1]=feval(@Gc,f2,Ntrans,Ncons,lengthu,ASTB(idc,:),C1MATB(idc,:),C2MATB(idc,:),DMATB(idc,:),LMATB(idc,:),N,true,false)
                [a2,c1mat2,c2mat2,dmat2,lmat2]=feval(@Gc,f1,Ntrans,Ncons,lengthu,ASTB(idc,:),C1MATB(idc,:),C2MATB(idc,:),DMATB(idc,:),LMATB(idc,:),N,false,false)
                
                obtrndp = OBTRNDB(idc,:);
                
                if (ita==1)||(sum(f2)==0)
                    if a1 <= ASTB(idc,:)
                        ASTB = [ASTB; a1];
                        F=[F ; f2];
                        C=[C ; c1];
                        C1MATB=[C1MATB ; c1mat1];
                        C2MATB=[C2MATB ; c2mat1];
                        DMATB=[DMATB ; dmat1];
                        LMATB=[LMATB ; lmat1];
                        obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                        OBTRNDB=[OBTRNDB ; obtrnd1];                    
                    end
                else
                    if a1 < ASTB(idc,:)
                        ASTB = [ASTB; a1];
                        F=[F ; f2];
                        C=[C ; c1];
                        C1MATB=[C1MATB ; c1mat1];
                        C2MATB=[C2MATB ; c2mat1];
                        DMATB=[DMATB ; dmat1];
                        LMATB=[LMATB ; lmat1];
                        obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                        OBTRNDB=[OBTRNDB ; obtrnd1];                    
                    end
                end
                    if sum(c1)>0
                        ASTB = [ASTB; a2];
                        F=[F ; f1];
                        C=[C ; c1];
                        C1MATB=[C1MATB ; c1mat2]; 
                        C2MATB=[C2MATB ; c2mat2]; 
                        DMATB=[DMATB ; dmat2]; 
                        LMATB=[LMATB ; lmat2]; 
                        obtrnd2=obtrndp; %obtrnd(:,find(f1,1,'last'))=a2;
                        OBTRNDB=[OBTRNDB ; obtrnd2];
                    end

                if isempty(F) && sum(c1)==0 && ASTB(idc,:)>max(Aic)
                    Aic=[Aic ; a2];
                    C1mat=[C1mat ; c1mat2];
                    C2mat=[C2mat ; c2mat2];
                    Dmat=[Dmat ; dmat2];
                    Lmat=[Lmat ; lmat2];
                    Obtrnd=[Obtrnd ; obtrnd2];
                end

                ASTB(idc,:)=[];
                C1MATB(idc,:)=[];
                C2MATB(idc,:)=[];
                DMATB(idc,:)=[];
                LMATB(idc,:)=[];
                OBTRNDB(idc,:)=[];
                
                [NumberofNodes,~]=size(F);
            end
            if ~isempty(F)  
                iterat=inf;
                if ~flagconst
                    if NumberofNodes<100*Nworker
                        FC=F;
                        F=[];
                        CC=C;
                        C=[];
                        ASTBT=ASTB;
                        ASTB=[];
                        C1MATBT=C1MATB;
                        C1MATB=[];
                        C2MATBT=C2MATB;
                        C2MATB=[];
                        DMATBT=DMATB;
                        DMATB=[];
                        LMATBT=LMATB;
                        LMATB=[];
                        OBTRNDBT=OBTRNDB;
                        OBTRNDB=[];
                    else
                        shf=randperm(NumberofNodes,100*Nworker);
                        FC=F(shf,:);
                        F(shf,:)=[];
                        CC=C(shf,:);
                        C(shf,:)=[];
                        ASTBT=ASTB(shf,:);
                        ASTB(shf,:)=[];
                        C1MATBT=C1MATB(shf,:);
                        C1MATB(shf,:)=[];
                        C2MATBT=C2MATB(shf,:);
                        C2MATB(shf,:)=[];
                        DMATBT=DMATB(shf,:);
                        DMATB(shf,:)=[];
                        LMATBT=LMATB(shf,:);
                        LMATB(shf,:)=[];
                        OBTRNDBT=OBTRNDB(shf,:);
                        OBTRNDB(shf,:)=[];
                    end
                else
                    [~,shf]=sort(sum(C,2));
                    F=F(shf,:);
                    C=C(shf,:);
                    ASTB=ASTB(shf,:);
                    C1MATB=C1MATB(shf,:);
                    C2MATB=C2MATB(shf,:);
                    LMATB=LMATB(shf,:);
                    OBTRNDB=OBTRNDB(shf,:);
                    FC=F(1:100*Nworker,:);
                    F(1:100*Nworker,:)=[];
                    CC=C(1:100*Nworker,:);
                    C(1:100*Nworker,:)=[];
                    ASTBT=ASTB(1:100*Nworker,:);
                    ASTB(1:100*Nworker,:)=[];
                    C1MATBT=C1MATB(1:100*Nworker,:);
                    C1MATB(1:100*Nworker,:)=[];
                    C2MATBT=C2MATB(1:100*Nworker,:);
                    C2MATB(1:100*Nworker,:)=[];
                    DMATBT=DMATB(1:100*Nworker,:);
                    DMATB(1:100*Nworker,:)=[];
                    LMATBT=LMATB(1:100*Nworker,:);
                    LMATB(1:100*Nworker,:)=[];
                    OBTRNDBT=OBTRNDB(1:100*Nworker,:);
                    OBTRNDB(1:100*Nworker,:)=[];
                end
                [NodestoWorker,~]=size(FC)
                Li=randperm(NodestoWorker); 
                FC=FC(Li,:);
                CC=CC(Li,:);  
                ASTBT=ASTBT(Li,:);
                C1MATBT=C1MATBT(Li,:);
                C2MATBT=C2MATBT(Li,:);
                DMATBT=DMATBT(Li,:);
                LMATBT=LMATBT(Li,:);
                OBTRNDBT=OBTRNDBT(Li,:);
                divid=floor(NodestoWorker/Nworker);
                for tt=1:Nworker
                    if tt==1
                        t=1;
                    else
                        t=t+divid;
                    end
                    if tt==Nworker
                        problems(tt).F=FC(t:NodestoWorker,:);
                        problems(tt).C=CC(t:NodestoWorker,:);
                        problems(tt).A=ASTBT(t:NodestoWorker,:);
                        problems(tt).C1MATB=C1MATBT(t:NodestoWorker,:);
                        problems(tt).C2MATB=C2MATBT(t:NodestoWorker,:);
                        problems(tt).DMATB=DMATBT(t:NodestoWorker,:);
                        problems(tt).LMATB=LMATBT(t:NodestoWorker,:);
                        problems(tt).OBTRNDB=OBTRNDBT(t:NodestoWorker,:);
                    else
                        problems(tt).F=FC(t:(divid*tt),:);
                        problems(tt).C=CC(t:(divid*tt),:);
                        problems(tt).A=ASTBT(t:(divid*tt),:);
                        problems(tt).C1MATB=C1MATBT(t:(divid*tt),:);
                        problems(tt).C2MATB=C2MATBT(t:(divid*tt),:);
                        problems(tt).DMATB=DMATBT(t:(divid*tt),:);
                        problems(tt).LMATB=LMATBT(t:(divid*tt),:);
                        problems(tt).OBTRNDB=OBTRNDBT(t:(divid*tt),:);
                    end
                end

                %This section is solved in parallel
                parfor z=1:Nworker
                  if ~isempty(problems(z).F)  
                    FT=problems(z).F
                    CT=problems(z).C;
                    AT=Inf;
                    AST=problems(z).A;
                    C1MATST=problems(z).C1MATB;
                    C2MATST=problems(z).C2MATB;
                    DMATST=problems(z).DMATB;
                    LMATST=problems(z).LMATB;
                    OBTRNDST=problems(z).OBTRNDB;
                    C1MATT=zeros(1,N*lengthu);
                    C2MATT=zeros(1,N*Ntrans);
                    DMATT=zeros(1,N);
                    LMATT=zeros(1,szlam);
                    OBTRNDT=zeros(1,Ntrans);
                    DUMBA=Aic;
                    timer=0;
                    termin=true;
                    toccers=0;
                    while termin
                        tic
                        %chooses which node to branch
                        timer=timer+1
                        if flagconst
                            [~,idn]=min(sum(CT,2));
                        else
                            if toccers<90
                                [~,idn]=min(AST);
                            else
                                [~,idn]=min(sum(CT,2));
                            end
                        end
                        f1=FT(idn,:)
                        c1=CT(idn,:)
                        AS=AST(idn,:);
                        C1MATS=C1MATST(idn,:);
                        C2MATS=C2MATST(idn,:);
                        DMATS=DMATST(idn,:);
                        LMATS=LMATST(idn,:);
                        OBTRNDS=OBTRNDST(idn,:);
                        f2=f1;

                        %chooses which variable to branch
                        h1k=h1j;
                        h1k(~c1)=Inf
                        [~,id]=min(h1k);
                        f2(id)=true
                        c1(id)=false
                        c2=c1;
                        OBTRNDST_idn=OBTRNDST(idn,:);

                        %removing node to be branched
                        FT(idn,:)=[];
                        CT(idn,:)=[];
                        AST(idn,:)=[];
                        C1MATST(idn,:)=[];
                        C2MATST(idn,:)=[];
                        DMATST(idn,:)=[];
                        LMATST(idn,:)=[];
                        OBTRNDST(idn,:)=[];
                        flag1=true;
                        flag2=true;

                        %checking whether node1 is terminal
                        cc=sum(c1);
                        if flag1 && ~isempty(OBTRNDST_idn)
                              [a,c1mat,c2mat,dmat, lmat]=feval(@Gc,f2,Ntrans,Ncons,lengthu,0,zeros(1,N*lengthu),zeros(1,N*Ntrans),zeros(1,N),zeros(1,szlam),N,true,true);
                              Abound=min([max(AT) max(DUMBA)]);
                              if a<AS %Abound
                                percent=(AS-a)/abs(AS)
                                    if (cc==0 || percent<= pct)
                                        AT=[a ; AT];
                                        C1MATT=[c1mat ; C1MATT];
                                        C2MATT=[c2mat ; C2MATT];
                                        DMATT=[dmat ; DMATT];
                                        LMATT=[lmat ; LMATT];
                                        DUMBA=[a ; DUMBA];
                                        obtrnd=OBTRNDS; obtrnd(:,id)=a;
                                        OBTRNDT=[obtrnd ; OBTRNDT];

                                        if size(AT,1)>ndes
                                            [~,IDDD]=max(AT);
                                            AT(IDDD,:)=[];
                                            C1MATT(IDDD,:)=[];
                                            C2MATT(IDDD,:)=[];
                                            DMATT(IDDD,:)=[];
                                            LMATT(IDDD,:)=[];
                                            OBTRNDT(IDDD,:)=[];
                                        end

                                        flag1=false;                           
                                    end
                                end
                        end

                        %checking whether node2 is terminal
                        cc=sum(c2)
                       
                            if flag2 && cc==0
                                    AT=[AS ; AT];
                                    C1MATT=[C1MATS ; C1MATT];
                                    C2MATT=[C2MATS ; C2MATT];
                                    DMATT=[DMATS ; DMATT];
                                    LMATT=[LMATS ; LMATT];
                                    DUMBA=[AS ; DUMBA];
                                    OBTRNDT=[OBTRNDS ; OBTRNDT];
                                    
                                    if size(AT,1)>ndes
                                        [~,IDDD]=max(AT);
                                        AT(IDDD,:)=[];
                                        C1MATT(IDDD,:)=[];
                                        C2MATT(IDDD,:)=[];
                                        DMATT(IDDD,:)=[];
                                        LMATT(IDDD,:)=[];
                                        OBTRNDT(IDDD,:)=[];
                                    end
                                    
                                    flag2=false;
                            end
                        


                        %Pruning of Branch 1
                        if (flag1 && cc>0)
                            [a,c1mat,c2mat,dmat, lmat]=feval(@Gc,f2,Ntrans,Ncons,lengthu,AS,C1MATS,C2MATS,DMATS,LMATS,N,true,false);
                            Abound=min([max(AT) max(DUMBA)]);
                            if a<AS %Abound
                                AST=[AST ; a];
                                FT=[FT ; f2];
                                CT=[CT ; c1];
                                C1MATST=[C1MATST ; c1mat];
                                C2MATST=[C2MATST ; c2mat];
                                DMATST=[DMATST ; dmat];
                                LMATST=[LMATST ; lmat];

                                obtrnd=OBTRNDS; obtrnd(:,id)=a;
                                OBTRNDST=[OBTRNDST ; obtrnd];
                            else
                                flag1=false;
                            end
                        end

                        %Pruning of Branch 2
                        if (flag2 && cc>0)
                                AST=[AST ; AS];                            
                                FT=[FT ; f1];
                                CT=[CT ; c2];
                                C1MATST=[C1MATST ; C1MATS];
                                C2MATST=[C2MATST ; C2MATS];
                                DMATST=[DMATST ; DMATS];
                                LMATST=[LMATST ; LMATS];

                                obtrnd=OBTRNDS; 
                                OBTRNDST=[OBTRNDST ; obtrnd];
                        else
                                flag2=false;
                        end
                        dtt=toc;
                        toccers=toccers+dtt;
                        if timer>iterat || isempty(FT) || toccers>300
                            termin=false;
                        end
                   end
                    %Storing existing nodes and solutions from iteration
                    [ttte,~]=size(FT);
                    if ttte>0
                        F=[F ; FT];
                        C=[C ; CT];
                        ASTB=[ASTB ; AST];
                        C1MATB=[C1MATB ; C1MATST];
                        C2MATB=[C2MATB ; C2MATST];
                        DMATB=[DMATB ; DMATST];
                        LMATB=[LMATB ; LMATST];
                        OBTRNDB=[OBTRNDB ; OBTRNDST];
                    end
                    SOLUTIONS(z).D2=AT;
                    SOLUTIONS(z).D3=C1MATT;
                    SOLUTIONS(z).D4=C2MATT;
                    SOLUTIONS(z).D5=DMATT;
                    SOLUTIONS(z).D6=LMATT;
                    SOLUTIONS(z).D7=OBTRNDT;
                  end  
                end

                %Post-processing of solutions
                [~,st]=size(SOLUTIONS)
                for i=1:st
                    Aic=[Aic ; SOLUTIONS(i).D2];
                    C1mat=[C1mat ; SOLUTIONS(i).D3];
                    C2mat=[C2mat ; SOLUTIONS(i).D4];
                    Dmat=[Dmat ; SOLUTIONS(i).D5];
                    Lmat=[Lmat ; SOLUTIONS(i).D6];
                    Obtrnd=[Obtrnd ; SOLUTIONS(i).D7];
                end
            end
            if size(Aic,1)>ndes
                [Aic,idpd]=sort(Aic);
                C1mat=C1mat(idpd,:);
                C2mat=C2mat(idpd,:);
                Dmat=Dmat(idpd,:);
                Lmat=Lmat(idpd,:);
                Obtrnd=Obtrnd(idpd,:);
                idsss=size(Aic);
                Aic((ndes+1):idsss,:)=[];
                C1mat((ndes+1):idsss,:)=[];
                C2mat((ndes+1):idsss,:)=[];
                Dmat((ndes+1):idsss,:)=[];
                Lmat((ndes+1):idsss,:)=[];
                Obtrnd((ndes+1):idsss,:)=[];
            end
            SOLUTIONS=[];
        end
    end


end

%This is your objective function
%flag2 defines rather the node is terminal or not, logical 1 or 0
%flag defines rather the node has been upwardly branched or downwardly branched, 1 == upward branched, 0 == downward branched
function [a,c1mat,c2mat,dmat, lmat]=Gc(f,cv,Ncons,lengthu,as,c1mats,c2mats,dmats,lmats,n,flag,flag2)
szlam = n^2 + (n-n^2)/2;
k=cv-sum(f);
warning('off', 'all')
if flag2
        [a, c1mat, c2mat, dmat, lmat] = ModelEstimation(f,Ncons,c1mats,c2mats,dmats,lmats);
else
    if k==cv
        a=inf; c1mat=zeros(1,n*lengthu); c2mat=zeros(1,n*cv);dmat=zeros(1,n);  lmat=zeros(1,szlam);
    elseif flag
        [a, c1mat, c2mat, dmat, lmat] = ModelEstimation(f,Ncons,c1mats,c2mats,dmats,lmats);
    else
        a=as; c1mat=c1mats; c2mat=c2mats; dmat=dmats; lmat=lmats;
    end
end
end