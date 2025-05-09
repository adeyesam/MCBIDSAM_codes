function [Aic,Cmat,Dmat,Lmat,Obtrnd]=BranchandBound(N,ndes,Ntrans,Ncons,Nworker)

% ndes is the number of solutions you desire, ie. 10 => returns the global top 10 solutions
% Nworker is the number of workers, ie processors, available to you to solve the problem

timer=0;
szlam = N^2 + (N-N^2)/2;
pct = 1e-6;

% h1j is a vector which you use to for your branching strategy
h1j=1:Ntrans;
flag1=true;
flag2=true;

% F is the fixed set of a node, C is the candidate set of a node, both vectors of logical 1s and 0s
F=false(1,Ntrans);
C=true(1,Ntrans);

% ASTB is your initial incumbent, ie perfection
ASTB=inf; CMATB=zeros(1,N*Ntrans); DMATB=zeros(1,N); LMATB=zeros(1,szlam); OBTRNDB = zeros(1,Ntrans);


ita=0;
Aic=Inf(ndes,1);
Cmat=zeros(ndes,N*Ntrans);
Dmat=zeros(ndes,N);
Lmat=zeros(ndes,szlam);
Obtrnd=zeros(ndes,Ntrans);
dumbsolve();

    function dumbsolve()
        while ~isempty(F)
            ita=ita+1;
            Abound=max(Aic);
            %Pruning of supernodes
            if ita>5
                [bumb,~]=size(F);
                discard=false(bumb,1);
                L1=ASTB>Abound;
                discard=find(L1);
                if any(discard)
                    F(discard,:)=[];
                    C(discard,:)=[];
                    ASTB(discard,:)=[];
                    CMATB(discard,:)=[];
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
                CMATBT=[];
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
                [a1,cmat1,dmat1,lmat1]=feval(@Gc,f2,Ntrans,Ncons,ASTB(idc,:),CMATB(idc,:),DMATB(idc,:),LMATB(idc,:),N,true,false)
                [a2,cmat2,dmat2,lmat2]=feval(@Gc,f1,Ntrans,Ncons,ASTB(idc,:),CMATB(idc,:),DMATB(idc,:),LMATB(idc,:),N,false,false)
                
                obtrndp = OBTRNDB(idc,:);
                
                if (ita==1)||(sum(f2)==0)
                    if a1 <= ASTB(idc,:)
                        ASTB = [ASTB; a1];
                        F=[F ; f2];
                        C=[C ; c1];CMATB
                        CMATB=[CMATB ; cmat1];
                        DMATB=[DMATB ; dmat1];
                        LMATB=[LMATB ; lmat1];
                        obtrnd1=obtrndp; obtrnd1(:,id)=a1;
                        OBTRNDB=[OBTRNDB ; obtrnd1];                    
                    end
                else
                    if a1 < ASTB(idc,:)
                        ASTB = [ASTB; a1];
                        F=[F ; f2];
                        C=[C ; c1];CMATB
                        CMATB=[CMATB ; cmat1];
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
                        CMATB=[CMATB ; cmat2]; 
                        DMATB=[DMATB ; dmat2]; 
                        LMATB=[LMATB ; lmat2]; 
                        obtrnd2=obtrndp; %obtrnd(:,find(f1,1,'last'))=a2;
                        OBTRNDB=[OBTRNDB ; obtrnd2];
                    end

                if isempty(F) && sum(c1)==0 && ASTB(idc,:)>max(Aic)
                    Aic=[Aic ; a2];
                    Cmat=[Cmat ; cmat2];
                    Dmat=[Dmat ; dmat2];
                    Lmat=[Lmat ; lmat2];
                    Obtrnd=[Obtrnd ; obtrnd2];
                end

                ASTB(idc,:)=[];
                CMATB(idc,:)=[];
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
                        CMATBT=CMATB;
                        CMATB=[];
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
                        CMATBT=CMATB(shf,:);
                        CMATB(shf,:)=[];
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
                    CMATB=CMATB(shf,:);
                    LMATB=LMATB(shf,:);
                    OBTRNDB=OBTRNDB(shf,:);
                    FC=F(1:100*Nworker,:);
                    F(1:100*Nworker,:)=[];
                    CC=C(1:100*Nworker,:);
                    C(1:100*Nworker,:)=[];
                    ASTBT=ASTB(1:100*Nworker,:);
                    ASTB(1:100*Nworker,:)=[];
                    CMATBT=CMATB(1:100*Nworker,:);
                    CMATB(1:100*Nworker,:)=[];
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
                CMATBT=CMATBT(Li,:);
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
                        problems(tt).CMATB=CMATBT(t:NodestoWorker,:);
                        problems(tt).DMATB=DMATBT(t:NodestoWorker,:);
                        problems(tt).LMATB=LMATBT(t:NodestoWorker,:);
                        problems(tt).OBTRNDB=OBTRNDBT(t:NodestoWorker,:);
                    else
                        problems(tt).F=FC(t:(divid*tt),:);
                        problems(tt).C=CC(t:(divid*tt),:);
                        problems(tt).A=ASTBT(t:(divid*tt),:);
                        problems(tt).CMATB=CMATBT(t:(divid*tt),:);
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
                    CMATST=problems(z).CMATB;
                    DMATST=problems(z).DMATB;
                    LMATST=problems(z).LMATB;
                    OBTRNDST=problems(z).OBTRNDB;
                    CMATT=zeros(1,N*Ntrans);
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
                        CMATS=CMATST(idn,:);
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
                        CMATST(idn,:)=[];
                        DMATST(idn,:)=[];
                        LMATST(idn,:)=[];
                        OBTRNDST(idn,:)=[];
                        flag1=true;
                        flag2=true;

                        %checking whether node1 is terminal
                        cc=sum(c1);
                        if flag1 && ~isempty(OBTRNDST_idn)
                              [a,cmat,dmat, lmat]=feval(@Gc,f2,Ntrans,Ncons,0,zeros(1,N*Ntrans),zeros(1,N),zeros(1,szlam),N,true,true);
                              Abound=min([max(AT) max(DUMBA)]);
                              if a<AS %Abound
                                percent=(AS-a)/abs(AS)
                                    if (cc==0 || percent<= pct)
                                        AT=[a ; AT];
                                        CMATT=[cmat ; CMATT];
                                        DMATT=[dmat ; DMATT];
                                        LMATT=[lmat ; LMATT];
                                        DUMBA=[a ; DUMBA];
                                        obtrnd=OBTRNDS; obtrnd(:,id)=a;
                                        OBTRNDT=[obtrnd ; OBTRNDT];

                                        if size(AT,1)>ndes
                                            [~,IDDD]=max(AT);
                                            AT(IDDD,:)=[];
                                            CMATT(IDDD,:)=[];
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
                                CMATT=[CMATS ; CMATT];
                                DMATT=[DMATS ; DMATT];
                                LMATT=[LMATS ; LMATT];
                                DUMBA=[AS ; DUMBA];
                                OBTRNDT=[OBTRNDS ; OBTRNDT];
                                
                                if size(AT,1)>ndes
                                    [~,IDDD]=max(AT);
                                    AT(IDDD,:)=[];
                                    CMATT(IDDD,:)=[];
                                    DMATT(IDDD,:)=[];
                                    LMATT(IDDD,:)=[];
                                    OBTRNDT(IDDD,:)=[];
                                end
                                flag2=false;
                       end
                        
                        %Pruning of Branch 1
                        if (flag1 && cc>0)
                            [a,cmat,dmat, lmat]=feval(@Gc,f2,Ntrans,Ncons,AS,CMATS,DMATS,LMATS,N, true,false);
                            Abound=min([max(AT) max(DUMBA)]);
                            if a<AS %Abound
                                AST=[AST ; a];
                                FT=[FT ; f2];
                                CT=[CT ; c1];
                                CMATST=[CMATST ; cmat];
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
                                CMATST=[CMATST ; CMATS];
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
                        CMATB=[CMATB ; CMATST];
                        DMATB=[DMATB ; DMATST];
                        LMATB=[LMATB ; LMATST];
                        OBTRNDB=[OBTRNDB ; OBTRNDST];
                    end
                    SOLUTIONS(z).D2=AT;
                    SOLUTIONS(z).D3=CMATT;
                    SOLUTIONS(z).D4=DMATT;
                    SOLUTIONS(z).D5=LMATT;
                    SOLUTIONS(z).D6=OBTRNDT;
                  end  
                end

                %Post-processing of solutions
                [~,st]=size(SOLUTIONS)
                for i=1:st
                    Aic=[Aic ; SOLUTIONS(i).D2];
                    Cmat=[Cmat ; SOLUTIONS(i).D3];
                    Dmat=[Dmat ; SOLUTIONS(i).D4];
                    Lmat=[Lmat ; SOLUTIONS(i).D5];
                    Obtrnd=[Obtrnd ; SOLUTIONS(i).D6];
                end
            end
            if size(Aic,1)>ndes
                [Aic,idpd]=sort(Aic);
                Cmat=Cmat(idpd,:);
                Dmat=Dmat(idpd,:);
                Lmat=Lmat(idpd,:);
                Obtrnd=Obtrnd(idpd,:);
                idsss=size(Aic);
                Aic((ndes+1):idsss,:)=[];
                Cmat((ndes+1):idsss,:)=[];
                Dmat((ndes+1):idsss,:)=[];
                Lmat((ndes+1):idsss,:)=[];
                Obtrnd((ndes+1):idsss,:)=[];
            end
            SOLUTIONS=[];
        end
    end


end

%This is your objective function
function [a,cmat,dmat, lmat]=Gc(f,cv,Ncons,as,cmats,dmats,lmats,n,flag,flag2)
szlam = n^2 + (n-n^2)/2;
k=cv-sum(f);
warning('off', 'all')
if flag2
        [a, cmat, dmat, lmat]=ModelEstimation(f,Ncons,cmats,dmats,lmats);
else
    if k==cv
        a=inf; cmat=zeros(1,n*cv); dmat=zeros(1,n);  lmat=zeros(1,szlam);
    elseif flag
        [a, cmat, dmat, lmat]=ModelEstimation(f,Ncons,cmats,dmats,lmats);
    else
        a=as; cmat=cmats; dmat=dmats; lmat=lmats;
    end
end
end