function PF=MORBF(nval,numcand,Nmax)
 
%nval is the number of dimensions
%numcand is the number of candidate points
%Nmax is the max number of iterations
%MORBF uses the RBF algorithm for a known function, ZDT4
%referenced in the Deb et al (2002) paper, but a different function can be %used in its stead (ZDT4 specific lines are noted below)

%generate initial parameters 
J=nval;
n0=2*J+1;   %doubling the number of starting points to cover more pts
n=0;
fmax=200;
rmax=3;   %this is because i have a higher starting ro_n
PF=zeros(1,3+J);    %PF, f1, f2, J dimensions of Xset
PF(1,2:3)=inf;
 
while n<Nmax
    cycle=n0;
    Cfail=0;
    rcount=0;
    ro_n=ones(J,1)+[0; 2*ones(J-1,1)-(ones(J-1,1))];  %ZDT4
    n=n+n0;
    NumPF=size(PF,1);
    F1=zeros(n0,1);
    F2=zeros(n0,1);
    
    %Generate initial values via LHS method (unconstrained)
    L=lhsu(zeros(J,1),ones(J,1),n0);
    Xset=L';
    Xset(2:J,:)=Xset(2:J,:)*10-5;
 
    %Generating initial Pareto filter
    for i=1:cycle
    %ZDT4 function	
        F1(i,1)=Xset(1,i);
        gx=1+10*(J-1)+sum(Xset(2:J,i).^2-10*cos(4*pi*Xset(2:J,i)));
        F2(i,1)=gx*(1-sqrt(Xset(1,i)/gx));
        %%update Pareto set
        PF(NumPF+1,:)=[1 F1(i,1) F2(i,1) Xset(:,i)']; 
        stop=0;
        nu=1;
        while and(nu<=NumPF,stop==0)
            if or(and(F1(i,1)<PF(nu,2),F2(i,1)<=PF(nu,3)),and(F1(i,1)<=PF(nu,2),F2(i,1)<PF(nu,3)))
                PF(nu,1)=0;
            elseif or(and(F1(i,1)>=PF(nu,2),F2(i,1)>PF(nu,3)),and(F1(i,1)>PF(nu,2),F2(i,1)>=PF(nu,3)))
                PF(NumPF+1,1)=0;
                stop=1;
            end
            nu=nu+1;
        end
        NumPF=size(PF,1);
    end
    oldPF=PF;
    PF=0;
    for i=1:NumPF
        if oldPF(i,1)==1
            if PF==0
                PF=oldPF(i,:);
            else
                PF=[PF; oldPF(i,:)];
            end
        end
    end
    NumPF=size(PF,1);

    %%Now do interpolation
    while and(rcount<=rmax,n<Nmax)
        %%let's try using thin-plate splines    
        if cycle==n0  
            phi=zeros(cycle,cycle);
            for i=1:cycle
                for j=1:cycle
                    if Xset(:,i)-Xset(:,j)==0
                        phi(i,j)=0;
                    else
                        phi(i,j)=sum((Xset(:,i)-Xset(:,j)).^2)*log(sqrt(sum((Xset(:,i)-Xset(:,j)).^2))+realmin);  
                    end
                end
            end
            poly=zeros(cycle,J+1);
            for i=1:cycle
                for j=1:J
                    poly(i,j)=Xset(j,i);
                end
                poly(i,J+1)=1;
            end    
 
        else
            phiprime=zeros(cycle,cycle);
            phiprime(1:cycle-1,1:cycle-1)=phi;
            phi=phiprime;
            for i=1:cycle-1
                if Xset(:,i)-Xset(:,cycle)==0
                    phi(i,cycle)=0;
                else
                    phi(i,cycle)=sum((Xset(:,i)-Xset(:,cycle)).^2)*log(sqrt(sum((Xset(:,i)-Xset(:,cycle)).^2))+realmin);  
                end
            end
            for j=1:cycle
                if Xset(:,cycle)-Xset(:,j)==0
                    phi(cycle,j)=0;
                else
                    phi(cycle,j)=sum((Xset(:,cycle)-Xset(:,j)).^2)*log(sqrt(sum((Xset(:,cycle)-Xset(:,j)).^2))+realmin);  
                end                
            end
            polyprime=zeros(cycle,J+1);
            polyprime(1:cycle-1,1:J+1)=poly;
            poly=polyprime;
            for j=1:J
                poly(cycle,j)=Xset(j,cycle);
            end
            poly(cycle,J+1)=1;
        end
 
        F1lambda=[phi poly; poly' zeros(J+1,J+1)]\[F1(1:cycle,1); zeros(J+1,1)];
        F1coeff=F1lambda(cycle+1:cycle+J+1);
        F1lambda=F1lambda(1:cycle);
 
        F2lambda=[phi poly; poly' zeros(J+1,J+1)]\[F2(1:cycle,1); zeros(J+1,1)];
        F2coeff=F2lambda(cycle+1:cycle+J+1);
        F2lambda=F2lambda(1:cycle);
        %then randomly generate candidate points using RBF
        candpt=zeros(J,numcand);
        for i=1:numcand
            samples=zeros(J,1);
            for j=1:J
                samples(j,1)=rand;
            end
            min_cycle=ceil(rand*NumPF);
            candpt(:,i)=PF(min_cycle,4:J+3)'+norminv(samples,zeros(J,1),ro_n);
            candpt(1,i)=min(max(0,candpt(1,i)),1); %in case it goes less than 0 or greater than Ymax
            candpt(2:J,i)=min(max(-5*ones(J-1,1),candpt(2:J,i)),5*ones(J-1,1)); %ZDT4
        end
 
        %Interpolate the candidate pts to obtain Einterpolate and Vinterpolate
        dist=zeros(numcand)+inf;
        F1interpolate=zeros(numcand,1);
        F2interpolate=zeros(numcand,1);
 
        for i=1:numcand
            for j=1:cycle
                F1interpolate(i,1)=F1interpolate(i,1)+F1lambda(j)*sum((candpt(:,i)-Xset(:,j)).^2)*log(sqrt(sum((candpt(:,i)-Xset(:,j)).^2))+realmin);
                F2interpolate(i,1)=F2interpolate(i,1)+F2lambda(j)*sum((candpt(:,i)-Xset(:,j)).^2)*log(sqrt(sum((candpt(:,i)-Xset(:,j)).^2))+realmin);
 
            end
            F1interpolate(i,1)=F1interpolate(i,1)+candpt(:,i)'*F1coeff(1:J)+F1coeff(J+1);
            F2interpolate(i,1)=F2interpolate(i,1)+candpt(:,i)'*F2coeff(1:J)+F2coeff(J+1);
        end
            %%update Pareto set
        CandPF=ones(numcand);
        for i=1:numcand-1
            for k=i+1:numcand
                if CandPF(k)==1
                    if or(and(F1interpolate(i,1)<=F1interpolate(k,1),F2interpolate(i,1)<F2interpolate(k,1)),and(F1interpolate(i,1)<F1interpolate(k,1),F2interpolate(i,1)<=F2interpolate(k,1)))
                        CandPF(k)=0;     %CandPF is the candidate pareto filter
                    elseif or(and(F1interpolate(i,1)>=F1interpolate(k,1),F2interpolate(i,1)>F2interpolate(k,1)),and(F1interpolate(i,1)>F1interpolate(k,1),F2interpolate(i,1)>=F2interpolate(k,1)))
                        CandPF(i)=0;
                    end
                end
            end
        end
        
        %Choose max dist within candidate Pareto frontier to obtain best fit
        for i=1:numcand
            for j=1:NumPF
                if CandPF(i)==1
                    dist(i)=min(dist(i),sqrt((F1interpolate(i)-PF(j,2))^2)+sqrt((F2interpolate(i)-PF(j,3))^2));
                end
            end
        end
        
        %candidate pt  
        score=0;
        min_i=1;
        for i=1:numcand
            if CandPF(i)==1     %this criteria guarantees that only Pareto optimal solutions are considered
                if dist(i)>score
                    score=dist(i);
                    min_i=i;
                end
            end
        end
 
 
        %Evaluate costly function and update information
 
        Xset=[Xset candpt(:,min_i)];
        %ZDT4
        F1(cycle+1,1)=Xset(1,cycle+1);
        gx=1+10*(J-1)+sum(Xset(2:J,cycle+1).^2-10*cos(4*pi*Xset(2:J,cycle+1)));
        F2(cycle+1,1)=gx*(1-sqrt(Xset(1,cycle+1)/gx));
        %%update Pareto set
        PF(NumPF+1,:)=[1 F1(cycle+1,1) F2(cycle+1,1) Xset(:,cycle+1)'];
        stop=0;
        nu=1;
        while and(nu<=NumPF,stop==0)
            if or(and(F1(cycle+1,1)<PF(nu,2),F2(cycle+1,1)<=PF(nu,3)),and(F1(cycle+1,1)<=PF(nu,2),F2(cycle+1,1)<PF(nu,3)))
                 PF(nu,1)=0;
            elseif or(and(F1(cycle+1,1)>=PF(nu,2),F2(cycle+1,1)>PF(nu,3)),and(F1(cycle+1,1)>PF(nu,2),F2(cycle+1,1)>=PF(nu,3)))
                PF(NumPF+1,1)=0;
                stop=1;
            end
            nu=nu+1;
        end
        if PF(NumPF+1,1)==1
            NumPF=size(PF,1);
            oldPF=PF;
            PF=0;
            for i=1:NumPF
                if oldPF(i,1)==1
                    if PF==0
                        PF=oldPF(i,:);
                    else
                        PF=[PF; oldPF(i,:)];
                    end
                end
            end
            NumPF=size(PF,1);            
            Cfail=0;
        else
            PF=PF(1:NumPF,:);   
            Cfail=Cfail+1;
        end
        if Cfail>fmax
            Cfail=0;
            rcount=rcount+1;
            ro_n=0.5*ro_n;
        end
        if cycle==2000
            rcount=rmax+1;
        end
        n=n+1
        cycle=cycle+1
    end
end


