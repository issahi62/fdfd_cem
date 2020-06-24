 %--- main parameters
    rhow = 650;     % density of wood, kg/m^3
    d = 0.02;       % wood particle diameter, m
    Ti = 300;       % initial particle temp, K
    Tinf = 673;     % ambient temp, K
    h = 60;         % heat transfer coefficient, W/m^2*K
    % A = pre-exponential factor, 1/s and E = activation energy, kJ/mol
    A1 = 1.3e8;     E1 = 140;   % wood -> gas
    A2 = 2e8;       E2 = 133;   % wood -> tar
    A3 = 1.08e7;    E3 = 121;   % wood -> char 
    R = 0.008314;   % universal gas constant, kJ/mol*K
    %--- initial calculations
    b = 1;          % shape factor, b = 1 cylinder, b = 2 sphere
    r = d/2;        % particle radius, m
    nt = 1000;      % number of time steps
    tmax = 840;     % max time, s
    dt = tmax/nt;   % time step spacing, delta t
    t = 0:dt:tmax;  % time vector, s
    m = 20;         % number of radius nodes
    steps = m-1;    % number of radius steps
    dr = r/steps;   % radius step spacing, delta r
    %--- build initial vectors for temperature and thermal properties
    i = 1:m;
    T(i,1) = Ti;    % column vector of temperatures
    TT(1,i) = Ti;   % row vector to store temperatures 
    pw(1,i) = rhow; % initial density at each node is wood density, rhow
    pg(1,i) = 0;    % initial density of gas
    pt(1,i) = 0;    % inital density of tar
    pc(1,i) = 0;    % initial density of char
    %--- solve system of equations [A][T]=[C] where T = A\C
    for i = 2:nt+1
        % kinetics at n
        [rww, rwg, rwt, rwc] = funcY(A1,E1,A2,E2,A3,E3,R,T',pw(i-1,:));
        pw(i,:) = pw(i-1,:) + rww.*dt;      % update wood density
        pg(i,:) = pg(i-1,:) + rwg.*dt;      % update gas density
        pt(i,:) = pt(i-1,:) + rwt.*dt;      % update tar density
        pc(i,:) = pc(i-1,:) + rwc.*dt;      % update char density
        Yw = pw(i,:)./(pw(i,:) + pc(i,:));  % wood fraction
        Yc = pc(i,:)./(pw(i,:) + pc(i,:));  % char fraction
        % thermal properties at n
        cpw = 1112.0 + 4.85.*(T'-273.15);   % wood heat capacity, J/(kg*K) 
        kw = 0.13 + (3e-4).*(T'-273.15);    % wood thermal conductivity, W/(m*K)
        cpc = 1003.2 + 2.09.*(T'-273.15);   % char heat capacity, J/(kg*K)
        kc = 0.08 - (1e-4).*(T'-273.15);    % char thermal conductivity, W/(m*K)
        cpbar = Yw.*cpw + Yc.*cpc;  % effective heat capacity
        kbar = Yw.*kw + Yc.*kc;     % effective thermal conductivity
        pbar = pw(i,:) + pc(i,:);   % effective density
        % temperature at n+1
        Tn = funcACbar(pbar,cpbar,kbar,h,Tinf,b,m,dr,dt,T);
        % kinetics at n+1
        [rww, rwg, rwt, rwc] = funcY(A1,E1,A2,E2,A3,E3,R,Tn',pw(i-1,:));
        pw(i,:) = pw(i-1,:) + rww.*dt;
        pg(i,:) = pg(i-1,:) + rwg.*dt;
        pt(i,:) = pt(i-1,:) + rwt.*dt;
        pc(i,:) = pc(i-1,:) + rwc.*dt;
        Yw = pw(i,:)./(pw(i,:) + pc(i,:));
        Yc = pc(i,:)./(pw(i,:) + pc(i,:));
        % thermal properties at n+1
        cpw = 1112.0 + 4.85.*(Tn'-273.15);
        kw = 0.13 + (3e-4).*(Tn'-273.15);
        cpc = 1003.2 + 2.09.*(Tn'-273.15);
        kc = 0.08 - (1e-4).*(Tn'-273.15);
        cpbar = Yw.*cpw + Yc.*cpc;
        kbar = Yw.*kw + Yc.*cpc; 
        pbar = pw(i,:) + pc(i,:);
        % revise temperature at n+1
        Tn = funcACbar(pbar,cpbar,kbar,h,Tinf,b,m,dr,dt,T);
        % store temperature at n+1
        T = Tn;
        TT(i,:) = T';
    end
    %--- plot data
    figure(1)
    plot(t./60,TT(:,1),'-b',t./60,TT(:,m),'-r')
    hold on
    plot([0 tmax/60],[Tinf Tinf],':k')
    hold off
    xlabel('Time (min)'); ylabel('Temperature (K)');
    sh = num2str(h);  snt = num2str(nt);  sm = num2str(m);
    title(['Cylinder Model, d = 20mm, h = ',sh,', nt = ',snt,', m = ',sm])
    legend('Tcenter','Tsurface',['T\infty = ',num2str(Tinf),'K'],'location','southeast')
    
    figure(2)
    plot(t./60,pw(:,1),'--',t./60,pw(:,m),'-','color',[0 0.7 0])
    hold on
    plot(t./60,pg(:,1),'--b',t./60,pg(:,m),'b')
    hold on
    plot(t./60,pt(:,1),'--k',t./60,pt(:,m),'k')
    hold on
    plot(t./60,pc(:,1),'--r',t./60,pc(:,m),'r')
    hold off
    xlabel('Time (min)'); ylabel('Density (kg/m^3)');