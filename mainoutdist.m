%MAINOUTDIST Implements the Self-Triggered Control described in Gleizer and
%Mazo Jr. (2018) for reproducing its results. This code simulates the
%scenario described in Section 5 of Gleizer and Mazo Jr. (2018).
% 
%   Author: Gabriel de A. Gleizer, 2018 (g.gleizer@tudelft.nl)
%
% REFERENCES:
% G. A. Gleizer and M. Mazo Jr. Self-triggered output feedback
%     control for perturbed linear systems. IFAC-PapersOnLine,
%     51(23):248--253, 2018.

%% Flags
DEBUG_PLOTS = false;
INCLUDE_DIST = true;  % Include disturbances

%% Load data
% Plant
Ap = [1.38, -0.208, 6.715, -5.676;...
    -0.581,  -4.29,  0,     0.675;...
     1.067,  4.273, -6.654, 5.893;...
     0.048,  4.273, 1.343, -2.104];
 
Bp = [0,   0;...
    5.679, 0;...
    1.136, -3.146;...
    1.136, 0];

E = [1; 0; 0; 0];

Cp = [eye(2), [1 -1; 0 0]];

% Controller
%K = [1 -4];
h = 0.01;

Ac = [1 0; 0 1];
Bc = h*[0 1; 1 0];
Cc = [-2 0; 0 8];
Dc = [0 -2; 5 0];

rho = 0.6;
gamma = 0.5;

% Triggering conditions
sigma = 0.1;  %0.25;

% Bound on disturbance
W = 0.1;

% Dimensions
np = size(Ap,1);  % states of the plant
nc = size(Ac,1);  % states of the controller
pp = size(Cp,1);  % measured plant outputs
mp = size(Bp,2);  % number of control inputs
nw = size(E,2);   % number of disturbances

%% Now put the system in the normal observable form
T = [Cp; Cp*Ap];
Tinv = inv(T);
Ap = T*Ap*Tinv;
Bp = T*Bp;
E = T*E;
Cp = Cp*Tinv;

%% Simulation data
x0 = [1; -1; -1; 1]*10;
x0 = T*x0;
xc0 = [0; 0];
USE_INPUT = true;
y0 = x0(1:mp);

TEND = 10;
kfinal = 20;  % For plots (avoids anoying outlier)

% Disturbance signal
omega = @(t) W*sin(pi*t).*((t >= 0) & (t <= 8));
omega = @(t) W*((t >= 0) & (t <= TEND/2));
odeplant = @(t,xp,u) Ap*xp + Bp*u + E*omega(t);

%% Q Matrix
Q1 = [(1 - sigma^2)*(Cp'*Cp), zeros(np,nc);...
    zeros(nc,np), (1 - sigma^2)*(Cc'*Cc)];
Q2 = [-Cp', zeros(np,mp);...
    (1-sigma^2)*Cc'*Dc, -Cc'];
Q3 = [eye(pp) + (1 - sigma^2)*(Dc'*Dc), -Dc'; ...
    -Dc, eye(mp)];

Q = [Q1, Q2; Q2', Q3];

if ~USE_INPUT
    YtoXI = [Cp, zeros(pp,nc+pp+mp)];
    YHATtoXI = [zeros(pp,np+nc),eye(pp),zeros(pp,mp)];
    EtoXI = YtoXI - YHATtoXI;
    Q = EtoXI'*EtoXI - sigma^2*(YtoXI'*YtoXI);
    Q1 = Q(1:(np+nc),1:(np+nc));
    Q2 = Q(1:(np+nc),(np+nc+1):end);
    Q3 = Q((np+nc+1):end,(np+nc+1):end);
end

% NEW: Ctheta
Ctheta = [eye(np); zeros(nc,np)];

%% Discrete behavior
CE = [Cp, zeros(pp,nc); Dc*Cp, Cc];
Phip = expm(Ap*h);
Mfixed = (Ap*[eye(np), zeros(np,nc)] + ...
    Bp*[Dc*Cp, Cc]);
I0 = [eye(np),zeros(np,nc)];
OI = [zeros(nc,np), eye(nc)];

% Loop to compute Mks
Phipk = Phip;
Ack = Ac;
aggr = Bc*[Cp, zeros(pp,nc)];
KMAX = 1000;

try
    clear M1 M2 M;
end
for k = 1:KMAX
    M1{k} = I0 +  Ap\(Phipk - eye(np))*Mfixed;
    M2{k} = Ack*OI + aggr;
    M{k} = [M1{k}; M2{k}];
    Phipk = Phip*Phipk;
    aggr = aggr + Ack*Bc*[Cp, zeros(pp,nc)];
    Ack = Ac*Ack;
end

%% Conic equation matrices
try
    clear Qq F FF
end

for k = 1:KMAX
   Qq{k} = [M{k}', CE']*Q*[M{k}; CE];
   lambda = eig(Qq{k});
   maxeig(k) = max(real(lambda));
   mineig(k) = min(real(lambda));
   if mineig(k) > 0
       break;
   end
end

kbeg = find(maxeig > 0, 1, 'first');
kend = find(mineig > 0, 1, 'first');
if isempty(kend)
    kend = KMAX;
end

dkmax = min(kend,KMAX);

%% Reshape from cell into array
QQ = reshape(cell2mat(Qq),np+nc,np+nc,dkmax);

%% Aproximation of the effect of bounded disturbance
% Bound on thetak
a = max(eig(Ap'+Ap))/2;
enorm = norm(E);
kv = 1:kend;
theta = abs((1-exp(a*kv*h))/a)*enorm*W;
% New: using Schur based norm - need to run expmattest.m first...
for k = kv
    %theta(k) = intschur(h*k)*enorm*W;
    F{k} = 2*Ctheta'*[Q1, Q2]*[M{k}; CE];
    c(k) = theta(k)^2*norm(Ctheta'*Q1*Ctheta);
end

%% Extract data for simulink
FF = reshape(cell2mat(F),np,np+nc,kend);

%% Try the one using some information on the norm of the state.

% factor for how bigger initial estimate of the unmeasured states is
CF = 1.5;  

%Initialization
x = [x0; xc0];
xnorm = norm([x0;xc0]);
H = eye(np-mp);
ybarc = zeros(np-mp,1);
ybarnorm = sqrt(xnorm^2 - norm([y0;xc0])^2);
ybarnorm = CF*ybarnorm;

y = y0;
xc = xc0;
k = 0;

% Update inv(H) instead of H
Hi = H;
Hi = Hi*ybarnorm^2;
ybarnorm = 1;

% Save variables
xlog = [];
klog = [];
dklog = [];
kslog = [];
Hilog = Hi;
yclog = ybarc';

% Run the closed loop system
if DEBUG_PLOTS
    figure(10);
    hold all;
end

while k <= TEND/h
    xlog = [xlog; x'];
    klog = [klog; k];
    y = Cp*x(1:np);
    xc = x(np+1:end); 
    
    xhat = [y; ybarc; xc];
    
    for dk = kbeg:kfinal
        % Q indexing
        Qk = QQ(:,:,dk);
        Qkj = Qk(:,mp+1:np);
        Qkjj = Qk(mp+1:np,mp+1:np);
        Qkjn = Qk(mp+1:np,1:np);
        
        % Worst case computation
        xQx = xhat'*Qk*xhat;
        Qx = Qkj'*xhat;
        
        xQHQx = Qx'*Hi*Qx;
        lH = eigs(Hi*Qkjj,1,'largestreal');
        
        xQe = 2*sqrt(ybarnorm*xQHQx);
        eQe = lH*ybarnorm^2;
        
        if INCLUDE_DIST
            xQw = norm(FF(:,:,dk)*xhat)*theta(dk);
            wQw = c(dk);
            eQw = 2*theta(dk)*sqrt(eigs(Qkjn'*Hi*Qkjn,1));
        else
            xQw = 0;
            wQw = 0;
            eQw = 0;
        end
        
        if 0
            % Alternative way: consolidate the uncertain factor
            % Note: it sucks and there may be a reason why
            Hic = zeros(np);
            Hic(mp+1:end,mp+1:end) = Hi;
            Hiw = theta(dk)^2*eye(np);
            pstar = sqrt(trace(Hic))/sqrt(trace(Hiw));
            if pstar == 0
                Hic = Hiw;
            elseif ~isinf(pstar)
                Hic = (1 + 1/pstar)*Hic + (1 + pstar)*Hiw;
            end

            Fx = FF(:,:,dk)*xhat;
            xFHFx = Fx'*Hic*Fx;
            llH = eigs(Hic*Qk(1:np,1:np),1,'largestreal');
            xQr = 2*sqrt(xFHFx);
            rQr = llH;
        end
        
        if xQx + xQe + eQe + xQw + wQw + eQw > 0
            dklog = [dklog; dk];
            kslog = [kslog; k];
            break;
        end
    end
    % Get transition matrix
    Mk = M{dk};
    
    % Retrieve controller data
    y = Cp*x(1:np);
    xc = x(np+1:end);
    
    % Walk system
    Mk1 = Mk([1:mp end-nc+1:end], [1:mp end-nc+1:end]);
    Mk2 = Mk([1:mp end-nc+1:end], [mp+1:np]);
    Mk3 = Mk([mp+1:np], [1:mp end-nc+1:end]);
    Mk4 = Mk([mp+1:np], [mp+1:np]);   
    
    % Considering disturbances, we need a proper solver
    nuhat = Cc*xc + Dc*y;
    [~,xpode] = ode45(@(t,x)odeplant(t,x(1:np),nuhat), h*[k, k+dk], x(1:np));
    xpnext = xpode(end,:)';
    
    % Get next values (x, xc, y)
    xnext = Mk*x;
    xnext(1:np) = xpnext;
    ynext = Cp*xnext(1:np);
    xcnext = xnext(np+1:end);
    
    % Compute observer gain
    L = zeros(np-mp,mp+nc);
    
    % Update H and ybarc
    a = (Mk3 + L*Mk1)*[y; xc] - L*[ynext; xcnext];
    %Mi = inv(Mk4);
    %H = Mi'*H*Mi;
    Hi = (Mk4 + L*Mk2)*Hi*(Mk4 + L*Mk2)';
    %H = (H+H')/2;
    Hi = (Hi+Hi')/2;

    ybarc = a + (Mk4 + L*Mk2)*ybarc;
    yclog = [yclog; ybarc'];
        
    if INCLUDE_DIST
        if 0
            % New: using the ellipsoidal toolbox
            Ed = ellipsoid(Hi);
            Ew = theta(dk)*ell_unitball(np-mp);
            % Get direction corresponding to the largest eigenvalue
            [vv,dd] = eigs(Hi,1);
            Em = minksum_ea([Ed,Ew],vv);
            [~, Hi] = parameters(Em);
            plot(Em)
        else
            % Alternative: min sum of squares of the semiaxes (chap 2.5
            % of the book Ellipsoidal Calculus for Estimation and
            % Control)
            Hball = theta(dk)^2*eye(np-mp);
            pstar = sqrt(trace(Hi))/sqrt(trace(Hball));
            if pstar == 0
                Hi = Hball;
            elseif ~isinf(pstar)
                Hi = (1 + 1/pstar)*Hi + (1 + pstar)*Hball;
            end
            if DEBUG_PLOTS
                plot(ellipsoid(ybarc,Hi));
            end
        end
    end
    Hilog(:,:,end+1) = Hi;

%     % Rescale H and ybarnorm because H tends to increase
%     lH = norm(H);
%     H = H/lH;
%     ybarnorm = ybarnorm/lH
%     
    % Next step
    x = xnext;
    k = k+dk;
end

%% Simulate with complete state information

%Initialization
x = [x0; xc0];
H = eye(np-mp);
ybarc = x(mp+1:np);

y = y0;
xc = xc0;
k = 0;

% Update inv(H) instead of H
Hi = H*0;
ybarnorm = 1;

% Save variables
xlognd = [];
klognd = [];
kslognd = [];
dklognd = [];

while k <= TEND/h
    xlognd = [xlognd; x'];
    klognd = [klognd; k];
    
    xhat = x;
    
    for dk = kbeg:kfinal
        % Q indexing
        Qk = QQ(:,:,dk);
        Qkj = Qk(:,mp+1:np);
        Qkjj = Qk(mp+1:np,mp+1:np);
        Qkjn = Qk(mp+1:np,1:np);
        
        % Worst case computation
        xQx = xhat'*Qk*xhat;
        Qx = Qkj'*xhat;
        
        xQHQx = Qx'*Hi*Qx;
        lH = eigs(Hi*Qkjj,1,'largestreal');
        
        xQe = 2*sqrt(ybarnorm*xQHQx);
        eQe = lH*ybarnorm^2;
        
        if INCLUDE_DIST
            xQw = norm(FF(:,:,dk)*xhat)*theta(dk);
            wQw = c(dk);
            eQw = 2*theta(dk)*sqrt(eigs(Qkjn'*Hi*Qkjn,1));
        else
            xQw = 0;
            wQw = 0;
            eQw = 0;
        end
        
        if xQx + xQe + eQe + xQw + wQw + eQw > 0
            dklognd = [dklognd; dk];
            kslognd = [kslognd; k];
            break;
        end
    end
    % Get transition matrix
    Mk = M{dk};
    
    % Retrieve controller data
    y = Cp*x(1:np);
    xc = x(np+1:end);
    
    % Considering disturbances, we need a proper solver
    nuhat = Cc*xc + Dc*y;
    [~,xpode] = ode45(@(t,x)odeplant(t,x(1:np),nuhat), h*[k, k+dk], x(1:np));
    xpnext = xpode(end,:)';
    
    % Get next values (x, xc, y)
    xnext = Mk*x;
    xnext(1:np) = xpnext;
    ynext = Cp*xnext(1:np);
    xcnext = xnext(np+1:end);
    
    x = xnext;
    k = k+dk;
end

%% Check the PETC
% Discrete model
s = [Ap, Bp; zeros(pp,np+pp)];
sd = expm(s*h);
Phip = sd(1:np,1:np);
Gammap = sd(1:np,np+1:end);

%Initialization
x = [x0; xc0];
y = y0;
xc = xc0;
k = 0;

xhat = x;
dk = 0;

% For the Ellipsoid-observer based prediction
xnorm = norm([x0;xc0]);
H = eye(np-mp);
ybarc = zeros(np-mp,1);
ybarnorm = sqrt(xnorm^2 - norm([y0;xc0])^2);
ybarnorm = CF*ybarnorm;
Hi = H;
Hi = Hi*ybarnorm^2;
ybarnorm = 1;

% Initialize hat variables
yhat = xhat(1:mp);
uhat = Cc*xc + Dc*yhat;
xchat = xc;

% Save variables
xloge = [];
kloge = [];
dkloge = [];
ksloge = [];
dklogs = dklog(1);  % Store STC-based worst case scenario

while k < TEND/h
    % Log x and k
    kloge = [kloge;k];
    xloge = [xloge;x'];
    
    % Walk system
    xp = x(1:np);
    y = xp(1:mp);
    xc = x(np+1:end);
    
    % Considering disturbances, we need a proper solver
    nuhat = Cc*xc + Dc*y;
    [~,xpode] = ode45(@(t,x)odeplant(t,x(1:np),nuhat), h*[k, k+1], x(1:np));
    xpnext = xpode(end,:)';
    xcnext = Ac*xc + Bc*yhat;
    xnext = [xpnext; xcnext];
    x = xnext;
    k = k+1;
    dk = dk+1;
    
    % Get current data
    xp = x(1:np);
    y = xp(1:mp);
    xc = x(np+1:end);
    
    %Qk = QQ(:,:,dk);
    
    u = Cc*xc + Dc*yhat;
    yerr = y - yhat;
    uerr = u - uhat;
    
    if k*h > 3
        macaco = 1;
    end
    if norm([yerr; uerr]) > sigma*norm([y;u]) || (dk == kfinal)  %xhat'*Qk*xhat > 0
        % LB: Update ellipsoid
        Mk = M{dk};
        
        Mk3 = Mk(mp+1:np, [1:mp end-nc+1:end]);
        Mk4 = Mk(mp+1:np, mp+1:np);
        ys = yhat;
        xcs = xchat;
        % Update H and ybarc
        a = Mk3*[ys; xcs];
        Hi = Mk4*Hi*Mk4';
        Hi = (Hi+Hi')/2;

        Hball = theta(dk)^2*eye(np-mp);
        pstar = sqrt(trace(Hi))/sqrt(trace(Hball));
        if pstar == 0
            Hi = Hball;
        elseif ~isinf(pstar)
            Hi = (1 + 1/pstar)*Hi + (1 + pstar)*Hball;
        end

        ybarc = a + Mk4*ybarc;
    
        
        % ETC: Update hat variables 
        xhat = x;
        yhat = xhat(1:mp);
        uhat = Cc*xc + Dc*yhat;
        dkloge = [dkloge; dk];
        ksloge = [ksloge; k];
        xchat = xc;  % For ellipsoid - store xc at last triggering
        
        % LB: estimate next triggering time
        xhats = [yhat; ybarc; xc];
        
        for dks = kbeg:kfinal
            % Q indexing
            Qk = QQ(:,:,dks);
            Qkj = Qk(:,mp+1:np);
            Qkjj = Qk(mp+1:np,mp+1:np);
            Qkjn = Qk(mp+1:np,1:np);
                
            % Worst case computation
            xQx = xhats'*Qk*xhats;
            Qx = Qkj'*xhats;
            
            xQHQx = Qx'*Hi*Qx;
            xQe = 2*sqrt(xQHQx);
            eQe = eigs(Hi*Qkjj,1,'largestreal');
            
            if INCLUDE_DIST
                xQw = norm(FF(:,:,dks)*xhat)*theta(dks);
                wQw = c(dks);
                eQw = 2*theta(dks)*sqrt(eigs(Qkjn'*Hi*Qkjn,1));
            else
                xQw = 0;
                wQw = 0;
                eQw = 0;
            end

            if xQx + xQe + eQe + xQw + wQw + eQw > 0
                dklogs = [dklogs; dks];
                break;
            end
        end
        dk = 0;
    end
end

%% Debug plots...
if DEBUG_PLOTS
    figure(1);
    plot(kslog,dklog);
    hold all;
    plot(kslognd, dklognd);
    legend('STC only output','STC full state');

    figure(2);
    plot(klog,xlog);

    figure(3);
    plot(ksloge,dkloge,'x-')
    hold all;
    plot(ksloge, dklogs(1:end-1),'o-')

    figure(4);
    plot(kloge,xloge);
end
