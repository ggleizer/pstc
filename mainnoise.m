%MAINNOISE Implements the Preventive Self-Triggered Control described in
%Gleizer and Mazo Jr. (2019) for reproducing the results of this paper.
%This code simulates the scenarios described in Section 5.
%
%   Requirements: Optimization Toolbox, Ellipsoidal Toolbox 
%   The Ellipsoidal Toolbox can be found in
%       https://nl.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
% 
%   Author: Gabriel de A. Gleizer, 2018--2019 (g.gleizer@tudelft.nl)
%
% REFERENCES:
% G. A. Gleizer and M. Mazo Jr. Self-triggered output-feedback control of
%     LIT systems subject to disturbances and noise. Submitted to
%     Automatica, 2019.
% A. B. Kurzhanski and I. Valyi. Ellipsoidal toolbox, 2006.
%     Technical Report.
% F. Schweppe. Recursive state estimation: Unknown but bounded errors and
%     system inputs. IEEE Transactions on Automatic Control, 13(1):22--28,
%     1968.

%% Flags
DEBUG = 0;
TEST_PLOTS = 0;
if exist('TRIG_LEVEL_GLOBAL','var')
    TRIG_LEVEL = TRIG_LEVEL_GLOBAL;
    disp('Using externally defined triggering "eps" level (TRIG_LEVEL_GLOBAL)');
else
    TRIG_LEVEL = 0.00;
end

fileId = 1;
if exist('FILE_ID_GLOBAL','var')
    fileId = FILE_ID_GLOBAL;
end


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

% Controller update rate / fundamental STC step size
h = 0.01;

% Dynamic controller data
Ac = [1 0; 0 1];
Bc = h*[0 1; 1 0];
Cc = [-2 0; 0 8];
Dc = [0 -2; 5 0];

% Triggering condition
sigma = 0.1;

% Dimensions
np = size(Ap,1);  % states of the plant
nc = size(Ac,1);  % states of the controller
pp = size(Cp,1);  % measured plant outputs
mp = size(Bp,2);  % number of control inputs
nw = size(E,2);   % number of disturbances

assert(rank([Cp; Cp*Ap; Cp*Ap^2; Cp*Ap^3]) == 4,...
       'System is not observable');


%% Assumption 4: Bound on disturbance
W_MAG = 0.1;

% Bound on noise
if exist('V_GLOBAL','var')
    V_EACH_ELEMENT = V_GLOBAL;
    disp('Using externally set noise value (V_GLOBAL)');
else
    V_EACH_ELEMENT = 0.01;
end

YFACTOR = 1.1;  % We never now the noise levels that precisely.
V = V_EACH_ELEMENT^2*YFACTOR^2*eye(pp)*pp;


%% If there is no noise, put the system in the normal observable form
if all(V(:) == 0)
    T = [Cp; Cp*Ap];
    Tinv = inv(T);
    Ap = T*Ap*Tinv;
    Bp = T*Bp;
    E = T*E;
    Cp = Cp*Tinv;
else
    T = eye(np);
end


%% Simulation data
x0 = [1; -1; -1; 1];
x0 = 10*x0;  % To compare with the NecSys paper...
x0 = T*x0;
xc0 = [0; 0];
USE_INPUT = true;
y0 = Cp*x0(1:np);

TEND = 10;
kfinal = 25;  % \bar{\kappa} in the paper

% Disturbance signal
%omega = @(t) W_MAG*sin(pi*t).*((t >= 0) & (t <= 8));
%omega = @(t) W_MAG*cos(pi*t).*((t >= 0) & (t <= 8));
omega = @(t) W_MAG*((t >= 0) & (t <= TEND/2));  % Like in NecSys
odeplant = @(t,xp,u) Ap*xp + Bp*u + E*omega(t);

% For reproducibility, pre-compute the noises
rng(1907);
noises = 2*V_EACH_ELEMENT*(rand(pp,TEND/h + kfinal + 1) - 0.5);


%% Timing
timePreProcess = 0;
timeWk = 0;
timeOfflineMatrices = 0;
timeInit = 0;
timesSTC = [];
timesFusion = [];
timesEta = [];
timesPrediction = [];

ticOfflineBegins = tic;


%% Ellipsoidal reachability for computation of \mathcal{X}_w
tic;

Wel = W_MAG*ell_unitball(nw);  % Ellipsoidal Toolbox command for a ball.
sys = linsys(Ap,E,Wel);  % Ellipsoidal Toolbox command
X0 = 1e-4*ell_unitball(np);  % Initial state set, should be the origin.
TINTV = [0,kfinal*h];  % Time interval to compute reachability
L0 = eye(np);  % Support vectors for tight approximation

reachoptions.approximation = 0;  % external (outer) approximation
RS = reach(sys, X0, L0, TINTV, reachoptions);  % Ellipsoidal Toolbox

% Extract matrices W_\kappa from the Ellipsoidal Toolbox RS structure
Wk = zeros(np,np,kfinal);
for kk = 1:kfinal
    RC = cut(RS,kk*h);  % Gets slice of the tube at instant kk*h
    EAC = get_ea(RC);  % Get the array of ellipsoids
    % Iteration to compute the intersection among EAC(i)
    EI = EAC(1);
    for ii = 2:size(EAC,1)
        EI = intersection_ea(EI,EAC(ii));
    end
    Wk(:,:,kk) = parameters(EI);  % Extract W_\kappa
end

timeWk = toc;


%% Q Matrix (\bar{Q} in the paper)
nz = pp+mp;  % z / zeta variable

Q1 = (1-sigma^2)*eye(nz);
Q2 = -eye(nz);
Q3 = eye(nz);
Q = [Q1, Q2; Q2', Q3];


%% Transition matrices
% Compute transition matrices:
%   M(\kappa) such that \xi_p(k+\kappa) = M(\kappa)[xp;xc;y]
%   N(\kappa) such that \zeta(k+\kappa) = M(\kappa)[xp;xc;y]

% First the more obvious CE: [y;u] = CE[xp;xc;y]
CE = [zeros(nz,np), [zeros(pp,nc), eye(pp); Cc, Dc]];

% Fundamental transition matrices
Abar = [Ap, Bp; zeros(mp,np+mp)];
Phibar = expm(Abar*h);
Phip = Phibar(1:np,1:np);
Gammap = Phibar(1:np,np+1:end);

% Auxiliary matrices
I0 = [eye(np),zeros(np,nc)];
OI = [zeros(nc,np), eye(nc)];

% Loop to compute Mks
Phipk = Phip;
Gammapk = Gammap;
Ack = Ac;
Bck = Bc;
KMAX = kfinal;

try
    clear M N;
end
for k = 1:KMAX
    M1 = [Phipk, Gammapk*Cc, Gammapk*Dc];
    M2 = [zeros(nc,np), Ack, Bck];
    N1 = Cp*[Phipk, Gammapk*Cc, Gammapk*Dc];
    N2 = [zeros(mp,np), Cc*Ack, Cc*Bck + Dc];
    M{k} = [M1; M2];
    N{k} = [N1; N2];

    Phipk = Phip*Phipk;
    Ack = Ac*Ack;
    
    Gammapk = Gammap + Phip*Gammapk;
    Bck = Bc + Ac*Bck;
end

%% Conic equation matrices (Q_\kappa)
try
    clear Qq
end

for k = 1:KMAX
   Qq{k} = [N{k}', CE']*Q*[N{k}; CE];
   lambda = eig(Qq{k});
   maxeig(k) = max(real(lambda));
   mineig(k) = min(real(lambda));
   % If Qk > 0, every state will have triggered up to here. Can stop.
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

% Turn cell array into a 3-array
nh = np+nc+pp;
QQ = reshape(cell2mat(Qq),nh,nh,dkmax);

timeOfflineMatrices = toc(ticOfflineBegins);

%% Get auxiliary matrices and offline computable vectors
% Eq. (21) in the paper

tic;

Cw = [Cp; zeros(nz+mp,np)];
Cv = [eye(pp); zeros(nz+mp,pp)];

% Pre-allocate kappa-dependent variables
Fv = zeros(nh,pp,kfinal);
Fw = zeros(nh,np,kfinal);
cvw = zeros(kfinal,1);
wQw = zeros(kfinal,1);
Rw = zeros(nh,nh,kfinal);  % Fw W Fw'
Rv = zeros(nh,nh,kfinal);  % Fv V Fv'

% Compute kappa-independent variables
Qv = Cv'*Q*Cv;
Qw = Cw'*Q*Cw;
cv = eigs(V*Qv,1);

% Auxiliary matrix
Qwv = Cw'*Q*Cv;

% Compute kappa-dependent variables
for kk = 1:kfinal
    Fv(:,:,kk) = [N{kk}; CE]'*Q*Cv;
    Fw(:,:,kk) = [N{kk}; CE]'*Q*Cw;  
    Rw(:,:,kk) = Fw(:,:,kk)*Wk(:,:,kk)*Fw(:,:,kk)';
    Rv(:,:,kk) = Fv(:,:,kk)*V*Fv(:,:,kk)'; 
    wQw(kk) = eigs(Wk(:,:,kk)*Qw,1);
    cvw(kk) = sqrt(eigs(V*(Qwv'*Wk(:,:,kk)*Qwv), 1));
end

timeOfflineMatrices = timeOfflineMatrices + toc;
timePreProcess = toc(ticOfflineBegins);

%%
%%%%%%%%%%%%%%%%%%%% SIMULATION BEGINS  %%%%%%%%%%%%%%%%%%%%%%%

% Initialization, log variables
x = [x0; xc0];
k = 0;

% Log variables
xlog = [];
klog = [];
dklog = [];
kslog = [];
Xlog = inf(np,np);
xtlog = [];
ylog = [];
dkzlog = [];

%% Initialization phase (Appendix C)

% Initialize intersection of ellyptical cylinders
% Determine kbar
Obs = [];
kbar = 0;
while true
    Obs = [Cp; Obs];
    if rank(Obs) == np
        break;
    end
    Obs = Obs*Phip;
    kbar = kbar + 1;
end

Obsbar = Obs*Phip^(-kbar);
    
psibar = [];
Vtilde = [];

while k <= TEND/h
    xlog = [xlog; x'];
    ylog = [ylog; (Cp*x(1:np))'];
    klog = [klog; k];
    xtlog = [xtlog; inf(1,np)];
    
    xhat = x;
    y = Cp*xhat(1:np);
    
    % Add noise
    y = y + noises(:, k + 1);
    xy = [x;y];
    
    % Collect outputs
    psibar = [psibar; y];
    
    % Stopping criteria: Obs is full rank at this point
    if k == kbar
        break
    end
    
    % Simulation goes on: during initialization, always trigger with dk = 1
    dk = 1;
    dklog = [dklog; dk];
    kslog = [kslog; k];
    
    % Get transition matrix
    Mk = M{dk};
    Mknn = Mk(1:np,1:np);
    Mkncy = Mk(1:np,np+1:end);
    Wkk = Wk(:,:,dk);
    
    % Retrieve controller data
    xc = x(np+1:end);
    u = Cc*xc + Dc*y;
    z = [y;u];
    
    % Run plant forward
    [~,xpode] = ode45(@(t,x)odeplant(t,x(1:np),u), h*(k:k+kfinal), x(1:np));
    xpnext = xpode(dk+1,:)';

    % Check when PETC would have triggered
    yz = xpode*Cp';
    yz = yz + noises(:,(k:k+kfinal) + 1)';

    uz = zeros(kfinal+1,mp);
    uz(1,:) = u';
    xcz = xc;
    for dkz = 2:kfinal+1
        xcz = Ac*xcz + Bc*y;
        uz(dkz,:) = (Cc*xcz + Dc*y)';
    end

    zetaz = [yz, uz];
    ezz = zetaz - repmat(z',kfinal+1,1);
    etaz = sum(ezz.^2,2) - sigma^2*sum(zetaz.^2,2);
    dkz = max(1, find(etaz > TRIG_LEVEL, 1, 'first') - 1);
    if isempty(dkz)
        dkz = kfinal;
    end
    assert(dkz >= dk, 'ETC would have triggered before!');
    dkzlog = [dkzlog, dkz];
    % End PETC calculations
    
    % Compute values for k+1 (x, xc, y)
    xnext = Mk*xy;
    xnext(1:np) = xpnext;
    ynext = Cp*xnext(1:np);
    xcnext = xnext(np+1:end);
    
    % Update psitilde values to current k (\tilde{\psi}(0,k))
    for ii = 0:k
        psibar(ii*pp+1:(ii+1)*pp) = psibar(ii*pp+1:(ii+1)*pp) +...
            Cp*Phip^(ii-k-1)*Gammap*u;
    end

    % Log state observer shape matrix just for time consistency
    Xlog(:,:,end+1) = inf(np,np);

    x = xnext;
    k = k+dk;
end

% Data gathered for initialization. Compute ellipsoid now
% -- this is an offline phase, could have been done above --
% 1. Build Vtilde(k) (Lemma 4 + minkowski sum)
tic;
for ii = 0:kbar-1  % No disturbance for kbar
    Vtildek = V;
    CC = Cp*Phip^(kbar-ii);
    Wkk = Wk(:,:,kbar-ii);
    CW = CC*Wkk*CC';
    
    % Minkowski sum outer-approximation
    pstar = sqrt(trace(Vtildek))/sqrt(trace(CW));
    if pstar == 0
        Vtildek = CW;
    elseif ~isinf(pstar)
        Vtildek = (1 + 1/pstar)*Vtildek + (1 + pstar)*CW;
    end

    Vtildek = (Vtildek+Vtildek')/2;
    Vtilde(ii*pp+1:(ii+1)*pp,:) = Vtildek;
end
Vtilde(kbar*pp+1:(kbar+1)*pp,:) = V;

thisTime = toc;
timeOfflineMatrices = timeOfflineMatrices + thisTime;
timePreProcess = timePreProcess + thisTime;
tic;

%2. Compute Vbar (Theorem 4)
Vbar = zeros(pp*(kbar+1));
for ii = 0:kbar
    Vbar(ii*pp+1:(ii+1)*pp,ii*pp+1:(ii+1)*pp)...
        = (kbar+1)*Vtilde(ii*pp+1:(ii+1)*pp,:);
end

% 3. Affine transformation: pinv(Obsbar)*(Obsbar*x)
xptilde = Obsbar\psibar;
X = Obsbar\Vbar/Obsbar';
X = (X+X')/2;
timeInit = toc;

% Save first state observer ellipsoid
E0 = ellipsoid(xptilde,X);

%% PSTC simulation
if TEST_PLOTS
    figure;
    hold all;
end

% Current plant output
y = Cp*x(1:np) + noises(:,k + 1);

timeThisStepSTC = timeInit; 

% Main simulation loop
while k <= TEND/h
    % Update log variables
    xlog = [xlog; x'];
    klog = [klog; k];
    ylog = [ylog; y'];
    xtlog = [xtlog; xptilde'];
    
    if all(V(:)==0) && TEST_PLOTS
        % Test plots for noiseless scenario (only unmeasured states matter)
        plot(ellipsoid(xptilde(3:4),(X(3:4,3:4)+X(3:4,3:4)')/2));
        plot(x(3),x(4),'b*');
    end
    
    % Current states
    xhat = x;
    xc = xhat(np+1:end);
    
    % PSTC triggering mechanism
    ticThisStepSTC = tic;
    
    ptilde = [xptilde; xc; y];  % Algorithm 1, line 3
    
    for dk = kbeg:kfinal  % Algorithm 1, lines 4--10
        tic;
        
        % Matrices indexing
        Qk = QQ(:,:,dk);
        Qkn = Qk(:,1:np);
        Qknn = Qk(1:np,1:np);
        
        Wdk = Wk(:,:,dk);
        Fxwk = Fw(1:np,:,dk);
        Fxvk = Fv(1:np,:,dk);
        Rwdk = Rw(:,:,dk);
        Rvdk = Rv(:,:,dk);
        
        % Compute etabar (Theorem 3)
        xQx = ptilde'*Qk*ptilde;
        Qp = Qkn'*ptilde;
        xQe = sqrt(Qp'*X*Qp);
        eQe = eigs(X*Qknn,1,'largestreal');
        xQw = sqrt(ptilde'*Rwdk*ptilde);
        eQw = sqrt(eigs(Rwdk(1:np,1:np)*X,1));
        xQv = sqrt(ptilde'*Rvdk*ptilde);
        eQv = sqrt(eigs(Rvdk(1:np,1:np)*X,1));
        
        etabar = xQx + 2*xQe + eQe + 2*xQw + wQw(dk) + 2*eQw + cv...
            + 2*cvw(dk) + 2*xQv + 2*eQv;
        
        % Time
        timesEta(end+1) = toc;
        
        if DEBUG  % Test what would be the triggering function using ETC
            udeb = Cc*xc + Dc*y; 
            
            % Run plant forward
            [~,xpodedeb] = ode45(@(t,x)odeplant(t,x(1:np),udeb),...
                                 h*[k, k+dk], x(1:np));
            xpnextdeb = xpodedeb(end,:)';
            
            % Get next values (x, xc, y)
            xnextdeb = M{dk}*[x;y];
            xnextdeb(1:np) = xpnextdeb;
            ynextdeb = Cp*xnextdeb(1:np);
            xcnextdeb = xnextdeb(np+1:end);
            
            % Add noise
            ynextdeb = ynextdeb + noises(:,k + dk + 1);
            
            z = [y;udeb];
            zeta = [ynextdeb; Cc*xcnextdeb + Dc*y]; 
            actualeta = norm(zeta-z)^2 - (sigma*norm(zeta))^2;
            assert(etabar > actualeta,...
                   'STC''s etabar not bigger than PETC''s eta');
        end

        if etabar > TRIG_LEVEL  % epsilon^2 in the paper
            dklog = [dklog; dk];
            kslog = [kslog; k];
            break;
        end
    end
    
    % Compute \mathcal{X}(k+dk|k)
    tic;
    
    % Get transition matrix
    Mk = M{dk};
    Wkk = Wk(:,:,dk);
    Mkn = Mk(1:np,:);  % the same as [Phip(dk), Gamma(dk)]
    Mknn = Mk(1:np,1:np); % the same as Phip(dk)
    
    % Algorithm 1, line 11
    xptilde = Mkn*ptilde;  
    X = Mknn*X*Mknn';
    
    % Algorithm 1, line 12
    pstar = sqrt(trace(X))/sqrt(trace(Wkk));
    if pstar == 0
        X = Wkk;
    elseif ~isinf(pstar)
        X = (1 + 1/pstar)*X + (1 + pstar)*Wkk;
    end
    X = (X+X')/2;

    timesPrediction(end+1) = toc;
    
    % STC of this step is done: final time
    timeThisStepSTC = timeThisStepSTC + toc(ticThisStepSTC);
    timesSTC(end+1) = timeThisStepSTC;
    
    % Run plant forward
    u = Cc*xc + Dc*y;
    z = [y; u];
    [~,xpode] = ode45(@(t,x)odeplant(t,x(1:np),u),...
                      h*(k:k+kfinal), x(1:np));
    xpnext = xpode(dk+1,:)'; 
   
    % Check when PETC would have triggered
    yz = xpode*Cp';
    yz = yz + noises(:, (k:k+kfinal) + 1)';

    uz = zeros(kfinal+1,mp);
    xcz = xc;
    uz(1,:) = u';
    for dkz = 2:kfinal+1
        xcz = Ac*xcz + Bc*y;
        uz(dkz,:) = (Cc*xcz + Dc*y)';
    end

    zetaz = [yz, uz];
    ezz = zetaz - repmat(z',kfinal+1,1);
    etaz = sum(ezz.^2,2) - sigma^2*sum(zetaz.^2,2);
    dkz = max(1, find(etaz > TRIG_LEVEL, 1, 'first') - 1);
    if isempty(dkz)
        dkz = kfinal;
    end
    assert(dkz >= dk, 'ETC would have triggered before!');
    dkzlog = [dkzlog, dkz];
    % End PETC calculations
    
    % ASSERT: did prediction work?
    assert((xptilde - xpnext)'*(X\(xptilde - xpnext)) <= 1,...
           'Prediction failed');
    
    % Get next values (x, xc, y)
    % ATTENTION: this is the next step, so PSTC actually starts HERE!!
    xnext = Mk*[x;y];
    xnext(1:np) = xpnext;
    ynext = Cp*xnext(1:np);
    xcnext = xnext(np+1:end);
    
    % Add noise
    ynext = ynext + noises(:, k + dk + 1);

    % Algorithm 1, line 2: fusion
    tic
    if any(V(:) > 0)  % if there is noise
        X = X(1:np,1:np);
        X = (X+X')/2;

        if TEST_PLOTS && k <= 200
            Ey = Cp*ellipsoid(xptilde,X);
            Ev = ellipsoid(ynext,V);
            plot([Ey,Ev]);
        end

        % Fusion function here
        [xptilde,X] = ellobserverintersection(xptilde,ynext,X,V,Cp);
        X = (X+X')/2;

        if TEST_PLOTS && k <= 200
            Eyy = Cp*ellipsoid(xptilde,X);
            plot(Eyy,'g');
        end
        
        assert((xptilde - xpnext)'*(X\(xptilde - xpnext)) <= 1,...
            'Prediction failed');
        
    else  % Noiseless: Schweppe, 1968 (Appendix D)
        X11 = X(1:pp,1:pp);
        if min(eig(X11)) > 1e-8
            X12 = X(1:pp,pp+1:end);
            X22 = X(pp+1:end,pp+1:end);
            e = ynext - Cp*xptilde;
            X11inv = inv(X11);
            beta2 = 1 - e'*X11inv*e;
            B22 = beta2*(X22 - X12'*X11inv*X12);
            X = zeros(np);
            X(pp+1:end,pp+1:end) = B22;
            xptilde = [ynext; xptilde(pp+1:end) + X12'*X11inv*e];
            
            assert((xptilde(pp+1:end) - xpnext(pp+1:end))'...
                   *(X(pp+1:end,pp+1:end)...
                   \(xptilde(pp+1:end) - xpnext(pp+1:end))) <= 1,...
                   'Prediction failed');
        end
    end
    timesFusion(end+1) = toc;
    timeThisStepSTC = timesFusion(end); 
    Xlog(:,:,end+1) = X;
    
    x = xnext;
    y = ynext;
    k = k+dk;
    if DEBUG
        fprintf('dk = %d, trace = %g\n', dk, trace(X));
    end
end

%% Time statistics
fprintf(fileId, '============== TIMING STATISTICS ==============\n');
fprintf(fileId, 'Pre-process:         %6.2f ms\n', timePreProcess*1000);
fprintf(fileId, '|--Reachability:     %6.2f ms\n', timeWk*1000);
fprintf(fileId, '|--Offline Matrices: %6.2f ms\n', timeOfflineMatrices*1000);
fprintf(fileId, 'Online\n');
fprintf(fileId, '|--Initialization:   %6.2f ms\n', timeInit*1000);
fprintf(fileId, '|--STC (ms):         |  mean  |  min  |  max  |\n');
fprintf(fileId, '   |--Cycle:         | %6.2f | %5.2f | %5.2f |\n',...
    mean(timesSTC)*1000, min(timesSTC)*1000, max(timesSTC)*1000);
fprintf(fileId, '      |--Fusion:     | %6.2f | %5.2f | %5.2f |\n',...
    mean(timesFusion)*1000, min(timesFusion)*1000, max(timesFusion)*1000);
fprintf(fileId, '      |--Eta:        | %6.2f | %5.2f | %5.2f |\n',...
    mean(timesEta)*1000, min(timesEta)*1000, max(timesEta)*1000);
fprintf(fileId, '      |--Prediction: | %6.2f | %5.2f | %5.2f |\n',...
    mean(timesPrediction)*1000, min(timesPrediction)*1000, max(timesPrediction)*1000);

%% Test plots...

if TEST_PLOTS
    try
        close(2); close(3); close(4);
    end
    
    figure(2);
    plot(kslog,dklog);
    hold on;
    plot(kslog, dkzlog);
    legend('STC', 'When PETC would have');

    figure(3);
    plot(klog,xlog);
    
    % Compare output with estimated output
    figure(4);
    plot(klog,ylog);
    hold on;
    plot(klog,xtlog(:,1:np)*Cp','--');

end
