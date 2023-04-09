

%Codice per simulare il prezzo di un future scadenza 30gg

T=30; %scadenza future

%Runge Kutta parametri
nint=2000; % numero di intervalli
h=1;%T/nint; %dt grandezza intervallo
npunti=nint+1; %numero punti griglia

%Future price simulation parameters
n = 500; %dati future 


%Parametri Schwartz

A = [0; 
    -0.0535;
    0.2940*2.7471; 
    0.1024*0.2829];

%dimensione vettore stati
d=length(A);

B = [-1.1374 1.1374 -1 2.1302;
     0  0   0           0;
     0  0   -2.7471     0;
     0  0   0       -0.2829];

omega0 = zeros(d,d);
omega1 = [1.0000 0.84000 -0.0405 0.02300;
          0.8400 1.1865  -0.3350 -0.2148;
          -0.0405 -0.335 0.4411  0.1733;
          0.02300 -0.2148 0.1733 0.1227];

%vettore stati
X = zeros(n,d);

k11=1.1374;mu2=-0.0535; x0=(k11*mu2);
theta0=mu2;
k33=0.2940; mu3=2.7471; qt0 = k33*mu3;
k44=0.1024; mu4=0.2829; vt0=k44 * mu4;

%valore iniziale vettore stati
X(1,:) = [ x0 theta0 qt0 vt0 ]; 
%moti browniani
Z = randn(n-1,d);
%discretizzazione sistema SDE
for i = 2 : n 
    X(i,:) = X(i-1,:) + h * (A + B * X( i-1 , : )')' + ( sqrt(h) * chol( omega0 + omega1 .* abs(X( i-1, d )) ) * Z(i-1,:)' )';  
    X(i,d) = abs(X( i, d ));
end

%Runge kutta risolve sistema ODE dei parametri affini del prezzo del future
[a,b]=rungekutta(nint,T,A,B,omega0,omega1);

%Log prezzo del future scadenza T
logFt = (  a + sum ( b .* X , 2 )); 

%Equazione osservazioni
H = 0.0001; %measurement error variance
yt = logFt + sqrt(H)*randn(n,1); %per avere likelihood gaussiana

