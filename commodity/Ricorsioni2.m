
%Verify recursion eq.11 Schwartz 1998

clear;

%Maturity
T=30; 
%Runge kutta
nint=2000; %densita' griglia
%Parametri comuni
h=T/nint; 

%Parametri simuazione
h=1;
L=10000; %nsimulazioni
n=5000; %npoint

%Parameters Schwartz

A = [0; 
    -0.0535;
    0.2940*2.7471; 
    0.1024*0.2829];

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

k11=1.1374;mu2=-0.0535; x0=(k11*mu2);
theta0=mu2;
k33=0.2940; mu3=2.7471; qt0 = k33*mu3;
k44=0.1024; mu4=0.2829; vt0=k44 * mu4;

%Initial prices
Xt0 = [x0 theta0 qt0 vt0];
X = repmat(Xt0,L,1); 

for i = 2 : n 
    Z = randn(L,d);
    drift = A + B * X';
    %attenzione sarebbe chol(omega0 + omega1 * vt) * Z ho semplificato
    %togliendo omega0 per velocizzare il codice
    %sqrt(2) * chol(A) = chol(2*A)
    if isnan(X)
        i
    end
    X = X + h * drift' + (sqrt(h) * sqrt(abs(X( : , d )))' .* (chol( omega1 )  * Z'))' ;  
    X(:,d) = abs(X( :, d ));
end

%prezzo future simulazione monte carlo
Fsim=mean(exp(X(1,:))); 

%prezzo future ricorsioni
[a,b]=rungekutta(nint,T,A,B,omega0,omega1);
Frecursion = exp(  a(1) + sum ( b(1,:) .* Xt0 ) ); 

[Fsim, Frecursion]
