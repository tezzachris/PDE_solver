
%Verify recursion eq.11 Schwartz 1998

clear;

T=0.5; %maturity, solve for [0,T]
nint=2000; %densita' griglia
h=T/nint; %dt step size 
L=1000; %nsimulazioni
n=nint+1; %npoint

%Parameters Schwartz

A = [0; 
    -0.0535;
    0.2940*2.7471; 
    0.1024*0.2829];

B = [-1.1374 1.1374 -1 2.1302;
     0  0   0           0;
     0  0   -2.7471     0;
     0  0   0       -0.2829];

omega0 = zeros(4,4);
omega1 = [1.0000 0.84000 -0.0405 0.02300;
          0.8400 1.1865  -0.3350 -0.2148;
          -0.0405 -0.335 0.4411  0.1733;
          0.02300 -0.2148 0.1733 0.1227];

d=length(A);
X = zeros(n,d);

k11=1.1374;mu2=-0.0535; x0=(k11*mu2);
theta0=mu2;
k33=0.2940; mu3=2.7471; qt0 = k33*mu3;
k44=0.1024; mu4=0.2829; vt0=k44 * mu4;

%Initial prices

X(1,:) = [ x0 theta0 qt0 vt0 ]; 
xTsim = zeros(1,L);

for j = 1 : L
    Z = randn(n-1,d);
    for i = 2 : n 
        media = A + B * X( i-1 , : )';
        X(i,:) = X(i-1,:) + h * media' + ( sqrt(h) * chol( omega0 + omega1 .* abs(X( i-1, d )) ) * Z(i-1,:)' )';  
        X(i,d) = abs(X( i, d ));
    end
    xTsim(j)=X(n,1);
end

%prezzo future simulazione monte carlo
Fsim=mean(exp(xTsim)); 


%Runge kutta

da = @(t,b)  b*A + 1/2 * b*omega0*b';
db = @(t,b) [b*B(:,1), b*B(:,2), b*B(:,3), b*B(:,4) + 1/2 * b*omega1*b'];

%Initial conditions

t=zeros(n,1);
t(1)=0;
a= zeros(n,1);
a(1)=0; 
b= zeros(n,d);
b(1,:)=[1 0 0 0];

for i = 1  : nint %non fino n

    t(i+1) = t(i) + h;
    
    a1 = da(t(i), b(i,:));
    b1 = db(t(i), b(i,:));

    a2 = da( t(i)+h/2, b(i,:)+h/2*a1 );
    b2 = db( t(i)+h/2, b(i,:)+h/2*b1 );

    a3 = da(t(i)+h/2, b(i,:)+h/2*a2 );
    b3 = db(t(i)+h/2, b(i,:)+h/2*b2 );
    
    a4 = da(t(i)+h,   b(i,:)+h*a3 );
    b4 = db(t(i)+h,   b(i,:)+h*b3 );

    a(i+1) = a(i) + h/6 * (a1 + 2*a2 + 2*a3 + a4);
    b(i+1,:) = b(i,:) + h/6 * (b1 + 2*b2 + 2*b3 + b4);

end 

%prezzo future ricorsioni
Frecursion = exp(  a(end) + sum ( b(end,:) .* X(1,:) ) ); 

[Fsim, Frecursion]
