function LL = likelihood_kalman(theta, yt, T, H)

%needs Xt0, Var(Xt0)=Pt0

n = length(yt);
nint = n - 1;
h = T / nint;

%xt0=[x0 theta0 qt0 vt0];
xt0=[ -0.0609   -0.0535    0.8076    0.0290];
d = size(xt0,2);

A = [0;
    theta(1);
    theta(2);
    theta(3)];
indice = [1 5 9 11 13 16];
B = zeros(d,d);
B(indice)=(theta(size(A,1):size(A,1)+length(indice)-1))';
omega0 = zeros(d,d);
omega1 = theta(size(A,1)+length(indice) : end);
omega1 = reshape(omega1,d,d);

%Runge kutta

[a,b]=rungekutta(nint,T,A,B,omega0,omega1);

%Kalman Filtering

%Initialize state vector
Xthat = zeros(n,d);
Xthat(1,:) = xt0; 
%Appears in prediction error
Xtcond=zeros(n,d);

%Log likelihood points
ll = zeros(n,1);
varxt0=h*(omega0+omega1*0.0290);
Pthat= varxt0;

for i = 2 : n 

    %Prediction equations
    Xtcond(i,:) = Xthat(i-1,:) + (A*h + h * B * Xthat(i-1,:)')';
    Z = b(i,:);
    ytcond = a(i) +  Z * Xtcond(i,:)';
    Ptcond = h^2*B*Pthat*B' + h*(omega0 + omega1*abs(Xthat(i-1,4)));
    %Prediction error
    vt = yt(i) - ytcond;
    Ft = Z * Ptcond * Z' + H;

    %Updating equations
    Xthat(i,:) =  Xtcond(i,:) + (Ptcond*Z'*inv(Ft)*vt)';  
    Pthat =  Ptcond - Ptcond*Z'*inv(Ft)*Z*Ptcond;

    ll(i)= -0.5*( log(2*pi) + log(Ft) + (vt)^2 / Ft );

end

LL = - sum(ll);
    
end