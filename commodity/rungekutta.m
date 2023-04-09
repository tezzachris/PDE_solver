
function [a,b] = rungekutta(nint,T,A,B,omega0,omega1)

d = size(A,1);
n = nint + 1;
h = T/nint;

da = @(t,b)  b*A + 1/2 * b*omega0*b';
db = @(t,b) [b*B(:,1), b*B(:,2), b*B(:,3), b*B(:,4) + 1/2 * b*omega1*b'];

%Initial conditions

t=zeros(n,1);
t(1)=0;
a= zeros(n,1);
a(1)=0; 
b= zeros(n,d);
b(1,:)=[1 0 0 0];

%Fourth order

for i = 1  : nint %perchè t=h*(1:nint)

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

%to order for time to maturity
a = flip(a);
b = flip(b);

end