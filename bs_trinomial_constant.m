%Solution of BS constant PDE using Trinomial tree
%Trick: convert to pde to heat eq. for having constant coeff.

function [abs_error] = bs_trinomial_constant( n ) %solution function for BS pde using trinomial setting 
  T = 0.5;
  h = T/n ; %intervals width
  xo = 5;
  d = length(xo);
  K = 5;
  r = .01;
  p = 1/3;
  sigma_0 = .03;
  %Initial price in log terms
  x = xfin(log(xo),n,sigma_0,p,h); %final vector of trinomial forward process
  y = zeros( n+1 , 2*n+1 ) ;  %creating option price solution matrix
  
  y(1,:)= exp(-r*T)*max(exp(x)-K,0); % term cond BLACK-SCHOLES, price of the option is the payoff of the option

  G= @(t,u,z,zz,sigma_0,r) ( 1/2*sigma_0^2 * zz  + (r-0.5*sigma_0^2)*z  );
  F = @(t,u,z,zz,sigma_0,r) ( G(t,u,z,zz,sigma_0,r) - 1/2 * trace( (sigma_0^2) * zz) );
  
  for i = 2:(n+1)
      cont=0; %counter for creating the option price matrix y
      subint=0:h:(T-h); %subinterval counter for t in BS formula
      t=subint( n + 2 - i ); % t for BS formula
      %xp = x(i:end-i+1); %stock prices for BS formula 
        for k = transpose(triplets(y(i-1,1:end-2*(i-2)))) %correct error if S=K
          cont=cont+1; %column selctor , xp(cont) gives the cont-th stock price at time t 
          y(i,cont) = expect(transpose(k),p,0,h,sigma_0,d)  + h*F(t,expect(transpose(k),p,0,h,sigma_0,d),expect(transpose(k),p,1,h,sigma_0,d),expect(transpose(k),p,2,h,sigma_0,d), sigma_0,r  ) ;
        end
       
  end
  sol=nonzeros(y(n+1,:));
  sol_exact=bs_price1d(r,sigma_0,0,T,0,K,xo);
  abs_error = abs(sol - sol_exact)/sol_exact;

end

%Generate forward process trinomial

function [y]=xfin(xo,n,sigma_0,p,h)
    y=flip([ flip(xo + (sqrt(h)*sigma_0*(-1/sqrt(p)))*[1:n]),xo,xo + (sqrt(h)*sigma_0*(1/sqrt(p)))*[1:n] ]);
end

%from a row vector (of odd size) construct a matrix of collection of 1x3
%elements

%Generate triplets pairs from a row vector

function ye =triplets(y) 
    ye=zeros((size(y,2)-2),3);
    for i=1:(size(y,2)-2)
        ye(i,:)=y(1,i:i+2);
    end
end

%Expectation

function e = expect(y,p,k,h,sigma_0,d) %y 1x3
    eps=[1/sqrt(p),0, -1/sqrt(p)]; 
    eps_prob=[p/2,(1-p),p/2]; %epsilon trinomial variables
    ke = eye(3,3); %diagonal matrix for multiply with y
    if k==1 
            ke(1:(3+1):end) = (transpose(inv(sigma_0)) * eps ) / sqrt(h); %substituting diagonal
    elseif k==2 
            num=zeros(1,3);
            for i = 1:3
                num(i) = (transpose(inv(sigma_0)) * (   (1-p)* eps(i)^2  - (1-3*p)*diag(eps(i)^2)  - 2*p* eye(d,d) ) * inv(sigma_0) ) / ((1-p)*h);
            end
            ke(1:(3+1):end)=num;
    end
          
    e=dot((y*ke),eps_prob); % dot product( (1x3 diag_3x3) , 1x3)
end

