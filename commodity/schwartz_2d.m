


function [abs_error] = schwartz_2d(  n  , L ) 
  
  T=.76 ; %maturity of contract
  h = T/n ; %step size
  p=1/3;
  x0=[21 .10]; %commodity price, convenience yield
  d=length(x0);
  sigma=[0.393  0.527]; %volatility 
  sigma_sim = sigma; %repelem( mean(sigma) , 2 ) ; %sigma_sim is like sigma bar in the paper,sigma instead of rep mean sigma; worse results
  rho = 0.766; %correlation
  rho_sim = rho; 
    
  rhomat = eye(d,d) + flip ( diag ( repelem( rho_sim  , d ) ) );
  chol_rho = chol (rhomat, 'lower');
 
  mu = 0.142; alpha=0.106; lambda=0.198; k=1.876; 

  %Generate forward process
  [x,epsil] = xsim(x0, n+1, sigma_sim, chol_rho, p, h, L); 
  
  %Terminal future price
  y = x(:,1,n+1) ; %future price at T = commodity price at T
 
  %Loop for computing backward the solution
  for t = n : -1 : 1
   %Kernel 1   
    Sinv = 1./[ sigma(1) .* x(:,1,t), sigma(2)*ones(L,1)];

    kernel1 = (1/sqrt(h)) * ( Sinv .* ( chol_rho' \ epsil(:,:,t)' )' )  ;

   %Exact
   % ( inv( eye(2,2) .* [0.3*2,0.2*3] * chol_rho  )' * epsil(:,:,1)' ) /sqrt(h)
   
   %Kernel2

   ondiag = epsil(:,:,t).^2 * ( (1-p) - (1-3*p) ) - 2*p ; %diagonal epsilon matrix
   offdiag = (1-p) * ( epsil(:,1,t) .* epsil(:,2,t) ); %offdiagonal epsilon matrix
   
   a1=( chol_rho' \ [ ondiag(:,1), offdiag ]')'; %central multiplication
   a2=( chol_rho' \ [ offdiag, ondiag(:,2) ]')'; %central multiplication 
   
   b1=[a1(:,1),a2(:,1)] / chol_rho; %first column elements 
   b2=[a1(:,2),a2(:,2)] / chol_rho; %second column elements
   
   kernel2on = Sinv.^2 .* [b1(:,1),b2(:,2)] ; %diagonal elements 
   kernel2off = Sinv(:,1) .* Sinv(:,2) .* b1(:,2);
   kernel2 = (1/((1-p)*h)) * [kernel2on, kernel2off];

%Exact 
%    1/(h-p*h)(eye(2,2) .* Sinv(1,:)) * inv_chol_rho' * ...
%    [ epsil(1,1,t).^2 * ( (1-p) - (1-3*p) ) - 2*p , epsil(1,1,t).*epsil(1,2,t) * (1-p) ; epsil(1,1,t).*epsil(1,2,t) * (1-p), epsil(1,2,t).^2 * ( (1-p) - (1-3*p) ) - 2*p ] * inv_chol_rho * (eye(2,2) .* Sinv(1,:))

  %Compute 'y part' for obtain beta of regression
        yreg0 = y  ; %ie observed value 
        yreg1 = y .* kernel1 ; 
        yreg2 = y .* kernel2 ; 
      
        if t==1 %compute the option price
              y = mean( yreg0 + h * ...
                F( x(:,:,t), yreg0, yreg1, yreg2, sigma,sigma_sim,rho,rho_sim,mu,alpha,lambda,k ) );
        else 
            %2- Compute basis functions - ex: polynomial 1+x+x.^2

            basis =  [ ones(L,1) , x(:,1,t), x(:,2,t), ...
                                    (x(:,1,t).^2),(x(:,2,t).^2),...
                                    (x(:,1,t) .* x(:,2,t)) ] ; 
            
            %3- Get beta from linear regression
            regressor = (basis'*basis) \ basis';
            alphaop0 = regressor * yreg0; 
            alphaop1 = regressor * yreg1;
            alphaop2 = regressor * yreg2; %two 3d array with second order derivatives
            
            %4- Get fitted values 
            op0 = basis * alphaop0; 
            op1 = basis * alphaop1; 
            op2 = basis * alphaop2; 

            y =  op0 + h * F( x(:,:,t), op0, op1, op2, sigma,sigma_sim,rho,rho_sim,mu,alpha,lambda,k) ;
        end 
  end 
 
  sol_exact = close_formula(x0,T,k,mu,alpha,lambda,sigma,rho);
  abs_error = abs( y - sol_exact );
  
end


%Simulate forward process x and the trinomials 

function  [x,epsil] = xsim(x0,n,sigma,chol_rho,p,h,L)
    
    %Generate forward process x 
    d = length(x0);
    x = zeros(L,d,n); 
    x(:,:,1)=repmat(x0,L,1); %x0 = [2,1,4,..d] vector of initial prices
  
    %Generate trinomiale epsilons : n x L
    epsil=rand(L,d,n-1);
    epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive
    
    for t=2:n 
         sigma_part = [ sigma(1) .* x(:,1,t-1) , sigma(2)*ones(L,1) ];
         x(:,:,t) =  x(:,:,t-1) + sqrt(h) * sigma_part .* ( chol_rho * epsil(:,:,t-1)' )'  ;        
    end
end 


%Input:
%Output: a vector of y
function [F] = F(x,u,du,ddu,sigma,sigma_sim,rho,rho_sim,mu,alpha,lambda,k)

  %PDE's
  G =  ( ( mu - x(:,2) ) .* x(:,1) .* du(:,1) ...
      + (k*(alpha-x(:,2))-lambda) .* du(:,2) ... %drift
      + 0.5 * ( sigma(1) * x(:,1) ).^2 .* ddu(:,1) ...  %variance
      + 0.5 * sigma(2)^2 .* ddu(:,2)  ...
      + rho * prod(sigma) * x(:,1) .* ddu(:,3) )  ; 

  %F= G - 0.5 * trace ( (sigma_0*sigma_0^T) * ddu^T ) 

  S12 = [ (sigma_sim(1) * x(:,1)).^2 .* ddu(:,1), sigma_sim(2)^2 .* ddu(:,2)  ] ;
  S22 = rho_sim * prod(sigma_sim) .* x(:,1) .* ddu(:, 3);
  
  tracevec  = sum( S12 + S22 , 2);

  F = G - 0.5 * tracevec ;  
end  


function Ft0 = close_formula(x0,T,k,mu,alpha,lambda,sigma,rho)
    
    S0 = x0(1);
    delta0 = x0(2);

    alphah = alpha - lambda/k;

    AT =  ( mu - alphah + 1/2 * (sigma(2)/k)^2 - (sigma(1)*sigma(2)*rho)/k )  * T ...
        + 1/4*sigma(2)^2 * ((1-exp(-2*k*T))/(k^3)) +...
        ( alphah*k + (sigma(1)*sigma(2)*rho) -  (sigma(2)^2)/k ) * ((1-exp(-k*T))/(k^2)) ;

    Ft0 = S0 * exp( -delta0*(1-exp(-k*T))/k + AT ) ;

end


