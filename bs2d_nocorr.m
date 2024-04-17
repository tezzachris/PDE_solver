%Compute numerical solution of classical BS PDE using monte carlo reg.
%non constant pde coeff. (make sigma_0 = sigma*x)

%2D Without Correlation

function [abs_error] = bs2d_nocorr( n  , L ) 
  
  T = 0.5; %time horizon
  rng(1); %set seed
  h = T/n ; %step size
  r= .0001; %interest rate
  p= 1/3; %probability of trinomial
  x0=[2 3]; %initial price
  d=length(x0); %dimension number
  sigma=[0.3 0.2]; % sigma of PDE
  sigmasim = sigma;% reshape ( repelem( mean(sigma) , 2 ) , 1, 1, d) ; % sigma for the forward
  K=2.5; %strike
  rhomat=eye(d);
  %Generate forward process
  [x,epsil] = xsim(x0,n+1,sigmasim,p,h,L); 

  %Terminal solution
  
  y = max( min(  x(:,:,end) , [], 2) - K , 0);
  
  for t = n : -1 : 1

  %1- Compute y for obtain beta of regression
        %Kernel1
        sigma0 = sigma .* x(:,:,t);
        kernel1 = epsil(:,:,t)./(sigma0*sqrt(h))  ;
        yreg1 = y .* (kernel1)  ; 

        %Kernel2
        kernel2 = ( sigma0.^2 .* ( epsil(:,:,t).^2 *(1-p) - (1-3*p)*epsil(:,:,t).^2 - 2*p ) )/(h*(1-p)) ;
        yreg2 = y .* (kernel2) ; 
        
        if t==1
            y = mean( y + h * F( x(:,:,t) , y, yreg1, yreg2, sigma0, sigma, r ) );
        else 
            %2- Compute basis functions - ex: polynomial 1+x+x.^2
            basis =  [ ones(L,1), x(:,1,t), x(:,2,t), ...
                                  x(:,1,t).^2, x(:,2,t).^2,...
                                 (x(:,1,t) .* x(:,2,t) )] ; %x dim(n x L), basis (L x 2)
            
            %3- Obtain beta from regression
            alphaop0 = (basis'*basis) \ basis' * y; 
            alphaop1 = (basis'*basis) \ basis' * yreg1 ;
            alphaop2 = (basis'*basis) \ basis' * yreg2 ;
            
            %4- Get regression operators
            op0 = basis * alphaop0; %ie fitted values
            op1 = basis * alphaop1;
            op2 = basis * alphaop2;
            y =  op0 + h * F( x(:,:,t) , op0, op1, op2, sigma0, sigma,r) ;
        end 
  end 
  sol_exact=JohnsonMinLuca(x0',sigma',r,K,T,rhomat);
  abs_error = abs( y - sol_exact );
end


%simulate forward process 
%output x is n x L, eps n-1 x L

function [x,epsil]=xsim(x0,n,sigmasim,p,h,L)

    d = length(x0);
    x = zeros(L,d,n); 
    x(:,:,1)=repmat(x0,L,1); %x0 = [2,1,4,..d] vector of initial prices
    %Generate trinomiale epsilons : n x L
    epsil=rand(L,d,n-1);
    epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive

    for t=2:n
        x(:,:,t) = x(:,:,t-1) + sqrt(h) * ( sigmasim .* x(:,:,t-1) ) .* epsil(:,:,t-1) ; %paper formula forward process
    end 
end 


function F = F(x,u,du,ddu,sigma0,sigma,r)
  %PDE u=f(t,x), sigma=diffusore

  G = sum( r * x .* du + 0.5 * sigma.^2 .* x.^2 .* ddu - r*u , 2 ) ;
     
  %F= G - 0.5 * sigma_0 * ddu
  F = G - 0.5 * sum (  sigma0.^2   .* ddu , 2 ) ; 
end



% function [cmin] = Johnson_2d(S0vet, sigmavet, r, K, T)
% 
%     d2=@(s,x,sigma) ((log(s/x)+ T*( r - 0.5*sigma^2) ) ) / (sigma*sqrt(T)) ;
%     d1=@(s,x,sigma) d2(s,x,sigma) + sigma*sqrt(T);
% 
%     d2prime=@(s,x,sigma) ((log(s) - T*(0.5*sigma^2) ) ) / (sigma*sqrt(T)) ;
%     d1prime=@(s1,s2,sigma12)  ((log(s1/s2) + (0.5*sigma12^2)*T ) ) / (sigma12*sqrt(T)) ;
% 
%     sigma12sq=@(sigma1,sigma2,rho12) (sigma1^2 + sigma2^2 - 2*rho12*sigma1*sigma2);
%     rho112=@(sigma1,rho12,sigma2) (sigma1 - rho12*sigma2 ) / sqrt(sigma12sq(sigma1,sigma2,rho12));
%     rho212=@(sigma1,sigma2,rho12,rho21,rho22) (sigma2^2-rho21*sigma2*sigma1-rho22*sigma2*sigma2+rho12*sigma1*sigma2)/(sqrt(sigma12sq(sigma2,sigma1,0))*sqrt(sigma12sq(sigma2,sigma2,0)));
%     
%     cmin= S0vet(1)*mvncdf( [ d1(S0vet(1), K, sigmavet(1)) , -d1prime(S0vet(1),S0vet(2),sigma12sq(sigmavet(1),sigmavet(2),0)) ],-rho112(sigmavet(1),0,sigmavet(2)) ) +...
%             +S0vet(2)*mvncdf( [d1(S0vet(2), K, sigmavet(2)) , -d1prime(S0vet(2), S0vet(2),sigma12sq(sigmavet(1),sigmavet(2),0)) ], -rho212(sigmavet(1),sigmavet(2),0,0,0))+...
%                 -K*exp(-r*T)*mvncdf([d2(S0vet(1),K,sigmavet(1)),d2(S0vet(2),K,sigmavet(2))],0);
% 
% end 
% 




