%Compute numerical solution of classical BS PDE using monte carlo reg.
%non constant pde coeff. (make sigma_0 = sigma*x)

function [abs_error] = schwartz_3d(  n  , L ) 
  
  T=0.5;
  h = T/n ; %step size
  r=0.001;
  p=1/3;
  x0=[2,3,4];
  d=length(x0);
  sigma=[0.3,0.2,0.1];
  sigma_sim = sigma;
  K=3;
  rho = [0.5,0.3,0.2]; %true correlation
  rho_sim = rho;

  %create covariance matrix
  A=tril(ones([d,d]),-1); %taglia la matrice dxd
  A(A>0) = rho_sim';
  rhomat = A + A' + eye(d);
  chol_rho = chol (rhomat, 'lower'); %check for upper 
  
  %Generate forward process
  [x,epsil] = xsim(x0, n+1, sigma_sim, chol_rho, p,h,L); 

  %Terminal solution
  y = max( min(  x(:,:,n+1) , [], 2) - K , 0);

  %Kernels
  %m = factorial(d) / ( factorial(d-2) * factorial(2) ); % lower triangular matrix elements of cov matrix, 2 is pairs
  n_combs = nchoosek( 1:d ,2);

  for t = n : -1 : 1 
        
        Sinv = 1./(sigma .* x(:,:,t));
        kernel1 = (1/sqrt(h)) * ( Sinv .* ( chol_rho' \ epsil(:,:,t)' )' )  ;
        yreg1 = y .* kernel1 ; 
        
        ondiag = epsil(:,:,t).^2 * ( (1-p) - (1-3*p) ) - 2*p ; %diagonal epsilon matrix
        offdiag = (1-p) * epsil(:,n_combs(:,1),t) .* epsil(:,n_combs(:,2),t) ;
        
        a1=(chol_rho' \ [ ondiag(:,1) , offdiag(:,1), offdiag(:,2)]')'; %central multiplication
        a2=(chol_rho' \ [ offdiag(:,1), ondiag(:,2), offdiag(:,3)]' )';
        a3=(chol_rho' \ [ offdiag(:,2), offdiag(:,3) , ondiag(:,3)]')';

        b1=[a1(:,1),a2(:,1),a3(:,1)] / chol_rho; %first column elements 
        b2=[a1(:,2),a2(:,2),a3(:,2)] / chol_rho;
        b3=[a1(:,3),a2(:,3),a3(:,3)] / chol_rho;
        
        kernel2on = Sinv.^2 .* [b1(:,1),b2(:,2),b3(:,3)] ;
        kernel2off = (Sinv(:,n_combs(:,1)) .* Sinv(:,n_combs(:,2))) .* [b1(:,2),b1(:,3),b2(:,3)]; %or b2(
        kernel2 = (1/(h-h*p)) * [kernel2on, kernel2off];

        yreg2 = y.* kernel2 ;
                 
        if t==1
             y = mean( y + h * ...
                F( x(:,:,t), y, yreg1, yreg2, r, sigma, rho, rho_sim, sigma_sim ) );
           else 
            %2- Compute basis functions - ex: polynomial 1+x+x.^2
            basis =  [ ones(L,1) , x(:,:,t), x(:,:,t).^2, ...
                        x(:,n_combs(:,1),t) .* x(:,n_combs(:,2),t) ] ; 
          
            %3- Obtain beta from regression
            regressors = (basis'*basis) \ basis';
            alphaop0 = regressors * y; 
            alphaop1 = regressors * yreg1;
            alphaop2 = regressors * yreg2 ;
            
            %4- Get regression operators
            op0 = basis * alphaop0; 
            op1 = basis * alphaop1;
            op2 = basis * alphaop2;

            y =  op0 + h * F( x(:,:,t), op0, op1, op2 , r, sigma, rho, sigma_sim, rho_sim) ;
        end 
  end 

  sol_exact=JohnsonMinLuca(x0',sigma',r,K,T,rhomat);
  abs_error = abs(y-sol_exact);

end



%Input: discretized inputs of PDE and parameters
%Output: discretized vector of y
function [F] = F(x,u,du,ddu, r, sigma, rho, sigma_sim, rho_sim )

  d = size(x,2); 
  n_combs = nchoosek( 1:d ,2);
  m = factorial(d) / ( factorial(d-2) * factorial(2) );

  G = sum (r * x .* du ...
      + 0.5 * (sigma .* x ).^2 .* ddu(:,1:d) ...
      + rho .* sigma(n_combs(:,1)) .* sigma(n_combs(:,2)) .* x(:,n_combs(:,1)) .* x(:,n_combs(:,2)) .* ddu(:, d+1 : d+m ) ...
      - r*u , 2); 

  Sdiag = (sigma_sim .* x).^2 .* ddu(:,1:d);
  Smix = rho_sim .* sigma_sim(n_combs(:,1)) .* sigma_sim(n_combs(:,2)) .* x(:,n_combs(:,1)) .* x(:,n_combs(:,2)) .* ddu(:, d+1 : d+m ) ;
  
  tracevec = Sdiag(:,1) + Smix(:,1) + Smix(:,2) + Sdiag(:,2) + Smix(:,1) + Smix(:,3)  + Sdiag(:,3) + Smix(:,2) + Smix(:,3) ;

  F = G - 0.5 * tracevec;
     
end  


  


function [x,epsil]=xsim(x0,n,sigma_sim,chol_rho,p,h,L )
    %Generate forward process x : n x L
    d = length(x0);
    x = zeros(L,d,n); 
    x(:,:,1)=repmat(x0,L,1); %x0 = [2,1,4,..d] vector of initial prices
  
    %Generate trinomiale epsilons : n x L
    epsil=rand(L,d,n-1);
    epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive

    for t=2:n 
         x(:,:,t) =  x(:,:,t-1) + sqrt(h) * ( sigma_sim .* x(:,:,t-1) ) .* ( chol_rho * epsil(:,:,t-1)' )'  ;        
    end
   
end 



