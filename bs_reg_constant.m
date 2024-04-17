%Compute numerical solution of BS PDE with constant coefficients
%one dimension (one underlying)
%regression and montecarlo
function [res] = bs_reg_constant( xo, n , T, sigma_1, K, r, L, p ) 

  h = T/n ; % dt
  %sigma_0 numerico, sigma_1 vero
  sigma_0=sigma_1;

  %PDE - uso sigma_0 anche per la pde
  %G = @(u,z,zz,sigma_1,r) ( 1/2*sigma_1^2 * zz  + (r-0.5*sigma_1^2)*z  );  
  %F = @(u,z,zz,sigma_1,r) ( G(u,z,zz,sigma_1,r) - 1/2 * trace( (sigma_1^2) * zz) );
  
  %Generate forward process
  [x,eps_sim]=xsim(log(xo),n,sigma_0,p,h,L); %note log(xo)
  %Terminal solution
  y =  exp(-r*T) * max( exp( x(end,:)) - K , 0 ) ;
 
  %Kernels
  kernel0= ones(L,1);
  kernel1 = (sigma_0') \ eps_sim  ./ sqrt(h);
  
  kernel1 = flip (kernel1); %start from n-1 to 1

  for i = 1:n
  %1- Compute y for obtain beta of regression
        yreg0 = y' .* kernel0 ;
        yreg1 = y' .* kernel1(i,:)' ; 
       
        if i==n
            y = mean( yreg0 + h .*(r- 0.5 *sigma_1.^2)*yreg1 );
        else 
             %2- Compute basis functions - ex: polynomial 1+x+x^2
            basis =  [ones(L,1), x(end-i,:)',transpose(x(end-i,:).^2)] ; %x dim(n x L), basis (L x 2)
            %3- Obtain beta from regression
            alphaop0 = (basis'*basis) \ basis' * yreg0; %
            alphaop1 = (basis'*basis) \ basis' * yreg1; % beta (2 x 1)

            %4- Get regression operators
            op0 = basis * alphaop0; 
            op1 = basis * alphaop1;
            y = transpose( op0 + h*(r-0.5*sigma_1.^2)*op1 );
        end 
  end 
  
  sol_exact=bs(r,sigma_1,0,T,0,K,xo);
  res=table(y,sol_exact,VariableNames={'Numerical';'Exact'});


end


%simulate forward process 
%output x is n x L, eps n-1 x L

function [x,epsil]=xsim(x0,n,sigma_0,p,h,L)
    %Generate trinomial rv's NOT WORKS WITH P>1
    n=n+1; %add one because x0 is at 1
    epsil=rand(n-1,L); 
    epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive
    %Generate forward process x
    x=zeros(n,L); x(1,:)=x0; % forward process nxL
    for j = 1:L
        for i=2:n
            x(i,j) = x(i-1,j) + sqrt(h)*sigma_0*epsil(i-1,j); %paper formula forward process
        end 
    end
end


