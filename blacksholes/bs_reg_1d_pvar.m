%Compute numerical solution of classical BS PDE using monte carlo reg.
%non constant pde coeff. (make sigma_0 = sigma*x)
%one dimension
function [abs_error] = bs_reg_1d_pvar( n , L ) 
  
  rng(1); %set seed
  x0=5;
  sigma=0.3; %true sigma
  sigma_sim = sigma ;  %numeric sigma
  r = 0.01;
  p = 1/3;
  K=x0;
  T=0.5;
  h = T/n; 
     
  %Generate forward proces
  [x,epsil] = xsim(x0,n+1,sigma_sim,h,L,r); 

  %Terminal solution
  y =  max( x(n+1,:) - K , 0 ) ;  %alternative use exp(-r*T) 

  for t = n   : -1 : 1
        
        %1- Compute response vector (kernel*y) to put in the regression
        yreg0 = y'  ; %Kernel0 is just ones

        %Kernel1
        sigma0 = x(t,:) * sigma;
        kernel1 = ( sigma0 .\ epsil(t,:) ) / sqrt(h);
        yreg1 = y' .* kernel1' ; 

        %Kernel2
        kernel2 = ( (sigma0.^2) .\ ( (1-p) * epsil(t,:).^2  - ...
            (1-3*p)*( epsil(t,:).^2)  - 2*p ) )  / ((1-p)*h) ;
 
        yreg2 = y' .* kernel2' ;
       
        if t==1
            y = mean( yreg0 + h * F( x(t,:)',yreg0,yreg1,yreg2,sigma0',sigma,r ) );
        else 
            %2- Compute basis functions - ex: polynomial 1+x+x^2 ie matrix
            %of regressors
            basis =  [ ones(L,1), x(t,:)', (x(t,:).^2)' ] ; 
            %3- Obtain beta coefficients from regression
            regressors = (basis'*basis) \ basis';
            alphaop0 = regressors * yreg0; 
            alphaop1 = regressors * yreg1; 
            alphaop2 = regressors * yreg2;
            %4- Get regression operators
            op0 = basis * alphaop0; %predicted values of yreg0
            op1 = basis * alphaop1;
            op2 = basis * alphaop2;
            y = transpose( op0 + h * F( x(t,:)',op0,op1,op2,sigma0',sigma,r) );
        end 
  end 
  
  sol_exact=bs_close_formula(r,sigma,0,T,0,K,x0);
  abs_error= abs( y-sol_exact );
  %res=table(y,sol_exact,abs_error,VariableNames={'Numerical';'Exact';'Absolute Error'});
   
end


%simulate forward process 
%output x is n x L, eps n-1 x L

function [x,epsil]=xsim(x0,n,sigma_sim,h,L,r)
    
    %Trinomials
    pmax = 1/3;
    q = r : (pmax-r)/(n-1) : pmax ; 
    epsil2 = rand(n-1,L); 
    epsil = zeros(n-1,L);
%     epsil( epsil <= r/sigma_sim ) = -2*r;
%     epsil( epsil > r & epsil <= 2*r ) = r;
%     epsil( epsil > 2*r ) = 2*r;

    %epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    %epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    %epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive
%     %Generate forward process 
    x=zeros(n,L); 
    x(1,:)=x0; 
    
    for i=2:n
        epsilon = epsil2(i-1,:);
        p = q(i);
        epsilon( epsilon <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
        epsilon( epsilon > p/2 & epsilon <= (1-p/2) ) = 0; %1-p prob to 0
        epsilon( epsilon > (1-p/2) ) = 1/sqrt(p); 
        epsil(i-1,:) = epsilon;
        sigma0 = sigma_sim * x(i-1,:);
        x(i,:) = x(i-1,:) + sqrt(h) * sigma0 .* epsil(i-1,:); %paper formula forward process
    end    
end

%PDE 
function F = F(x,u,du,ddu,sigma0,sigma,r)

    G = r * x .* du  + 0.5 * sigma^2 * x.^2 .*ddu - r*u;  
    F = G  - 0.5 * sigma0.^2  .* ddu ;

end





