%Compute numerical solution of classical BS PDE using monte carlo reg.
%non constant pde coeff. (make sigma_0 = sigma*x)
%The default absolute error tolerance for these cases is 1e-4


%whos A 
%where A is a matrix tells size in byte of matrix 
%1megabyte = 1e+6 bytes

%use sparse(epsil) to convert matrix into sparse structure

function [abs_error] = bs2d(  n  , L ) 
  
  T=0.5;
  %rng(1); %set seed
  h = T/n ; %step size
  r= 0.1;
  p=1/3;
  x0=[20  33];
  d=length(x0);
  sigma=[.2  .15]; %volatility 
  sigma_sim = sigma; %repelem( mean(sigma) , 2 ) ; %sigma_sim is like sigma bar in the paper,sigma instead of rep mean sigma; worse results
  rho = 0.66; %true correlation
  rho_sim = rho ; %numerical used in PDE function
  K= 5; %strike
  %cov_mat = diag(sigma) + flip ( diag ( repelem( rho * prod(sigma) ,2 ) ) );
  
  rhomat = eye(d,d) + flip ( diag ( repelem( rho_sim  , d ) ) );
  chol_rho = chol (rhomat, 'lower'); %check for upper 
  
  %Generate forward process
  [x,epsil] = xsim(x0, n+1, sigma_sim, chol_rho, p,h,L); 
  
  %Terminal solution
  y = max( min(  x(:,:,n+1) , [], 2) - K , 0);
 
  %Loop for computing backward the solution
  for t = n : -1 : 1
 
   %Kernel 1   
    Sinv = 1./(sigma .* x(:,:,t));
    kernel1 = (1/sqrt(h)) * ( Sinv .* ( chol_rho' \ epsil(:,:,t)' )' )  ;

   %Exact
   % ( inv( eye(2,2) .* [0.3*2,0.2*3] * chol_rho  )' * epsil(:,:,1)' ) /sqrt(h)
   
   %Kernel2
       
   ondiag = epsil(:,:,t).^2 * ( (1-p) - (1-3*p) ) - 2*p ; %diagonal epsilon matrix
   offdiag = (1-p) * ( epsil(:,1,t) .* epsil(:,2,t) ); %offdiagonal epsilon matrix
   
   a1=(chol_rho' \ [ ondiag(:,1), offdiag ]')'; %central multiplication
   a2=(chol_rho' \ [ offdiag, ondiag(:,2) ]')'; %central multiplication 
   
   b1=[a1(:,1),a2(:,1)] / chol_rho; %first column elements 
   b2=[a1(:,2),a2(:,2)] / chol_rho; %second column elements
   
   kernel2on = Sinv.^2 .* [b1(:,1),b2(:,2)] ; %diagonal elements 
   kernel2off = prod ( Sinv , 2) .* b1(:,2) ; %or b2(:,1)
   kernel2 = (1/((1-p)*h)) * [kernel2on, kernel2off];
%Exact 
%    1/(h-p*h)(eye(2,2) .* Sinv(1,:)) * inv_chol_rho' * ...
%    [ epsil(1,1,t).^2 * ( (1-p) - (1-3*p) ) - 2*p , epsil(1,1,t).*epsil(1,2,t) * (1-p) ; epsil(1,1,t).*epsil(1,2,t) * (1-p), epsil(1,2,t).^2 * ( (1-p) - (1-3*p) ) - 2*p ] * inv_chol_rho * (eye(2,2) .* Sinv(1,:))

  %Compute 'y part' for obtain beta of regression

        yreg1 = y .* kernel1 ; 
        yreg2 = y .* kernel2 ; 
      
        if t==1 %compute the option price
              y = mean( y + h * ...
                F( x(:,:,t), y, yreg1, yreg2, r, sigma, rho, rho_sim, sigma_sim ) );
        else 
            %2- Compute basis functions - ex: polynomial 1+x+x.^2

            basis =  [ ones(L,1) , x(:,1,t), x(:,2,t), ...
                        (x(:,1,t).^2),(x(:,2,t).^2),...
                        (x(:,1,t) .* x(:,2,t)) ] ; 
            
            %3- Get beta from linear regression
            regressors = (basis'*basis) \ basis';
            alphaop0 = regressors * y; 
            alphaop1 = regressors * yreg1;
            alphaop2 = regressors * yreg2; %two 3d array with second order derivatives
            
            %4- Get fitted values 
            op0 = basis * alphaop0; 
            op1 = basis * alphaop1; 
            op2 = basis * alphaop2; 

            y =  op0 + h * F( x(:,:,t), op0, op1, op2 , r, sigma, rho,rho_sim, sigma_sim) ;
        end 
  end 

  %results
  sol_exact = JohnsonMinLuca(x0',sigma',r,K,T,rhomat);
  abs_error = abs( y - sol_exact );
 
end


%Simulate forward process x and epsilon's trinomials 

function  [x,epsil] = xsim(x0,n,sigma,chol_rho,p,h,L)
    
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
         x(:,:,t) =  x(:,:,t-1) + sqrt(h) * ( sigma .* x(:,:,t-1) ) .* ( chol_rho * epsil(:,:,t-1)' )'  ;        
    end
end 
%Exact     
%(eye(2,2) .* [0.3*2,0.2*3]) * chol_rho *  epsil(1,:,t-1)'


%get diagonal elements from 3d array ie a page

%Input: A is a square matrix 3d array
%Output: diagonal elements of A for each page

% function [out]  = pagediag( A  ) 
% 
%         d1=size(A,1);d2=size(A,2);d3=size(A,3); 
%         index=sub2ind( size(A) , repmat(1:d1, 1, d3), repmat(1:d2 ,1, d3), repelem(1:d3, d2));
%         out = reshape(A(index),[1,d2,d3]);
%         
% end


%Input: discretized inputs of PDE and parameters
%Output: discretized vector of y
function [F] = F(x,u,du,ddu, r, sigma, rho, rho_sim, sigma_sim)

  %PDE's for Black scholes
  G = (r * x(:,1).* du(:,1) ...
      + r * x(:,2) .* du(:,2) ...
      + 0.5 * sigma(1)^2 * x(:,1).^2 .* ddu(:,1) ...
      + 0.5 * sigma(2)^2 * x(:,2).^2 .* ddu(:,2) ...
      + rho * prod(sigma) * x(:,1) .* x(:,2) .* ddu(:,3) ...
      - r*u  ); 

  %F= G - 0.5 * trace ( (sigma_0*sigma_0^T) * ddu^T ) 
  
  %To make the trace:

  S12 = (sigma_sim .* x).^2 .* ddu(:,[1,2]) ;
  S22 = rho_sim * prod(sigma_sim) .* prod(x,2) .* ddu(:, 3);
  
  tracevec  = sum( S12 + S22, 2);

  F = G - 0.5 * tracevec;
       
end  



% Code for kernel2

% (eye(2,2) .* Sinv(1,:)) * inv_chol_rho' * [ ep1(1,1), ep1(1,2); ep1(1,2), ep2(1,2)  ] * inv_chol_rho *  (eye(2,2) .* Sinv(1,:)) 
% 
% inv_chol_rho' * [ ep1(1,1), ep1(1,2); ep1(1,2), ep2(1,2)  ] * inv_chol_rho 
% 
% 
% ep=rand(L,d);
% ep1 = ep; ep2 = [ep1(:,2),rand(L,1)];
% 
% a1=(inv_chol_rho' * ep1')';
% a2=(inv_chol_rho' * ep2')';
% 
% b1=[a1(:,1),a2(:,1)] * inv_chol_rho; %first column elements 
% b2=[a1(:,2),a2(:,2)] * inv_chol_rho; %second column elements
% 
% 
% Sinv.^2 .* [b1(:,1),b2(:,2)] ;
% Sinv(:,1).* Sinv(:,2) .* [b1(:,2),b2(:,1)]



%% Code for sigma0 in pde
% 
% S = sigma_sim .* x;
%   %2. Matrix product eye (S1,S2) * gamma
%   S1 = S(:,1) .* ddu(:, [1,3]); %first row of matrix product
%   S2 = S(:,2) .* ddu(:, [4,2]); %second row of matrix product
%   %3. Make second product to the right
%   c = chol_rho * chol_rho';
%   B1=(c * [S1(:,1), S2(:,1)]')';
%   B2=(c * [S1(:,2), S2(:,2)]')';
%   %4. Make third product from right
%   prod1 = S(:,1) .* [B1(:,1), B2(:,1)];
%   prod2 = S(:,2) .* [B1(:,2), B2(:,2)];
%   %5. Make trace sum for every L
%   tracevec = sum ( [prod1(:,1), prod2(:,2)] , 2 );
 
% 
% 
%  %exact for first L
% (eye(2,2) .* S(1,:)) * chol_rho * chol_rho' * (eye(2,2) .* S(1,:))' * ...
% [ ddu(1,1) ddu(1,3); ddu(1,3) ddu(1,2)]
% 
