%Compute numerical solution of classical BS PDE using monte carlo reg.
%non constant pde coeff. (make sigma_0 = sigma*x)
%The default absolute error tolerance for these cases is 1e-4


%whos A 
%where A is a matrix tells size in byte of matrix 
%1megabyte = 1e+6 bytes

%use sparse(epsil) to convert matrix into sparse structure

function [abs_error] = schwartz_2d_old(  n  , L ) 
  
  T=0.5;
  h = T/n ; %step size
  p=1/3;
  x0=[20 .1 ];
  d=length(x0);
  sigma=[0.3  0.15]; %volatility 
  sigma_sim = sigma; %repelem( mean(sigma) , 2 ) ; %sigma_sim is like sigma bar in the paper,sigma instead of rep mean sigma; worse results
  sigma_sim = reshape(sigma_sim,1,1,d);
  rho = 0.1; %correlation
  rho_sim = rho; 
    
  rhomat = eye(d,d) + flip ( diag ( repelem( rho_sim  , d ) ) );
  chol_rho = chol (rhomat, 'lower');
 
  r=0.01; alpha=0.02; lambda=.002; k=0.01; delta=0.04;


  %Generate forward process
  [x,epsil] = xsim(x0, n+1, sigma_sim, chol_rho, p,h,L); 
  
  %Terminal solution
  y = pagetranspose( x(end,:,1)  ); %takes the min(S1,S2) if S01<S02 is only S1
 
  %Loop for computing backward the solution
  for t = n : -1 : 1
    %Sigma0  
    sigma0= times( repmat(eye(d,d),[1,1,L]), permute([sigma_sim(1)*x(t,:,1);sigma_sim(2)*ones(1,L)],[1 3 2])  );
    sigma0= pagemtimes ( sigma0 , repmat(chol_rho,[1,1,L]));  
    epsil_temp = pagetranspose (permute( epsil( t , : ,:), [1,3,2]));
   %Kernel 1
    kernel1 = permute ( pagemtimes ( pagetranspose ( pageinv( sigma0 ) ) , epsil_temp ) * 1/sqrt(h), [3,2,1]);
   %Kernel 2
    cross_prod = pagemtimes(epsil_temp, pagetranspose(epsil_temp));
    kernel2 = ( 1 / ((1-p)*h) ) * ...
        pagemtimes ( ...
        pagemtimes ( pagetranspose ( pageinv(sigma0) ) , ...
        (1-p) * cross_prod - ...
            (1-3*p) * times( repmat( eye(d,d),[1,1,L] ) , pagediag(cross_prod) ) - ... %pagediag is my function
               2*p*eye(d,d) ) ...
                    ,  pageinv(sigma0) ) ;


  %Compute 'y part' for obtain beta of regression
        yreg0 = y  ; %ie observed value 
        yreg1 = y .* kernel1 ; 
        yreg2 = reshape(y,1,1,L) .* kernel2 ; 
      
        if t==1 %compute the option price
            ddumat = yreg2;
            yreg2vec(:,:,1)=reshape(yreg2(1,:,:),size(yreg2(1,:,:),2),[])';
            yreg2vec(:,:,2)=reshape(yreg2(2,:,:),size(yreg2(2,:,:),2),[])';
            y = mean( yreg0 + h * ...
                F( pagetranspose(x(t,:,:)), yreg0, yreg1, yreg2vec, ddumat, ...
                                                sigma0, sigma,rho,r,alpha,lambda,k,delta ) );
        else 
            %2- Compute basis functions - ex: polynomial 1+x+x.^2

            basis =  [ ones(L,1), x(t,:,1)', x(t,:,2)', ...
                        (x(t,:,1).^2)',(x(t,:,2).^2)',...
                        (x(t,:,1) .* x(t,:,2))' ] ; 
            
            %3- Get beta from linear regression
            op0 = (basis'*basis) \ basis' * yreg0; 
            op1 = pagemtimes( (basis'*basis) \ basis' , yreg1 );
            op2 = pagemtimes( (basis'*basis) \ basis' , permute( yreg2, [3,2,1]) ); %two 3d array with second order derivatives
            
            %4- Get fitted values 
            op0 = basis * op0; 
            op1 = pagemtimes( basis , op1 );
            op2 = pagemtimes( basis , op2 ); %dduvec
            ddumat = [];
            ddumat(1:2:2*L,:) = op2(:,:,1);
            ddumat(2:2:2*L,:) = op2(:,:,2);
            ddumat = reshape(ddumat' , [d,d,L]);
            y =  op0 + h * F( pagetranspose(x(t,:,:)), op0, op1, op2, ddumat, sigma0,...
                sigma,rho,r,alpha,lambda,k,delta) ;
        end 
  end 

  %results
 
  sol_exact=close_formula(x0(1),T,k,r,alpha,lambda,sigma,rho,delta);
  abs_error = abs( y-sol_exact );
  
end


%Simulate forward pr ocess x and epsilon's trinomials 

function  [x,epsil] = xsim(x0,n,sigma_sim,chol_rho,p,h,L)
    %Generate trinomial rv's NOT WORKS WITH p>1
    
    %Generate forward process x : n x L
    d = length(x0);
    x = zeros(n,L,d); 
    x(1,:,:)=repmat(x0,L,1); %x0 = [2,1,4,..d] vector of initial prices

    %Generate trinomiale epsilons : n x L

%% (L,d,n) 

    epsil=rand(n-1,L,d);
    epsil( epsil <= p/2 ) = -1/sqrt(p); %p/2 prob to negative
    epsil( epsil > p/2 & epsil <= (1-p/2) ) = 0; %1-p prob to 0
    epsil( epsil > (1-p/2) ) = 1/sqrt(p); %p/2 prob to positive

    for t=2:n 
        sigma0= times( repmat(eye(d,d),[1,1,L]), permute([sigma_sim(1)*x(t-1,:,1);sigma_sim(2)*ones(1,L)],[1 3 2]) );
        sigma0= pagemtimes ( sigma0 , repmat(chol_rho,[1,1,L]));  
        epsil_temp = pagetranspose ( permute( epsil( t-1, : ,:), [1,3,2]) );
        prodotto = pagemtimes ( sigma0, epsil_temp); %product sigma_0 * epsilon
        prodotto = pagetranspose(permute( prodotto , [3,2,1])); 
        %Forward process 
        x(t,:,:) = x(t-1,:,:) + sqrt(h) *  prodotto;

    end
end 
      

%get diagonal elements from 3d array ie a page

%Input: A is a square matrix 3d array
%Output: diagonal elements of A for each page

function [out]  = pagediag( A  ) 

        d1=size(A,1);d2=size(A,2);d3=size(A,3); 
        index=sub2ind( size(A) , repmat(1:d1, 1, d3), repmat(1:d2 ,1, d3), repelem(1:d3, d2));
        out = reshape(A(index),[1,d2,d3]);
        
end


%Input:
%Output: a vector of y
function [F] = F(x,u,du,ddu,ddumat,sigma0,sigma,rho,r,alpha,lambda,k,delta)
  L = size(x,1);
  %PDE's
  G = (r - delta ) * x(:,:,1) .* du(:,:,1) + (k*(alpha-x(:,:,2))-lambda) .* du(:,:,2) ... %drift
      + 0.5 * sigma(1)^2 * (x(:,:,1).^2) .* ddu(:,1,1) ...  %variance
      +  0.5 * sigma(2)^2 * (x(:,:,2).^2) .* ddu(:,2,2) + ...
      rho * prod(sigma) * x(:,:,1) .* x(:,:,2) .* ddu(:,2,1)  ; 

  %F= G - 0.5 * trace ( (sigma_0*sigma_0^T) * ddu^T ) 
  cross_prod = pagemtimes( sigma0 , pagetranspose (sigma0));
  traccia = sum( pagediag( pagemtimes ( cross_prod , ddumat ) ));
  traccia = reshape(traccia,[L,1]);

  F = G - 0.5 * traccia ;  
end  


function Ft0 = close_formula(S0,T,k,r,alpha,lambda,sigma,rho,delta)

    alphah = alpha - lambda/k;
    AT = ( r - alphah+1/2*(sigma(2)/k)^2 - (sigma(1)*sigma(2)*rho)/k)*T ...
        + 1/4*sigma(2)^2 * (1-exp(-2*k*T))/k^3 +...
        ( alphah*k + (sigma(1)*sigma(2)*rho) -  sigma(2)^2/k )*(1-exp(-k*T))/k^2 ;
    Ft0 = S0 * exp( -delta*(1-exp(-k*T))/k +AT ) ;

end


