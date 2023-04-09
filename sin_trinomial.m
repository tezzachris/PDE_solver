%compute the backward solution 
%example 6 paper figlio
%one dimension
function [sol,sol_exact] =sin_trinomial( xo, n , T, d, sigma_up, sigma_down ) 
  h = T/n ;
  v_rovescia = sigma_up^2 / sigma_down^2; 
  p = min( 1/(2*(v_rovescia-1)) , 1/3); %prob of epsilon (if 1/3 we have trinomial)
  alfa_p = 2 ;%fixed in paper
  alfa_down = 1/(2*p*v_rovescia + alfa_p);
  sigma_0 = ( sigma_down / sqrt(2*alfa_down) ) * eye(d); %interval width
  x=xfin(xo,n,sigma_0,p,h); %final vector of trinomial forward process
  y = zeros( n+1 , 2*n+1 ) ; 
  y(1,:)= arrayfun(@(x) sin(T+x),x) ; %terminal cond example 2.6.1
  f = @(u,z,d,sigma_up,sigma_down) (1/d * z - d/2 * condition(sigma_down,sigma_up,u)^2 * u);
  G = @(u,z,zz,d,sigma_up,sigma_down) ( 1/2 * trace(condition(sigma_down,sigma_up,zz)^2*eye(d,d) *zz) - f(u,z,d,sigma_up,sigma_down)); %parabolic if nondecreasing in zz
  F = @(u,z,zz,d,sigma_up,sigma_down,sigma_0) (G(u,z,zz,d,sigma_up,sigma_down) - 1/2 * trace((sigma_0.^2) * zz));
 for i = 2:(n+1)
      cont=0;
      for k = transpose(triplets(y(i-1,1:end-2*(i-2))))%transpose(triplets(transpose(nonzeros(y(i-1,:)))))
          cont=cont+1; %column selector
          y(i,cont)=expectation(transpose(k),p,0,h,sigma_0,d)+h*F(expectation(transpose(k),p,0,h,sigma_0,d),expectation(transpose(k),p,1,h,sigma_0,d),expectation(transpose(k),p,2,h,sigma_0,d),d,sigma_up,sigma_down,sigma_0) ;
      end
  
  end
  sol=nonzeros(y(n+1,:));
  sol_exact=sin(0+xo);
end

%compute ending values of trinomial tree

function [y]=xfin(xo,n,sigma_0,p,h)
    y=flip([ flip(xo + (sqrt(h)*sigma_0*(-1/sqrt(p)))*[1:n]),xo,xo + (sqrt(h)*sigma_0*(1/sqrt(p)))*[1:n] ]);
end

%from a row vector (of odd size) construct a matrix of collection of 1x3
%elements

function ye =triplets(y) %y vettore riga
    ye=zeros((size(y,2)-2),3);
    for i=1:(size(y,2)-2)
        ye(i,:)=y(1,i:i+2);
    end
end

function e = expectation(y,p,k,h,sigma_0,d) %y 1x3
    eps=[1/sqrt(p),0, -1/sqrt(p)]; 
    eps_prob=[p/2,(1-p),p/2]; %epsilon trinomial variables
    ke=eye(3,3); %diagonal matrix for multiply with y
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

%condition to verify the sup and inf of the function

function y=condition(inf,sup,u) %check function for sup and inf
if u>=0
    y=sup;
elseif u<0
    y=inf;
end

end 