
function [CallMin] = JohnsonMinLuca(S0vet,sigmavet,r,E,T,rhomat)

d=length(S0vet);
qdivvet=zeros(d,1);
N = max(size(rhomat));

for i=1:1:N
    m(i) = r-qdivvet(i)-0.5*sigmavet(i)^2;
end

for i = 1:1:N
    for j=1:1:N
        sig2(i,j) = sigmavet(i)^2 + sigmavet(j)^2 -2*rhomat(i,j)*sigmavet(i)*sigmavet(j);
        TermToAdd = (qdivvet(i)-qdivvet(j))*T;
        d1primo(i,j) = (log(S0vet(i)/S0vet(j))+0.5*sig2(i,j)*T-TermToAdd)/(sqrt(sig2(i,j)*T));
    end
end

sig2square = sqrt(sig2);

for i = 1:1:N
    for j=1:1:N
        for k=1:1:N
            % i = 3; j = 2; k=1
            if (j==k)
                rho(i,j,k) = 1;
            else
                if (i==j)
                    rho(i,j,k) = (sigmavet(i)-rhomat(i,k)*sigmavet(k))/sig2square(i,k);
                else
                    if (i==k)
                        rho(i,j,k) = (sigmavet(i)-rhomat(i,j)*sigmavet(j))/sig2square(i,j);
                    else
                        num = sigmavet(i)^2 -rhomat(i,j)*sigmavet(i)*sigmavet(j)- ...
                        rhomat(i,k)*sigmavet(i)*sigmavet(k)+rhomat(j,k)*sigmavet(j)*sigmavet(k);
                        den = sig2square(i,j)*sig2square(i,k);
                        rho(i,j,k) = num/den;
                    end;
                end;
            end;
        end;
    end;
end;

rhopermutata = permute(rho,[3,2,1]);

for i=1:1:N
    d2(i) = (log(S0vet(i)/E) + m(i)*T)/(sigmavet(i)*sqrt(T));
    d1(i) = d2(i) + sigmavet(i)*sqrt(T);
end

for i=1:1:N
    argd(i,1) = d1(i);
    cont = 1;
    aggiungid = zeros(1,N-1);
    %aggiungirho = zeros(1,N-1);
    for j=1:1:N
        if (j==i)
            %cont = cont + 1;
        else
            aggiungid(cont) = d1primo(i,j);
     %       aggiungirho(cont) = rho(i,j,i);
            cont = cont + 1;
        end;
    end
    
    argd(i,2:N) = aggiungid;
    % argrho = aggiungirho;
    
  %  enne(i) = mvncdf(argd(i,1:N),0*argd(i,1:N),rhopermutata(:,:,i));
  enne(i) = 0;  
    
end;





%%%%%% RICALCOLA RHO TENENDO CONTO CHE rho_i,j,k è la correlazione tra 
% log(S_j/S_i) e log(S_k,S_i) se k e j sono diversi da i
% se k è uguali a i allora rho_i,j,k è la correlazione tra log(S_i) e 
% log(S_j/S_i) 
% se j è uguali a i allora rho_i,j,k è la correlazione tra log(S_i) e 
% log(S_k/S_i)

for i=1:1:N
    for j=1:1:N
        for k=1:1:N
            if (i==1) && (j==1) && (k==3)
                k=k;
            end;
            if (j==k)
                rho2(i,j,k) = 1;
            else
            if (j==i)
                rho2(i,j,k) = (sigmavet(i)-rhomat(i,k)*sigmavet(k))/sig2square(i,k);
            else
                if (k==i)
                    rho2(i,j,k) = (sigmavet(i)-rhomat(i,j)*sigmavet(j))/sig2square(i,j);
                else
                num = sigmavet(i)^2 -rhomat(i,j)*sigmavet(i)*sigmavet(j)- ...
                rhomat(i,k)*sigmavet(i)*sigmavet(k)+rhomat(j,k)*sigmavet(j)*sigmavet(k);
                        den = sig2square(i,j)*sig2square(i,k);
                        rho2(i,j,k) = num/den;
                end;
            end;
            end;
        end;
    end;
end;

    



% rho2 è tale per cui, per i fissato:
% rho2(i,j,j) = 1
% rho2(i,i,j) = corr(log(S_i), log(S_i/S_j))
% se i diverso da k, rho2(i,j,k) = corr(log(S_i/S_j), log(S_i/S_k))

% rho2permutata pone le correlazioni in una matrice gestibile da matlab,
% ovvero rho2permutata(:,:,i) è la matrice di correlazione che serve per i fissato
% nella formula. Occorre permutare altrimenti matlab gestisce male
% rho2(i,:,:), ovvero non lo considera una matrice 

%%%%%%%%% inverto rho2
%for i = 1:1:N
%   for j=1:1:N
%       for k=1:1:N
%           if (j == k)
%           else 
%               if (j == i) || (k == i) 
%                       rho2(i,j,k) = -rho2(i,j,k);
%               end;
%           end;
%       end;
%   end;
%end;

for i = 1:1:N
   for j=1:1:N
       for k=1:1:N
           if (j == k)
           else 
               if (j == i) || (k == i) 
                       rho2(i,j,k) = -rho2(i,j,k);
               end;
           end;
       end;
   end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho2permutata = permute(rho2,[3,2,1]);

% ricalcola argd
for i=1:1:N
    for j=1:1:N
    if (i==j);
        argd(i,j) = d1(i);
    else
        % ho aggiunto il meno
    argd(i,j) = -d1primo(i,j);
    end;
    end;
end;

zeroN = zeros(1,N);
% CALCOLA GLI N
for i=1:1:N
    rr = rho2permutata(:,:,i);
    vettore = argd(i,:);
    enne(i) = mvncdf(vettore,zeroN,rr);
    %if (i==2)
    %    norm(argd(i,:));
    %    norm(rr);
    %    enne(i);
    %end
end;

% sotto ho messo il segno piu' dove prima era -d2
ennex = mvncdf(d2,0*d2,rhomat);


% sotto invece di "1-enne" ho messo "enne"
CallMin = sum(S0vet'.*enne.*exp(-qdivvet'*T)) -E*exp(-r*T)*(ennex);
end
  