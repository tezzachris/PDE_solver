function [ycall,yput] = bs_price1d(r,sig,qdiv,T,t0,E,S0)

% CALCOLA IL PREZZO DI UNA CALL E DI UNA PUT VANIGLIA
% r: TASSO DI INTERESSE
% sig: VOLATILITA'
% qdiv: RATEO DI DIVIDENDO
% T: SCADENZA
% t0: DATA DI VALUTAZIONE (INIZIALE, A PATTO DI DEFINIRE OPPORTUNAMENTE T PUO' ESSERE SETTATA A 0)
% E: STRIKE
% S0: PREZZO DEL SOTTOSTANTE AL TEMPO t0 

d1 = (log(S0/E) + (r - qdiv + 0.5d0*sig*sig)*(T-t0))/(sig*sqrt(T-t0));
d2 = d1-(sig*sqrt(T-t0));

ycall = S0.*normcdf(d1)*exp(-qdiv*(T-t0)) - E*(exp(-r*(T-t0)))*normcdf(d2);

%yput = ycall + E*exp(-r*(T-t0)) - S0*exp(-qdiv*(T-t0));

% la riga sopra calcola il prezzo della put usando la cosiddetta "Put-Call
% parity", ultimo argomento del corso.
% In alternativa, si poteva anche usare la formula di Black-Scholes, come abbiamo già
% visto a lezione:

yput =  E*(exp(-r*(T-t0)))*normcdf(-d2) - S0.*normcdf(-d1)*exp(-qdiv*(T-t0));


end