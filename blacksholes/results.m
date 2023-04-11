
xaxis=5:5:60;
sols=[]; sols2=[];
%do the variance of solutions
for n=xaxis
    L=n^2*600;
    num=bs_trinomial_constant( n );
    num2=bs_reg_1d( n , L );
    sols=[sols,num];

    sols2=[sols2,num2];
    n
end
plot(xaxis,sols,'r');yline( 0 , 'g' )
hold on 
plot(xaxis,sols2,'b') 
hold off





for n=1:20
    num=bs2d( 20 , 1e6 );
    sols=[sols,num];
    n
end

var(sols)
