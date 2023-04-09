

function param = fmincon_svm(yt)

    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-20);
    options  =  optimset(options , 'TolX'        , 1e-10);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'off');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'Maxiter'      , 50000);
    options  =  optimset(options , 'MaxFunEvals' , 400*6);
    options = optimset(options, 'HessUpdate', 'bfgs');
    options = optimset(options, 'DiffMinChange', 0.00001);
    options = optimset(options, 'DiffMaxChange', 0.001);
    options = optimset(options, 'UseParallel', false);
    
    A = [0; 
        -0.0535;
        0.2940*2.7471; 
        0.1024*0.2829];

    B = [-1.1374 1.1374 -1 2.1302;
         0  0   0           0;
         0  0   -2.7471     0;
         0  0   0       -0.2829];

    omega1 = [1.0000 0.84000 -0.0405 0.02300;
              0.8400 1.1865  -0.3350 -0.2148;
              -0.0405 -0.335 0.4411  0.1733;
              0.02300 -0.2148 0.1733 0.1227];

    
    theta = [nonzeros(A); nonzeros(B); nonzeros(omega1)];


 
    param = fminunc( 'likelihood_kalman', theta , options, yt,T,nint);
  


%     A = []; %beta + alpha < 1
%     b = [];
%     Aeq=[];
%     beq=[];
%     lb =[];
%     ub= [];
%     nonlcon = [];
%     fun = @(p0) likelihood_kalman(p0,yt,T,nint);
% 
%     param=fmincon(fun,p0, A,b,Aeq,beq,lb,ub,nonlcon,options);
    
end