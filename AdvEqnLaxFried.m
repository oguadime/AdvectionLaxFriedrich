function AdvEqnLaxFried=AdvEqnLaxFried(M,nu)

format long 

%% set up initial data
A = 0; 
B = 1;
dx = (B-A)/M;
x = (A:dx:B); uprev=uinit(x); unew = 0*uprev; t = 0; n = 0; a = 1;
dt = nu*dx/a;
%%
plot(x,uprev); 
title(sprintf('Initial condition at t=%g',t)); 
pause;

while t < 0.5
    %% advance to new time step
    t = t + dt;
    exvec = exfun(x,t);
    unew(1) = exfun(0,t); 
    for j = 2:M
         unew(j) = 0.5*(uprev(j-1)+uprev(j+1))-0.5*nu*(uprev(j+1)-uprev(j-1));
    end
    unew(M+1) = 0.5*(uprev(M)+uprev(2))-0.5*nu*(uprev(2)-uprev(M)); %% Periodic B.C
    
    %% plot, compare with true solution ...?
    plot(x,unew,'bo-',x,exvec,'r'); 
    title(sprintf('Numerical solution at t=%g',t));
    axis([-1 3 -0.1 1.5]);
    legend('numerical','exact');
    pause;

    %% calculate error 
    uprev = unew;
    gridL2 = sqrt(dx)*norm(exvec-unew,2);
    gridL1 = dx*norm(exvec-unew,1);
    gridLinf = norm(exvec-unew,inf);
end

fprintf('Grid Error L2=%g\n',gridL2);
fprintf('Grid Error L1=%g\n',gridL1);
fprintf('Grid Error Linf=%g\n',gridLinf);
end

%% also need functions 


function v = uinit(x)
    v = 0*x;
    for j = 1:length(x)
    if (x(j)>0) && (x(j)<1)
        v(j)=exp(-100*(x(j)-0.3).^2);
    else 
        newvalue=x(j)-(floor(x(j)));
        v(j)=exp(-100*(newvalue-0.3).^2);
    end
    end
end


function v = exfun(x,t)
    v=uinit(x-t);
end

