function   [rambda, lik, ac ] = findparameter_RKLR(X, T, kernel, ksigma,range)

v=10.^[range(1):0.2:range(end)];
for n=1:length(v)
    
     [alpha, beta, pp, w, Z, lik(n), ac(n), K ] = ...
                     Regu_Kernel_Logistic_Regress(X, T, kernel, ksigma, v(n));
end
subplot(2,1,1); plot(v,ac,'.-');set(gca,'xscale','log')
subplot(2,1,2); plot(v,lik,'.-');set(gca,'xscale','log')

[f,m]=max(ac);
rambda=v(m);
