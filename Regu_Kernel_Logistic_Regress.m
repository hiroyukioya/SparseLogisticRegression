function [alpha, beta, pp, w, Z, lik, ac, K ] = Regu_Kernel_Logistic_Regress(X, T, kernel, ksigma, lambda)
%%
%  /////////   Regularized Kernel Logistic Regression    /////////
%
%  Refs: M. Maalouf, et. al.,     Comput. Stat. Data. Anal (2011)
%            Cawley and Talbot    Mach. Learn. (2008) 71; 243-264    
%            Zhu and Hastie         J. Comp. Graph. Stat. (2005), 14, 185-205 
%             
%
%  [ INPUTS ]
%     X = input data [ n x p ]  (n=observations, p=dimension)
%     T = target (response) vector [ n x1 ] (Binary)
%     kernel = "linear", "rbf" or "polynomial"
%     ksigma= kernel parameter for non-linear kernel 
%         if kernel="polynomial" , specify scaling parameter ksigma(2) .
%                         e.x., ksigma=[3 2] ..... degree of poly=3
%                         scaling .........   paramete=2
%     lambda= inverse of regularization parameter (scalar)
%                 small lambda --> strong regularization 
%
%  [OUTPUTS]
%     alpha, beta:  Coefficients
%     pp: LOO estimate of posterior probability
%     w:  Primal space weights
%     Z: LOO outputs(eta)
%     lik:  Negative log likelihood
%     p: trained data p

%                                                 H.Oya (2010)
%%   Initialize data matrices...

T=T(:);          %T is  a target vector
[n,p]=size(X);
[c,d]=size(T);

if n~=c;
    error('X and Y matrices must have the same number of rows.....')
end

% Select kernel ...............
o=strncmpi(kernel,'rbf',3);
o3=strncmpi(kernel,'gauss',4);
o1=strncmpi(kernel,'polynomial',4);
o2=strncmpi(kernel,'linear',3);


% Construct Kernel Matrix --------------
if o2==1      % Liner kernel
    mode=1;
    K=X*X';
elseif o==1 | o3==1     % RBF kernel
    mode=2;
    [K]=KernelMatrix(X,ksigma,'gaus');
elseif o1==1         % Polynomial kernel
    mode=3;
    [K]=KernelMatrix(X,ksigma,'poly');
end

%% MAIN LOOP
maxiter=150;
alpha=zeros(n,1);    % [n x 1]
eta=zeros(n,1);
eta=randn(n,1);         % [n x 1]
mu=exp(eta)./(1+exp(eta));
y= exp(eta)./(1+exp(eta)).^2;    % = p*(1-p)
old_lik= 10^8;
tol=10^-6;

%--------- Cowley's  loop (Efficient IRWLS) ---------
for i=1:maxiter

    % Weighting matrix
    W=y;      % [n x 1]
    IW = diag(1./W);    % [n x n]
    
    % Adjusted psudo-response
    z = eta+(T-mu)./W;    % [n x 1]

    % Construct M
    M=K+lambda*IW;    % [ n x n]
    old_M=M;
    % Cholesky factorization
    R=chol(M) ;
    % Solve first equation
    zeta=R\(R'\ones(n,1));     % [n x n]
    % Solve second equation
    xi=R\(R'\z);
    % Calculate beta and alpha
    new_beta=sum(xi)/sum(zeta);
    new_alpha=xi-zeta*new_beta;
    
    % Compute eta
    eta = K*new_alpha+new_beta;
    % Compute mu
    mu = exp(eta)./(1+exp(eta));
    % Compute y , this is needed for weighting (W)......
    % Below is equal to y*(1-y)
    y= exp(eta)./(1+exp(eta)).^2;
    
    % Negative log-likelihood (cross-entropy)
    e1 = -(T.*eta-log(1+exp(eta)));
    lik=(0.5*lambda*alpha'*K*alpha)+sum(e1);
    
    % Evaluate the stopping criterion
     D=(old_lik-lik);
%       D=(alpha-new_alpha);
%       figure(1);plot(i,lik,'s');hold on
    
    if D<tol & i<=maxiter
%           fprintf(' Iterarion for convergence   %2.0f :  ',i-1);
          alpha=new_alpha;
          beta=new_beta;
          lik=old_lik;
          p=old_y;
          M=old_M;
        break
    elseif D>=tol & i>=maxiter
        error(' Iteration did not converge .....')
    end
    old_y=y;
    old_lik=lik;
    alpha=new_alpha;  % Update alpha (Dual-space weights)
    i=i+1;
end

% Feature-space weights (Valid for linear kernel case) 
w = alpha'*X;

%% LOO cross validation

%-------------------- 
iR = inv(R);
cii=sum(iR.^2,2)-zeta.^2/sum(zeta);
% LOO outputs
Z=z-alpha./cii;
% LOO post. prob.
pp=exp(Z)./(1+exp(Z));   

% LOO regularized negative log-likelihood (cross-entropy).
%          <<<--- log(1+exp(-tz))=log(exp(1+tz))-tz;
e1 = -(T.*Z-log(1+exp(Z)));
lik=(0.5*lambda*alpha'*K*alpha)+sum(e1)/2;
% 
%     f3=find(T==0);   f4=find(T==1);    
%     ttt=length(f4)/(length(T)); % ttt is the decision threshold
%     f1=find(pp<=ttt);    
%     f2=find(pp>ttt);
%     % LOO accuracy..........
%     ac=(length(intersect(f1,f3))+length(intersect(f2,f4)))/length(T); 

ip=sign(pp-0.5);
tt=T; f=find(tt==0); tt(f)=-1;
L=find(tt.*ip==1);
%LOO accuracy
ac=length(L)/length(T);
fprintf(' LOO accuracy   %2.4f \n',ac);


%%
function [K,D,CK]=KernelMatrix(data,sigma,kernel)
%%  Kernel matrix 
%  Compute kernel matrix using Gaussian kernel with sigma value
 
if sigma<=0;
    error('       sigma must be positive     ');
end

[a,b]=size(data);
K=zeros(a,a);

o=strncmpi(kernel,'gaussian',4);
o1=strncmpi(kernel,'polynomial',4);

if (o==1) & (o1==0)
    for n=1:a
            for m=n:a
                  [K(n,m)]=gaussian_kernel_function(data,n,m,sigma);
            end
    end
elseif (o1==1) & (o==0)
    for n=1:a
            for m=n:a
                  [K(n,m)]=polynomial_kernel_function(data,n,m,sigma);
            end
    end
end
    
% Kernel matrix
for n=1:a
      m=n:a;
      K(m,n)=K(n,m)';
end

% Normalize
D=diag(1./sqrt(diag(K)));
NK=D*K*D;

%Centering
d=sum(K)/a;
e=sum(d)/a;
J=ones(a,1)*d;
CK=K-J-J'+e*ones(a,a);

%%
function [a]=gaussian_kernel_function(data,ind1,ind2,sigma)

k1=dot(data(ind1,:),data(ind1,:));
k2=dot(data(ind2,:),data(ind2,:));
kb=dot(data(ind1,:),data(ind2,:));

% Gaussian Kernel Function
a=exp(-(k1+k2-2*kb)/(2*sigma(1)^2));

%%
function [a]=polynomial_kernel_function(data,ind1,ind2,sigma)
eta=1/(sigma(2)^2); 
kb=data(ind1,:)*eta*data(ind2,:)';
% inhomogenous polynomial kernel of degree sigma
a=(kb+1).^sigma(1);
