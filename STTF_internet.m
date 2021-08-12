function STTF=STTF_internet( X,Omega)

%% produce data

data=X(Omega);
known=Omega;
%% our method 
opts = [];
opts.maxIter=150;
opts.tol = 1e-5; % run to maxit by using negative tolerance
opts.Mtr = X; % pass the true tensor to calculate the fitting
Nway=size(X);

opts.alpha = [1,1,1];opts.alpha = opts.alpha/sum(opts.alpha);
opts.rank_min = {20*ones(1,Nway(1)),7*ones(1,Nway(2)),7*ones(1,Nway(3))};
opts.coreNway ={60*ones(1,Nway(1)),7*ones(1,Nway(2)),7*ones(1,Nway(3))};
%[2 2 0.8 0.6 0.4 5e-02 5e-02 5e-02 5e-02];
opts.beta1=0;opts.beta2=2;opts.beta3=0;
STTF=STTF_solver(data,known,Nway,X,opts);

end