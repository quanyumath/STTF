function STTF=STTF_video( X,Omega)

%% produce data

data=X(Omega);
known=Omega;
%% our method 
opts = [];
opts.maxIter=300;
opts.tol = 1e-5; % run to maxit by using negative tolerance
opts.Mtr = X; % pass the true tensor to calculate the fitting
Nway=size(X);

opts.alpha = [0.1,0.1,1];opts.alpha = opts.alpha/sum(opts.alpha);
opts.rank_min = {1*ones(1,Nway(1)),1*ones(1,Nway(2)),1*ones(1,Nway(3))};
opts.coreNway ={3*ones(1,Nway(1)),3*ones(1,Nway(2)),50*ones(1,Nway(3))};
opts.beta1=0.001;opts.beta2=0.01;opts.beta3=20;
STTF=STTF_solver(data,known,Nway,X,opts);

end



