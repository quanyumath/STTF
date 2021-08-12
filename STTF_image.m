function STTF=STTF_image( X,Omega)

%% produce data

data=X(Omega);
known=Omega;
%% our method 
opts = [];
opts.maxIter=150;
opts.tol = 1e-5; % run to maxit by using negative tolerance
opts.Mtr = X; % pass the true tensor to calculate the fitting
Nway=size(X);

%% 
opts.alpha = [1,2,3]/6;
opts.rank_min = {1*ones(1,Nway(1)),1*ones(1,Nway(2)),5*ones(1,Nway(3))};
opts.coreNway ={3*ones(1,Nway(1)),3*ones(1,Nway(2)),50*ones(1,Nway(3))};
% p=0.1
 opts.beta1=1;opts.beta2=0.5;opts.beta3=0;
% p=0.2
%opts.beta1=0.5;opts.beta2=0.1;opts.beta3=0;
% p=0.3
%opts.beta1=0.1;opts.beta2=0.05;opts.beta3=0;
STTF=STTF_solver(data,known,Nway,X,opts);

end



