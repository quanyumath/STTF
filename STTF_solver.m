function TC = STTF_solver(data,known,Nway,Mtr,opts)
%% our method 
n1=Nway(1);n2=Nway(2);n3=Nway(3);
rank_min = opts.rank_min;maxit=opts.maxIter;tol=opts.tol;
alpha=opts.alpha;
coreNway=opts.coreNway;
[known,id] = sort(known);data = data(id);normdata=norm(data(:));
beta1=opts.beta1;beta2=opts.beta2;beta3=opts.beta3;
%% initialization X_j and Y_j
for j=1:3
    TC{j}=randn(Nway);
    TC{j}(known)=data;
end
for j=3                     
    TC{j}=fft(TC{j},[],j);
    for i = 1:Nway(j)
        if j==1
            [Uk,Sigmak,Vk]=svds(squeeze(TC{j}(i,:,:)),coreNway{j}(i));
        elseif j==2
            [Uk,Sigmak,Vk]=svds(squeeze(TC{j}(:,i,:)),coreNway{j}(i));
        else
            [Uk,Sigmak,Vk]=svds(TC{j}(:,:,i),coreNway{j}(i));
        end
        X{i,j} = Uk*Sigmak;
        Y{i,j} = Vk';
    end
end

  for i = 1:Nway(1)
      X{i,1} = rand(Nway(2),coreNway{1}(i));
      Y{i,1} = rand(coreNway{1}(i),Nway(3));
  end
  for i = 1:Nway(2)
      X{i,2} = rand(Nway(1),coreNway{2}(i));
      Y{i,2} = rand(coreNway{2}(i),Nway(3));
  end
%     for i = 1:Nway(3)
%       X{i,3} = rand(Nway(1),coreNway{3}(i));
%       Y{i,3} = rand(coreNway{3}(i),Nway(2));
%   end

%% Generate \mathcal{H},\mathcal{G},\mathcal{F},
diagT_3=eye(n3,n3-1);  %H
for i=n3:-1:2
    diagT_3(i,:)=diagT_3(i,:)-diagT_3(i-1,:);
end
T3(1,:,:)=diagT_3;
for i=2:n1
    T3(i,:,:)=zeros(size(diagT_3));
end
T3_bar=fft(T3,[],1);
for n=1:n1
    T3sq{n}=squeeze(T3_bar(n,:,:));T3SQ{n}=T3sq{n}*T3sq{n}';
end
%-------------------------------------
diagT_2=eye(n2,n2-1);    %G
for i=n2:-1:2
    diagT_2(i,:)=diagT_2(i,:)-diagT_2(i-1,:);
end
T2(:,:,1)=diagT_2;
for i=2:n3
    T2(:,:,i)=zeros(size(diagT_2));
end
T2_bar=fft(T2,[],3);
for n=1:n3
    T2sq{n}=T2_bar(:,:,n);T2SQ{n}=T2sq{n}*T2sq{n}';
end
%-------------------------------------
diagT_1=eye(n1,n1-1);  %F
for i=n1:-1:2
    diagT_1(i,:)=diagT_1(i,:)-diagT_1(i-1,:);
end
T1(:,1,:)=diagT_1';
for i=2:n2
    T1(:,i,:)=zeros(size(diagT_1'));
end
T1_bar=fft(T1,[],2);
for n=1:n2
    T1sq{n}=squeeze(T1_bar(:,n,:));T1SQ{n}=T1sq{n}'*T1sq{n};
end












%% compute the initialization residual
%C有三种展开形式
for j=1:3
    C{j}=zeros(Nway);
end
for j=1:3                  
    for n = 1:Nway(j)
        if j==1
            C{j}(n,:,:)=X{n,j}*Y{n,j};
        elseif j==2
            C{j}(:,n,:)=X{n,j}*Y{n,j};
        else
            C{j}(:,:,n)=X{n,j}*Y{n,j};
        end
    end
    TTC{j}=ifft(C{j},[],j);
end
TC=alpha(1)*TTC{1}+alpha(2)*TTC{2}+alpha(3)*TTC{3}; 
TC(known) = data;TC=real(TC);
C{3}=fft(TC,[],3);C{2}=fft(TC,[],2);C{1}=fft(TC,[],1);
rho=0.95;mu=1.01;lamda=0.1;
fprintf('Iteration:     ');
tic
for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);
    oldTC=TC;
    %% update (X_j,Y_j)
    for j=1:3                            
        for n = 1:Nway(j)
            if j==3
                Ysq{n,j}=Y{n,j}*Y{n,j}';YT2=Y{n,j}*T2sq{n};
                X{n,j}=alpha(3)*C{j}(:,:,n)*Y{n,j}'*inv(alpha(3)*Ysq{n,j}+beta2*(YT2*YT2')+lamda*eye(size(Ysq{n,j})));     
                Xsq{n,j} = X{n,j}'*X{n,j};
                Y{n,j} = inv(lamda*eye(size(Xsq{n,j}))+Xsq{n,j})*alpha(3)*X{n,j}'*C{j}(:,:,n)*inv(alpha(3)*eye(size(T2SQ{n}))+beta2*T2SQ{n});
                C{j}(:,:,n)=X{n,j}*Y{n,j};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                              
%                 Ysq{n,j}=Y{n,j}*Y{n,j}';
%                 X{n,j}=alpha(3)*C{j}(:,:,n)*Y{n,j}'*inv(alpha(3)*Ysq{n,j}+lamda*eye(size(Ysq{n,j})));     
%                 Xsq{n,j} = X{n,j}'*X{n,j};
%                 Y{n,j} = inv(lamda*eye(size(Xsq{n,j}))+Xsq{n,j})*X{n,j}'*C{j}(:,:,n);
%                 C{j}(:,:,n)=X{n,j}*Y{n,j};
            elseif j==2
                Ysq{n,j}=Y{n,j}*Y{n,j}';                
                X{n,j}=inv(alpha(2)*eye(size(T1SQ{n}))+beta1*T1SQ{n})*alpha(2)*squeeze(C{j}(:,n,:))*Y{n,j}'*inv(Ysq{n,j}+lamda*eye(size(Ysq{n,j})));
                Xsq{n,j} = X{n,j}'*X{n,j};T1X=T1sq{n}*X{n,j};
                Y{n,j} = inv(lamda*eye(size(Xsq{n,j}))+alpha(2)*Xsq{n,j}+beta1*(T1X'*T1X))*alpha(2)*X{n,j}'*squeeze(C{j}(:,n,:));
                C{j}(:,n,:)=X{n,j}*Y{n,j};
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Ysq{n,j}=Y{n,j}*Y{n,j}';
%                 X{n,j}=squeeze(C{j}(:,n,:))*Y{n,j}'*inv(Ysq{n,j}+lamda*eye(size(Ysq{n,j})));
%                 Xsq{n,j} = X{n,j}'*X{n,j};
%                 Y{n,j} = alpha(2)*inv(lamda*eye(size(Xsq{n,j}))+alpha(2)*Xsq{n,j})*X{n,j}'*squeeze(C{j}(:,n,:));
%                 C{j}(:,n,:)=X{n,j}*Y{n,j};
            else
                Ysq{n,j}=Y{n,j}*Y{n,j}';YT3=Y{n,j}*T3sq{n};
                X{n,j}=alpha(1)*squeeze(C{j}(n,:,:))*Y{n,j}'*inv(alpha(1)*Ysq{n,j}+beta3*(YT3*YT3')+lamda*eye(size(Ysq{n,j})));
                Xsq{n,j} = X{n,j}'*X{n,j};
                Y{n,j} = inv(lamda*eye(size(Xsq{n,j}))+Xsq{n,j})*alpha(1)*X{n,j}'*squeeze(C{j}(n,:,:))*inv(alpha(1)*eye(size(T3SQ{n}))+beta3*T3SQ{n});
                C{j}(n,:,:)=X{n,j}*Y{n,j};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Ysq{n,j}=Y{n,j}*Y{n,j}';
%                 X{n,j}=alpha(1)*squeeze(C{j}(n,:,:))*Y{n,j}'*inv(alpha(1)*Ysq{n,j}+lamda*eye(size(Ysq{n,j})));
%                 Xsq{n,j} = X{n,j}'*X{n,j};
%                 Y{n,j} = inv(lamda*eye(size(Xsq{n,j}))+Xsq{n,j})*X{n,j}'*squeeze(C{j}(n,:,:));
%                 C{j}(n,:,:)=X{n,j}*Y{n,j};
            end
        end
    end
    %% adjust the rank of (X_j,Y_j)
    for j=1:3                                         
        if  rho<1
        max_k=max(coreNway{j});
        sum_k=sum(coreNway{j});
        sigmas=zeros(max_k*Nway(j),1);
        for i=1:Nway(j)
            s = svd(Xsq{i,j});
            sigmas((i-1)*max_k+1:(i-1)*max_k+length(s))=s;
        end
        [dR,id]=sort(sigmas,'descend');
        drops = dR(1:sum_k-1)./dR(2:sum_k);
        [dmx,imx] = max(drops);
        rel_drp = (sum_k-1)*dmx/(sum(drops)-dmx);
        if rel_drp>10
            thold=rho*sum(dR);
            iidx=0;ss=0;
            len=length(dR);
            for i=1:len
                ss=ss+dR(i);
                if(ss>thold)
                    iidx=i;
                    break;
                end
            end
            if(iidx>sum(rank_min{j}(length(Nway(j)))))  %length(N(j))
                idx=floor((id(iidx+1:sum_k)-1)/max_k);
                for n=1:Nway(j)
                    num=length(find(idx==n-1));
                    if(num>0)
                        if coreNway{j}(n)-num>rank_min{j}(n)
                            coreNway{j}(n) = coreNway{j}(n)-num;
                        else
                            coreNway{j}(n) = rank_min{j}(n);
                        end
                        [Qx,Rx] = qr(X{n,j},0);
                        [Qy,Ry] = qr(Y{n,j}',0);
                        [U,S,V] = svd(Rx*Ry');
                        sigv = diag(S);
                        X{n,j} = Qx*U(:,1:coreNway{j}(n))*spdiags(sigv(1:coreNway{j}(n)),0,coreNway{j}(n),coreNway{j}(n));
                        Y{n,j} = (Qy*V(:,1:coreNway{j}(n)))';
                        if j==3
                            C{j}(:,:,n)=X{n,j}*Y{n,j};
                        elseif j==2
                            C{j}(:,n,:)=X{n,j}*Y{n,j};
                        else
                            C{j}(n,:,:)=X{n,j}*Y{n,j};
                        end    
                    end
                end
            end
            rho=rho*mu;
        end
        end
    end
  %% update TC  
    TC=alpha(1)*ifft(C{1},[],1)+alpha(2)*ifft(C{2},[],2)+alpha(3)*ifft(C{3},[],3);   
    %% judge whether converges 
    ress=norm(TC(:)-oldTC(:))/norm(oldTC(:));
    if  ress<tol || k==maxit
        TC(known) = data;
        TC=real(TC);
        TC = max(TC,0);TC=min(TC,1);
        break;
    end
    TC(known) = data;TC=real(TC);TC = max(TC,0);TC=min(TC,1);
%    RSE(k)=norm(TC(:)-Mtr(:))/norm(Mtr(:))
    TC=1.2*TC-0.2*oldTC;
    %% update C
    C{3}=fft(TC,[],3);C{2}=fft(TC,[],2);C{1}=fft(TC,[],1);
end
end