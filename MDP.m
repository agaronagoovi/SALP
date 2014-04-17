function [ V, w ] = MDP( Pssa,L )
%MDP Summary of this function goes here
%   Detailed explanation goes here
mqlength = 5;
nu = 3;
nx = (mqlength+1)^3;
sz = mqlength+1;
surrogateV = @(q) norm(q)^2;
V = zeros(nx,1);
alpha = 0.98;


for q1=0:mqlength
    for q2=0:mqlength
        for q3=0:mqlength
            i = addr(q1+1,q2+1,q3+1,mqlength);
            V(i) = surrogateV([q1,q2,q3]);
        end
    end
end
    
% compute Hamiltonian for current v
H = zeros(nx,nu);
for iu = 1:nu
    H(:,iu) = L + alpha*Pssa(:,:,iu)*V;
end

% compute policy
[~, policy] = min(H,[],2);

init_dstrbt = ones(nx,1)./nx;

P = zeros(nx,nx);
for i = 1:nx
    P(i,:) = Pssa(i,:,policy(i));
end

stat_dstrbt = (1-alpha)*init_dstrbt'/(eye(nx)-alpha*P);
S = zeros(200,2);
for i=1:size(S,2)
    S(:,i) = randsample(nx,size(S,1),true,stat_dstrbt);
end

%%basis matrix%%
basis = @(q) [1, q(1)^2, q(2)^2, q(3)^2];
Phi = zeros(nx,4);
for x1=1:sz
    for x2=1:sz
        for x3=1:sz
            i = addr(x1,x2,x3,mqlength);
            Phi(i,:) = basis([x1-1,x2-1,x3-1]);
        end
    end
end


implicit = 1;
w = cell(size(S,2),2);
options = optimoptions('linprog','Display','iter');
weights = zeros(size(Phi,2),1);
%%Linear Program explicit budget%%
if ~implicit
    theta = 1;
    for e = 1:size(S,2)
        uS = unique(S(:,e));
        w{e,1} = uS;
        f=zeros(size(Phi,2)+size(uS,1),1);
        for i=1:size(uS,1)
            f = f - [Phi(uS(i,1),:)';zeros(size(uS,1),1)];
        end
        f = f./size(uS,1);
        
        A = zeros(size(uS,1)*nu,size(Phi,2)+size(uS,1));
        b = zeros(size(uS,1)*nu, 1);
        for iu = 1:nu
            starti = (iu-1)*size(uS,1) + 1;
            endi = iu*size(uS,1);
            T = Phi-alpha*Pssa(:,:,iu)*Phi;
            s=1;
            for i=starti:endi
                A(i,1:size(Phi,2)) = T(uS(s,1),:);
                A(i, size(Phi,2)+s) = 1;
                b(i,1) = L(uS(s,1));
                s = s+1;
            end
        end
        A = [A;[zeros(1,size(Phi,2)),ones(1,size(uS,1))]];
        b = [b;size(uS,1)*theta];
        lb = [-Inf*ones(size(Phi,2),1);zeros(size(uS,1),1)];
        w{e,2} = linprog(f,A,b,[],[],lb,[],[],options);
        weights = weights + w{e,2}(1:size(Phi,2));
    end


%%Linear Program implicit budget%%
else
    for e = 1:size(S,2)
        uS = unique(S(:,e));
        w{e,1} = uS;
        f=zeros(size(Phi,2)+size(uS,1),1);
        for i=1:size(uS,1)
            f(1:size(Phi,2)) = f(1:size(Phi,2)) - Phi(uS(i,1),:)';
        end
        f(size(Phi,2)+1:end) = ones(size(uS,1),1)*(2/(1-alpha));
        f = f./size(uS,1);

        A = zeros(size(uS,1)*nu,size(Phi,2)+size(uS,1));
        b = zeros(size(uS,1)*nu, 1);
        for iu = 1:nu
            starti = (iu-1)*size(uS,1) + 1;
            endi = iu*size(uS,1);
            T = Phi-alpha*Pssa(:,:,iu)*Phi;
            s=1;
            for i=starti:endi
                A(i,1:size(Phi,2)) = T(uS(s,1),:);
                A(i, size(Phi,2)+s) = 1;
                b(i,1) = L(uS(s,1));
                s = s+1;
            end
        end
        lb = [-Inf*ones(size(Phi,2),1);zeros(size(uS,1),1)];
        w{e,2} = linprog(f,A,b,[],[],lb,[],[],options);
        weights = weights + w{e,2}(1:size(Phi,2));
    end      
end
weights = weights./size(S,2);

V = Phi*weights;

    
end




