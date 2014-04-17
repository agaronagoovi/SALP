function [exp_cost, trajs, policy] = evalPolicy(V, Pssa, L,policy)
%EVALPOLICY Summary of this function goes here
%   Detailed explanation goes here
mqlength = 5;
nu = 3;
nx = (mqlength+1)^3;
alpha = 0.98;
if nargin<4
    H = zeros(nx,nu);
    for iu = 1:nu
        H(:,iu) = L + alpha*Pssa(:,:,iu)*V;
    end
    % compute policy
    [~, policy] = min(H,[],2);
end


START = 1;
s = START;
MAXSTEPS = 5000;
MAXITER = 10;
trajs = cell(MAXITER,2);
exp_cost = 0;

for e=1:MAXITER
    traj = zeros(MAXSTEPS+1,1);
    cost = 0;
    for i=1:MAXSTEPS
        traj(i) = s;
        u = policy(s);
        p = Pssa(s,:,u);
        ns = randsample(nx,1,true,p);
        cost= cost + (alpha^(i-1))*L(ns);
        s = ns;
    end
    traj(MAXSTEPS+1) = ns;
    trajs{e,1} = traj;
    %cost = cost/MAXSTEPS;
    exp_cost = exp_cost + cost;
    trajs{e,2} = cost;
end

exp_cost = exp_cost/MAXITER;

end

