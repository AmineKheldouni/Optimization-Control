funcprot(0);
// Initialization
Tf = 10 ;
Nx = 9 ; // x#
xref = 5;
state = [1:Nx];
// Initialization of the matrix used to store Bellman values
V = ones(Tf,Nx) * %inf;
U = ones(Tf-1,Nx) * %inf;


// Compute B(x) subset :
function [res]=B(x)
    if x == 1 then
        res = [1];
    elseif x == Nx then
        res = [-1];
    elseif x == 2 then
        res = [0,1];
    elseif x == Nx-1 then
        res = [-1,0];
    else
        res = [-1,0,1];
    end
endfunction

// Compute Bellman function at final time
function [res]=K(x)
    res = (x-xref).^2;
endfunction

function [res]=p(w)
    if w == 1 then
        res = 1/2;

    else
        res = 1/2;
    end     
endfunction


X = zeros(1,Nx);
for i=1:Nx
    X(1,i) = i;
end
V(Tf,:) = K(X);

// make a time backward loop
Ti = 1;

for t = Tf-1:-1:Ti do
    for x = 1:Nx do
        cost_x = %inf;
        u_opt = 0;
        Bset = B(x);
        for k = 1:length(Bset) do
           u = Bset(k);
            cost_x_u_w = 0;
            for w = -1:2:1 do
                cost_x_u_w = cost_x_u_w + V(t+1,x+u+w)*p(w);
            end
            cost_x_u = cost_x_u_w;
            if cost_x_u < cost_x then
                u_opt = u;
                cost_x = cost_x_u;
            end
        end
        V(t,x) = cost_x;
        U(t,x) = u_opt;
    end
end

disp(V);
disp(U);


function [ud]=optControl(t,x)
    cost_x = %inf;
    ud = 0;
    Bset = B(x);
    for k = 1:length(Bset) do
       u = Bset(k);
        cost_x_u_w = 0;
        for w = -1:2:1 do
            cost_x_u_w = cost_x_u_w + V(t+1,x+u+w)*p(w);
        end
        cost_x_u = cost_x_u_w;
        if cost_x_u < cost_x then
            ud = u;
            cost_x = cost_x_u;
        end
    end
endfunction


ud = optControl(3,5);
disp(ud);


function [cost]=simulation_mc(x0,policy,N)
    cost = 0;
    W = [1,-1];
    for i = 1:N
     x = x0;
        for t = 1:Tf-1
            w = W(grand(1,1,'uin',1,2));
            x = x + policy(t,x) + w;
        end
        cost = cost + (x-xref).^2
    end
    cost = cost/N;
endfunction

timer()
for x0=1:Nx do
    disp(simulation_mc(x0,optControl,10000));
end
timer()

for x0=1:Nx do
    disp(simulation_mc(x0,optControl,1000));
end
function [cost]=simulation_ex(x0,policy)
    // Exact computation with the law of W
    Wa=all_w(Tf-1);
    cost = 0;
    for i = 1: size(Wa,'r')
        x = x0;
        for t = 1:Tf-1
            x = x + policy(t,x) + Wa(i,t);
        end
        cost = cost + (x-xref).^2
    end
    cost = cost/size(Wa,'r');
endfunction

function W=all_w(n)
    // generated all the possible (W_1,...,W_(TF-1))
    if n == 1 then W=[-1;1]
        else
            Wn = all_w(n-1);
            W=[-1*ones(size(Wn,'r'),1),Wn;
        +1*ones(size(Wn,'r'),1),Wn];
    end
endfunction;

timer()
for x0=1:Nx do
    disp(simulation_ex(x0,optControl));
end
timer()

function costs = simulation_dp(policy)
    // evaluation by dynamic programming with fixed policy
    Vs = ones(Tf,length(X))*%inf;
    // Bellman function at time TF
    Vs(Tf,:) = (X-xref).^2;
    W = [-1,1];
    // Compute final value functions
    // Loop backward over time:
    for t = (Tf-1):-1:1
        for x= 1:Nx
            // loop on noises
            EV=0;
            for iw = 1:size(W,"*")
                next_state = x + policy(t,x) + W(iw);
                EV = EV + p(iw)*Vs(t+1,next_state);
            end
            Vs(t,x) = EV;
        end
    end
    costs= Vs(1,:);
endfunction

timer()
simulation_dp(optControl)

timer()
