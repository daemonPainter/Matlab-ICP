function [ T , Y ] = ICP( P , X )
%ICP Iterative Closest Point algorithm
%   [T,Y] = ICP(P,X) applies the Iterative Closest Point algorithm to the
%   matrix P so that it is registered to the model X. Since points in P may
%   not have the same order they have in model X, they are re-arranged so
%   that they match their best closest neighbour. Re-ordered and registered
%   P set is saved in Y, while T contains the homogeneous transform that
%   can register P to X reference frame (i.e. X = T*P)
%   P,X (and Y) must be 4-by-n matrices with points in homogeneous
%   coordinates (i.e. [x;y;z;1]).

if size(P,1)~=4 || size(X,1)~=4
    fprintf('Input matrices must be 4-by-N\n');
    return
end
    for i=1:size(P,2)
        if P(4,i)~=1
            fprintf('Input matrices must follow this convention:\n');
            fprintf('dimension:4-by-n; rows 1 to 3 are cartesian coordinates, row 4 contains only ones (homogeneous coordinates)\n');
            return
        end
    end
    for i=1:size(X,2)
        if X(4,i)~=1
            fprintf('Input matrices must follow this convention:\n');
            fprintf('dimension:4-by-n; rows 1 to 3 are cartesian coordinates, row 4 contains only ones (homogeneous coordinates)\n');
            return
        end
    end    
        
I_max = 1000;       % max num of iteration
threshold = 1e-99;  % accuracy threshold
T0 = eye(4);        % initial guess
T_k = T0;           % guess at iter k
iter = 0;           % # of the current iteration
fcost = 1000;       % function cost to be minimized
delta_fcost = 1000; % delta fcost between two iterations

k = 0;
P_0 = P;
T = eye(4);

while( delta_fcost > threshold && k < I_max )
    if k>0  % can we remove this if?
        Y_k = closest_points(P_k,X);
    else
        Y_k = closest_points(P_0,X);
    end
    T_k = quaternion_matching(Y_k,X);
        T = T_k*T;  % update wrt previous step
    P_k = T*P_0;
    
    Y = closest_points(P_k,X);  % need to re-arrange points again
    Y = [Y;ones(1,size(Y,2))];
    % the cost function is the sum of the RMS errors
    e = rms(Y(1:3,:)-P_k(1:3,:));
    fcost(k+1) = sum(e);
    % compute error riduction between two successive steps
    if k==0
        if T == eye(4)
            fprintf('P and X are the same set\n');
            return
        end
        delta_fcost = 1000;
    else
        delta_fcost = abs(fcost(k+1) - fcost(k));
    end
    k = k+1;
end

% Transform
display(['Stopped at iteration: ' num2str(k)])
display(['Cost function: ' num2str(fcost(end))])

end

