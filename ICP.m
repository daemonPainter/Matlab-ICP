function [ T , R ] = ICP( S , M )
%ICP Iterative Closest Point algorithm
%   [T,R] = ICP(S,M) applies the Iterative Closest Point algorithm to the
%   matrix S (source) trying to re-arrange its components so they match the best
%   order wrt the order of the points in the M (model) matrix. S,M must be
%   4-by-n matrices with points in homogeneous coordinates (i.e., 4th row
%   must be a 1s row).
%   T is the homogeneous transform that registers R to M. R is a 4-by-n
%   matrix with the same convention as above that contains points taken
%   from S that match the best order wrt to M.


if size(S,1)~=4 || size(M,1)~=4
    fprintf('Input matrices must be 4-by-N');
    return
end
    for i=1:size(S,2)
        if S(4,i)~=1
            fprintf('Input matrices must follow this convention:\n');
            fprintf('dimension:4-by-n; rows 1 to 3 are cartesian coordinates, row 4 contains only ones (homogeneous coordinates)\n');
            return
        end
    end
    for i=1:size(M,2)
        if M(4,i)~=1
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

while( delta_fcost > threshold && iter < I_max )
    iter = iter + 1;
    % transform the source points with the guessed T_k transform
    if iter>1
        source_t = T_k(:,:,iter-1) * S;
    else
        source_t = T_k * S;
    end
    % Compute corrispondence between the transformed source points (source)
    % and the transformed points (model)
    [source_tc] = closest_points(source_t,M);
    % Inputs as 4 x n
    source_tc = [source_tc(1:3,:); ones(1,size(source_tc,2))];
    % Compute Corrispondent point Registration
    T_new = quaternion_matching(source_tc,M);
    % Update the guessed Transform Matrix
    if iter>1
        T_new = T_new * T_k(:,:,iter-1);
    else
        T_new = T_new * T_k;
    end
    % Compute updated source transformed points
    source_tct = T_new * S;
    T = T_new;  % OUTPUT
    % Compute corrispondence and find the total error
    [model_c] = closest_points(source_tct,M);
    e = rms(model_c(1:3,:)-source_tct(1:3,:));

    % The cost function is the sum of the RMS errors
    fcost(iter) = sum(e);
    % Save the matrix T_iter
    T_k(:,:,iter) = T;
    
    % saves output
    R = [model_c; ones(1,size(model_c,2))];
    
    % Compute the delta between the current and the previous
    % iteration - if the delta is lower than threshold
    % the algorithm stops
    if(iter <= 2)
        delta_fcost = 1000;
    else
        delta_fcost = abs(fcost(iter) - fcost(iter-1));
    end
end


end

