function [L,Lsym, A] = adjacency2(locations, neighMin, neighMax)
% locations = matrix of [x-coordinate, y-coordinate]
% neighMin = minimum number of neighbors permitted for each node
% neighMax = maximum number of neighbors permitted for each node

numStations = size(locations, 1);         % number of stations
threshold = 0.01;                          % threshold of initial pruning
theta = -5;                            % constant for the Gaussian distribution

C_Weighted = zeros(numStations);          % weighted adjacency matrix of complete graph
C = ones(numStations) - eye(numStations); % unweighted adjacency matrix of complete graph

% Form the weighted adjacency matrix of the complete graph
for i = 1:numStations
    for j = i+1:numStations
        % C_Weighted(i,j) = exp(theta*(norm(abs(locations(i,2)) - abs(locations(j,2))))^2);
        C_Weighted(i,j) = EdgeWeight(locations(i,:),locations(j,:));
        C_Weighted(j,i) = C_Weighted(i,j);
    end
end

% Use 'removeAndInsertEdges' to give a pruned version of 'C_Weighted' and 'C'. This uses 
% 'threshold' to delete edges, then 'neighMin' to re-insert certain edges in case it leads 
% to a node having too few neighbors.
[A_Weighted, A] = removeAndInsertEdges(C_Weighted, C, threshold, neighMin);

% The pruned graph may not be connected, so lower the threshold for connections until it 
% is connected.
while ~isConnected(A_Weighted)
    threshold = threshold - 0.05;
    [A_Weighted, A] = removeAndInsertEdges(C_Weighted, C, threshold, neighMin);
end

% Use 'reRemoveEdges' to do some more pruning. This uses 'neighMax' to remove low-weight 
% edges surrounding highly connected nodes, while making sure the removals neither lead to 
% a disconnected graph nor violate 'neighMin'.
[A_Weighted, A] = reRemoveEdges(A_Weighted, A, neighMin, neighMax);

% Make Laplacian matrix of the graph
D = diag(sum(A_Weighted, 2));
L = D - A_Weighted;
%Lsym = L;
Lsym = sqrt(D)\L/(sqrt(D));

end


function [A_Weighted, A] = removeAndInsertEdges(C_Weighted, C, threshold, neighMin)
% First remove edges with values below 'threshold'. Then add edges to certain nodes if 
% the removal leaves those nodes with below 'neighMin' neighbors.

% Prune away edges with weights below 'threshold'
[A_Weighted, A] = pruneEdgesUnderThreshold(C_Weighted, C, threshold);
neighNumA = findNeighNumA(A);    % number of neighbors of each node
lows = find(neighNumA < neighMin);    % find nodes with too few neighbors

% Loop through nodes with fewer than neighMin neighbors. Add highest valued possible 
% edges for each one (even if they are below 'threshold').
while ~isempty(lows)
    n = lows(1);                                  % get the first lows node, 'n'
    neighNum = size( find( A(n,:) ), 2 );         % n's number of neighbors
    
    % Add edges to the node n until it has at least neighMin neighbors.
    while neighNum < neighMin
        maxWeight = max( C_Weighted(n,:) );      % maximum possible edge weight for n
        highestNeigh = find( C_Weighted(n,:) == maxWeight );   % the neighbor from that edge
        A_Weighted(n,highestNeigh) = maxWeight;      % changing the adjacency matrices
        A_Weighted(highestNeigh,n) = maxWeight;
        A(n,highestNeigh) = 1;        A(highestNeigh,n) = 1;
        
        neighNum = neighNum + 1;
    end
    
    lows(1) = [];
end

end


function [A_Weighted, A] = pruneEdgesUnderThreshold(C_Weighted, C, threshold)
% Delete edges between nodes that are below 'threshold' in (complete) weighted adjacency 
% matrix 'C_Weighted' and unweighted adjacency matrix 'C'

numStations = size(C,1);
A_Weighted = repmat(C_Weighted, 1);
A = repmat(C, 1);

% Delete entries in D_Weighted and D that are below threshold
for i = 1:numStations
    for j = i+1:numStations
        if C_Weighted(i,j) < threshold
            A_Weighted(i,j) = 0;     A(i,j) = 0;
            A_Weighted(j,i) = 0;     A(j,i) = 0;
        end
    end
end

end


function [A_Weighted, A] = reRemoveEdges(C_Weighted, C, neighMin, neighMax)
% 'C_Weighted' and 'C' have too many edges. This subroutine returns 'A_Weighted' and 'A' 
% as pruned versions of C_Weighted and C, in accordance with 'neighMin' and 'neighMax', as 
% well as with being connected.

A_Weighted = C_Weighted; A = C;
numStations = size(A,1);

% Loop through the stations to see how to prune away edges up to neighMax while staying in
% accordance with neighMin and connectedness.
for i = 1:numStations
    neigh = find(A(i,:));              % neighbors of i
    neighNum = size(neigh,2);          % number of neighbors of i (subject to change)
    
    neighWeights = sort(A_Weighted(i,neigh));   % original neighbor weights in sorted order
    neighNum0 = neighNum;              % original number of neighbors
    counter = 1;         % keeps track of how many original neighbors have been seen
    
    % Make temporary copies of A_Weighted and A to check if an edge removal leads to a 
    % disconnected graph or too few neighbors to the neighbor
    B_Weighted = A_Weighted; B = A;
    
    % While i's number of neighbors is above neighMax and there are neighbors that could 
    % be removed, experiment with removing edges in order of increasing weights
    while neighNum > neighMax && counter < neighNum0
        % Find the edge with minimum weight, take it away from B_Weighted and B.
        smallNeigh = find(B_Weighted(i,:) == neighWeights(counter));
        B_Weighted(i,smallNeigh) = 0; B(i,smallNeigh) = 0;
        B_Weighted(smallNeigh,i) = 0; B(smallNeigh,i) = 0;
        
        % If the resulting graph is disconnected or has a station with too few neighbors, 
        % retrace to the original graph and move on to the next smallest edge.
        if ~isConnected(B_Weighted) || size(find(B(smallNeigh,:)), 2) < neighMin
            B_Weighted = A_Weighted; B = A;
            
        % Otherwise, copy B_Weighted into A_Weighted.
        else
            A_Weighted = B_Weighted; A = B;
            neighNum = neighNum-1;
        end
        
        counter = counter+1;    % increment how many original neighbors we have seen
    end
end

end


function neighNumA = findNeighNumA(A)
% Find number of neighbors of each node based on adjacency matrix 'A' (could be weighted 
% or unweighted)

numStations = size(A,1);
neighNumA = zeros(1,numStations);

for i = 1:numStations
    neighNumA(i) = size(find(A(i,:)), 2);
end

end


function bool = isConnected(A)
% Check if the adjacency matrix 'A' corresponds to a connected graph

% Make Laplacian matrix of the graph
D = diag(sum(A, 2));
Lsym = D - A;
%Lsym = sqrt(D)\L/(sqrt(D));

% If there are more than 1 eigenvalues equal to 0, then the graph is not connected
E = eig(Lsym);
[numZero, ~] = size(find(E<0.0001));
bool = ~(numZero > 1);

end