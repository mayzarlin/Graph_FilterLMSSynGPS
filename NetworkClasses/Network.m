% Node class contains properties and functions relating to each node in the network.
classdef Network < handle
    properties
        %% Network properties
        LMS       = 1;    %%% 1 runs the LMS algorithm
        GGF       = 1;    %%% 1 runs GGF alg.
        DGGF      = 0;    %%% 1 runs LGGF alg.
        GGFEV     = 1;    %%% 1 runs GGF alg with EV tracking
        LGF       = 1;    %%% 1 LGF ACG, 2 LGF AGC 0 Do not run LGF
        LGFEV     = 1;    %%% 1 LGF ACG with EV tracking
                          %%% 2 LGF AGC with EV tracking and 0 Do not run
                          %%% LGF EV
       EVtrackingIteration = 1;
       TVGraph   = 0;
       A            % Adjacency Matrix
       Lsym        % Symmetric Laplacian Matrix
       Xtrue       % True Solution for X
       Y            % Observed Data Value
       N            % Number of Node in the network       
       D            % Eigen Value of Laplacian (Sorted from Low to High)
       S            % Transpose of Eigen Vector Matrix of Laplacian (Sorted from Low to High)
       EstS         % Estimated Eigen Vector Matrix;a
       Adouble     % double Stochastic Matrix
       Local
       Stepsize     % Node.Stepsize
       theta        % TV recursion
       C            % Covariance Matrix
    end
    
    methods 
        
function Graph = Network(Lap,Adj,SS,Nm)    
            Graph.A                  = Adj; %% Adjacency matrix
            Graph.Xtrue              = zeros(Nm,1); %% True Graph Signal 
            Graph.N                  = Nm;   %% Number of node in the graph
            Graph.Lsym               = Lap;  %% Normalized Laplacian
            Graph.Y                  = zeros(Graph.N,1); %% Graph Noisy Measurement
            [EV,Eigen]               = eig(Lap);  %% Find the eigendecomposition of Normalized Laplacian
            [Graph.D,V]              = Graph.SortingEV(Eigen,EV); %% Sorting the EigenValue and Eigen Vector Matrices
            Graph.S                  = V';   %% Sorted Graph Foruier Transform Matrix 
            Graph.EstS               = rand(Graph.N); %% Estimated Graph Foruier Transform Matrix
            Graph.Stepsize           = SS;     %% Learning Stepsize
            Graph.theta              = 0;       %% Recurssive parameter for time varying signal model 
            for i = 1:Graph.N
                Lk(i)           = LocalNode(Nm); %% Contain class of LocalNode X,L,D,S,EstS and C for each node
            end
            Graph.Local         = Lk;   
            Graph.C             = zeros(Graph.N); %% Global Covariance Matrix
            Graph.findLocalL                       %% Find each node k Local laplacian and its eigen decomposition
            for i = 1:Graph.N
                Graph.Local(i).EstS = rand(Graph.N).*(abs(Graph.Local(i).S)>0); %% Estimated Local GFT matrix initialize with some random value at each column of true Local GFT  
            end
            Graph.UpdateAdouble;        %% Weighted Adjacency Matrix for Combining steps in Distributed Algorithm
            
     end
function UpdateL(Graph,C)
Link = [];
Ltemp = Graph.Lsym; %% To check the difference between old and new L
while(C>size(Link,1))  %% C is the max number of update link.
    ED   = randn;      %% ED decides whether to add an edge or remove an edge
    if(ED >0) %%% Adding the edge or Weight changes on the existing edge
            Link = Graph.addingtheEdge(Link);
    else %% Removing the edge
            Link = Graph.removingtheEdge(Link);
    end
end
            [EigenV,Eigen]              = eig(Graph.Lsym); %%the Eigen Decomposition of New Laplacian
            [Graph.D,EigenV]            = Graph.SortingEV(Eigen,EigenV); %% Sorting the new Eigen Decomposition 
            Graph.S                     = EigenV';
end
function findLocalL(Graph)  
        for  i = 1:Graph.N %% Extracking node k local L from Global Laplacian
            Graph.Local(i).L(:,i)       =(0.5)* Graph.Lsym(:,i); 
            Graph.Local(i).L(i,:)       = (0.5)* Graph.Lsym(i,:);
            for j = 1:Graph.N
                Graph.Local(i).L(j,j)   = (-0.5)* Graph.Lsym(i,j);
            end
            Graph.Local(i).L(i,i)       = 0;
            Graph.Local(i).L(i,i)        = 0.5+0.5*sum( Graph.Lsym(:,i)); % Update the diagonal value for normalzied L
        end
Graph.LocalPermuatationMat;     %% Find the Permuatation matrix for each local L
Graph.LocalEigenDecomposition;  %% Find the eigendecomposition for each local L
end
function UpdateAdouble(Graph)
NghMat          = Graph.Lsym ~= 0;  %% Logical Matrix with 1 at the neighboring element
Graph.Adouble   = 1./sum(NghMat,2).*NghMat; %% find the Weighted Adjacency
end
function UpdateXtrue(Graph)
    EigenV = Graph.D;
    EigenV(1,1) = 0;
    q = sqrt(0.01)*randn(Graph.N,1); %% random walk vector
    Graph.theta = 0.9*Graph.theta + Graph.S'*pinv(sqrt(EigenV))*q; %% Update the TV recurssive 
    Graph.Xtrue = 5*ones(Graph.N,1) + Graph.theta;  %% Update the Graph Signal
end
function LocalPermuatationMat(Graph)
    for i = 1:Graph.N %% Find the permuataion matrix for each node k
           CurrL    = Graph.Local(i).L;
           k        = 2;
           tempE    = zeros(Graph.N);
                  for j = 1:Graph.N
                       tempE(i,1) = 1;
                        if( j~=i && CurrL(j,j)>0)
                            tempE(j,k) = 1;
                            k=k+1;
                        end
                  end
          Graph.Local(i).E = tempE;
    end
end
function LocalEigenDecomposition(Graph)
    for i = 1:Graph.N
        ngh         = sum(Graph.Local(i).L(i,:)~=0); %% Number of node k neighbor 
        BigL        = Graph.Local(i).E'*Graph.Local(i).L*Graph.Local(i).E; %% permute node k laplaian to star network
        SmallL      = BigL(1:ngh,1:ngh); %% get the small Local laplacian only with node k's neighbor
        [SV,SD]     = eig(SmallL);       %% eigen decomposition for small local laplacian
        [SD,SV]     = Graph.SortingEV(SD,SV); %% Sorting the eigen vlaue
        Graph.Local(i).D(1:ngh,1:ngh)   = SD;  %% Permuted EigenValue  Matrix of local L 
        Graph.Local(i).S(1:ngh,1:ngh)   = SV'; %% Permuted EigenVector Matrix of local L 
    end
    
end
function [Link] = addingtheEdge(Graph,Link)
    loop = 1; Ltemp = Graph.Lsym;
    while(loop<Graph.N) % Link Selection Loop
          
                    % l = randi(n);
                    l = randi(Graph.N);
                    k = randi(Graph.N);
                   
                    %% Update the weight if the select link is an existing link
                    if( l~= k && findLK(Link,l,k)) % && PrevL(l,k)<0)
                        if (isempty(Link))
                            Link = [l k 1];
                        else
                            Link = [Link;l k 1];
                        end
                        wtemp =(-1)* exp(-1*(randn-randn)^2); %% New weight 
                        Ltemp(l,k) = wtemp;
                        Ltemp(k,l) = wtemp;
                        Ltemp(l,l) = 0;
                        Ltemp(l,l) = abs(sum(Ltemp(l,:))); %% update the diagonal element in Laplacian 
                        Ltemp(k,k) = 0;
                        Ltemp(k,k) = abs(sum(Ltemp(k,:)));
                        E = sort(eig(Ltemp));
                        if(abs(E(2)) > 0.0001 || abs(Ltemp(l,k))>0.01) % check the new graph is connected.
                           Graph.Lsym = Ltemp; %%
                           loop = Graph.N;
                        end
                     end
                    loop = loop +1;
    end
      
end
function [Link] = removingtheEdge(Graph,Link)
    Ltemp = Graph.Lsym; loop = 1;
    while(loop<Graph.N)
     l = randi(Graph.N);
     k = randi(Graph.N);
        if(l~=k && findLK(Link,l,k) && Ltemp(l,k)<0)
            if(isempty(Link))
                Link = [l k -1];
            else
                Link = [Link;l k -1];
            end
            Ltemp(l,k) = 0;
            Ltemp(k,l) = 0;
            Ltemp(k,k)  = 0;
            Ltemp(l,l)  = 0;
            Ltemp(k,k)  = abs(sum(Ltemp(k,:)));
            Ltemp(l,l)  = abs(sum(Ltemp(l,:)));
            E = sort(eig(Ltemp));
            if(abs(E(2))>0.01)
                Graph.Lsym = Ltemp; 
                loop = Graph.N;
            end
        end
        loop = loop +1;
    end
    
end
    end

methods(Static)
    function [Es,Vs] = SortingEV(E,V)
[es,List] = sort(diag(E),'ascend');
Es        = diag(es);
n = length(List);
Vs        = zeros(n);

        for i = 1:n
            Vs(:,i) = V(:,List(i));       
        end

    end
end
end
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
 function t =  findLK(Link,l,k)
     n = size(Link,1);
     t = true;
     for i = 1:n
        if(k == Link(i,1))
            if(l == Link(i,2))
                t = false;
            end
        end
        if(l == Link(i,1))
            if(k == Link(i,2))
                t = false;
            end
        end
     end
 end
%%%%%%%%%%%%
%%%%%%%%%%%%


% 
%  [U,E]               = eig(Graph.Local(i).L);
%         [ Graph.Local(i).D,V] = Graph.SortingEV(E,U);
%         V = V.*(abs(V)>eps(20^10));
%         
%  %%%% Remove 1 eigenvector 
%       %  Nnumber = find(Graph.Local(i).L(i,:));
%         for k = 1:Graph.N
% %             nonzeroEV = find(abs(V(:,k))>0);
% %             if(length(nonzeroEV) < length(Nnumber))
% %                 V(:,k) = zeros(Graph.N,1);
% %             end
%         s = (V(:,k)~=0);
%         ngh = (Graph.Local(i).L(i,:) ~= 0 );
%         if(abs(Graph.Local(i).D(k,k))<0.00001 && sum(s) < sum(ngh))
%             V(:,k) = zeros(Graph.N,1);
%         end
%         end
%         Graph.Local(i).S  = V'; 


% function UpdateL01(Graph,NnodeChange)
%     Ltemp     = Graph.Lsym;
%     PrevX = Graph.Xtrue; EdgeList = zeros(NnodeChange);
%    for i = 1:NnodeChange
%        while(1)
%         k = randi(Graph.N);
%            if(isempty(find(EdgeList==k,1)))
%                 PrevX(k) =  Graph.Xtrue(k) + randn;    
%                 Nneighbor = find(Graph.Lsym(k,:));
%                 Ltemp(k,k)= 0;
%                 for l = 1:length(Nneighbor)
%                     if(k ~= Nneighbor(l))
%                     Ltemp(k,Nneighbor(l)) = (-1)*EdgeWeight(PrevX(k),PrevX(Nneighbor(l)));
%                     Ltemp(Nneighbor(l),k) =  Ltemp(k,Nneighbor(l));
%                     Ltemp(Nneighbor(l),Nneighbor(l)) = 0;
%                     Ltemp(Nneighbor(l),Nneighbor(l)) = abs(sum(Ltemp(Nneighbor(l),:)));
%                     end
%                 end
%                 Ltemp(k,k)   = abs(sum(Ltemp(k,:)));
%                 EdgeList(i) = k;
%                 break;
%            end
%        end
%    end
%    Graph.Xtrue = PrevX;
%    Graph.Lsym = Ltemp;
% %%% Update Weight in Laplacian
%    
% end
