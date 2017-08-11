% Node class contains properties and functions relating to each node in the Graph.
classdef GF < handle
    
    methods (Static)
function GlobalAdapting(node,Graph)
        Stepsize = Graph.Stepsize/Graph.N;
        node.GlobalX =  (1-Stepsize)*node.GlobalX + Stepsize*Graph.Y; 
end
function GlobalFilter(node,Graph,EV)
I  = eye(Graph.N);
     if(EV == 1) 
        GF.GlobalCovUpdate(Graph);
        S = Graph.EstS;
        D = I;
     else
        S     = Graph.S;
        D     = Graph.D;
     end
     TrueS       = Graph.S; %(end:-1:1,:);
     F_est       = TrueS*node.GlobalX;
     F_true      = S*Graph.Xtrue;
    
     Theta  = zeros(Graph.N);
      for i = 1:Graph.N-1 %% No frequency response at frequency 0.
         if(Graph.D(i,i)>0.00001)
             filtercoef = 2*(F_est(i)^2-F_true(i)^2)/(D(i,i)*F_est(i)^2); %% Find the filter coefficient value
             filtercoef = max(filtercoef,0);
             filtercoef = min(filtercoef,1);
             Theta(i,i) = filtercoef*D(i,i);
         end
     end
     L = S'*Theta*S;
 %% Filtering
     node.GlobalX = (I-L)*node.GlobalX;
end
function GlobalCovUpdate(Graph)
        I =eye(Graph.N);
        Graph.C = Graph.C +(Graph.Y*Graph.Y');
        V = GF.EVtracking(Graph.EstS',Graph.C+I);
        Graph.EstS = V';
end
 %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
function Adapting(node,Graph)
     % Adapting
     X = zeros(Graph.N);
     for i = 1:Graph.N
         e          = zeros(Graph.N,1);
         e(i)       = 1;
         node.X(:,i) =  node.X(:,i) - Graph.Stepsize*node.X(:,i).*e....
                    + Graph.Stepsize*Graph.Y(i).*e;
     end
    
 end
function Combining(node,Graph)
%% Combining
X_inter = zeros(Graph.N);
  for  i = 1:Graph.N
      neighbor = find(Graph.A(i,:));
      n = length(neighbor);
      X = zeros(Graph.N,1);
     
      for j = 1:n
        X = X + (Graph.Adouble(i,neighbor(j)))*node.X(:,neighbor(j));
       
      end
      X_inter(:,i) = X + (Graph.Adouble(i,i))*node.X(:,i);
     
  end
  node.X = X_inter;
end
function LocalFilter(node,Graph,EV) %% This performs on the permutted local Laplacian 
I = eye(Graph.N);
     for i = 1:Graph.N     
         if(EV == 1)
            GF.LocalCovUpdate(Graph,i);
            LocalS 	= Graph.Local(i).EstS;
            D       = I;
         else
             LocalS = Graph.Local(i).S;
             D      = Graph.D;
         end
                 h              = zeros(Graph.N,1);
                 LocalFtrue     = Graph.Local(i).S*Graph.Local(i).E'*Graph.Xtrue; %% Find the energy of each node local true freqeucny (permuted).
                 LocalF         = LocalS*Graph.Local(i).E'*node.X(:,i); %% Find the estimated frequency of node k (permuted)
                 Nfreq          = sum(diag(Graph.Local(i).L)>0);
                 for j = 1:Nfreq
                     h(j) =  2*(LocalF(j)^2-LocalFtrue(j)^2)/(D(j,j)*LocalF(j)^2);
                     h(j) = min(h(j),1);
                     h(j) = max(h(j),0);
                     h(j) = h(j)*D(j,j);
                 end
                node.X(:,i) = node.X(:,i) - Graph.Local(i).E*LocalS'*diag(h)*LocalS*Graph.Local(i).E'*node.X(:,i);
     end
end
function LocalCovUpdate(Graph,i)
             nghList                = abs(Graph.Lsym(i,:))>0; %% find the list of node k neighbor 
             sngh                   = sum(nghList);
             NeighborY              = Graph.Local(i).E'*(nghList' .*Graph.Y);  %% get the nieghbor's measurement data
             Graph.Local(i).C       = Graph.Local(i).C + NeighborY*NeighborY';  %% update permuted local covariance matrix
            %SortL                   = Graph.Local(i).E'*Graph.Local(i).L*Graph.Local(i).E;
             SmallC                 = Graph.Local(i).C(1:sngh,1:sngh) + eye(sngh); %% find the samll version of permuted local covariance matrix
             SmallV                 = Graph.Local(i).EstS(1:sngh,1:sngh)'; %% extract small version of estimated LGFT matrix (Nk times Nk)
             SmallV                 = GF.EVtracking(SmallV,SmallC);     %% Update estimated small LGFT matrix (Nk times Nk)
             Graph.Local(i).EstS(1:sngh,1:sngh) = SmallV';             %% expand  estimated LGFT matrix to size of N
end
function UpdateV = EVtracking(V,L)
n = size(V,2);
Stepsize =1/norm(L);
I = eye(n);
UpdateV = V;
    for i = 1:n
        if(norm(V(:,i))>0)
         W = zeros(n); 
       for j = 1:i-1
            W = W + V(:,j)*V(:,j)';
        end
        UpdateV(:,i) = V(:,i) + Stepsize*(I - V(:,i)*V(:,i)')*(L)*(I - W)*V(:,i);
        UpdateV(:,i) = UpdateV(:,i)/norm(UpdateV(:,i));
        else
            UpdateV(:,i) = V(:,i);
        end
    end
end
end
end



