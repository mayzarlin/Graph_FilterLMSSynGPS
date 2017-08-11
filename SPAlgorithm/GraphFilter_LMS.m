function [node] = GraphFilter_LMS(node)
global Graph
if(Graph.GGF == 1) 
        
        %% Global Filtering
        GF.GlobalAdapting(node,Graph);
        GF.GlobalFilter(node,Graph,0);
        node.Globalmsd = (1/Graph.N)*norm(Graph.Xtrue - node.GlobalX)^2;
end
RMS = zeros(Graph.N,1);
if(Graph.LGF == 1)
%% LOCAL ALG ACG 
    GF.Adapting(node,Graph);
    GF.Combining(node,Graph);
    GF.LocalFilter(node,Graph,0);
    for i = 1:Graph.N
     RMS(i) = (1/Graph.N)*norm(Graph.Xtrue - node.X(:,i))^2;
    end
elseif(Graph.LGF == 2)
  
    GF.Adapting(node,Graph);
    GF.LocalFilter(node,Graph,0);
    GF.Combining(node,Graph);
    for i = 1:Graph.N
     RMS(i) = (1/Graph.N)*norm(Graph.Xtrue - node.X(:,i))^2;
    end
end
    node.Localmsd = mean(RMS);  
end





