function [node] = GFEVtracking_LMS(node)
global  Graph
if(Graph.GGFEV == 1)
        GF.GlobalAdapting(node,Graph);
        %% Global Filtering
        GF.GlobalFilter(node,Graph,1);
        node.Globalmsd = (1/Graph.N)*norm(Graph.Xtrue - node.GlobalX)^2;
end
RMS = zeros(Graph.N,1);
if(Graph.LGFEV == 1)
%% LOCAL ALG ACG 
    GF.Adapting(node,Graph);
    GF.Combining(node,Graph);
    GF.LocalFilter(node,Graph,1);
    for i = 1:Graph.N
     RMS(i) = (1/Graph.N)*norm(Graph.Xtrue - node.X(:,i))^2;
    end
elseif(Graph.LGFEV == 2)
  
    GF.Adapting(node,Graph);
    GF.LocalFilter(node,Graph,1);
    GF.Combining(node,Graph);
    for i = 1:Graph.N
     RMS(i) = (1/Graph.N)*norm(Graph.Xtrue - node.X(:,i))^2;
    end
end
    node.Localmsd = mean(RMS);  
end














