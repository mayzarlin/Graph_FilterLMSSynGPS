function node = DRegularized_LMS(node)
global  Graph
%% Regularized Global LMS
%% Adapting
if(Graph.GRLMS == 1)
     node = GF.GlobalAdapting(node,Graph);
     node = RLMS.GlobalReg(node,Graph);
    node.Globalmsd = (1/Graph.N)*norm(Graph.Xtrue-node.GlobalX)^2;
end

RMS = zeros(Graph.N,1);
 %%%%%%%%%%%%%% New Distributed 
        if(Graph.DRLMS == 1) %%%%ACG
            node = GF.Adapting(node,Graph);
            node = GF.Combining(node,Graph);
            node = RLMS.LocalReg(node,Graph);
            for i= 1:Graph.N
                RMS(i) = (1/Graph.N) * norm(Graph.Xtrue -node.X(:,i))^2;
            end
        elseif(Graph.DRLMS == 2) %%% AGC
            node = GF.Adapting(node,Graph);
            node = RLMS.LocalReg(node,Graph);
            node = GF.Combining(node,Graph);
            for i = 1:node.N
                RMS(i) = (1/node.N)*norm(Graph.Xtrue - node.X(:,i))^2;
            end
        end
    node.Localmsd = mean(RMS);
end







