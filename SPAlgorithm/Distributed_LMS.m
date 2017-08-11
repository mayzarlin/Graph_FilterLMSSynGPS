function node = Distributed_LMS(node)
global Graph
%% Global LMS
if(Graph.LMS == 1)
    GF.GlobalAdapting(node,Graph);
    node.Globalmsd = (1/Graph.N)*norm(Graph.Xtrue - node.GlobalX)^2;
end

if(Graph.DGGF == 1)

 for i = 1:node.N
         node.F_inter(:,i) =node.F(:,i) -node.Stepsize*node.X(i,i)*node.S(:,i)....
                         + node.Stepsize* node.Y(i).*node.S(:,i);
 end
%%%%% Combining Frequency
F_comb = zeros(node.N);
 for  i = 1:node.N
      neighbor = find(node.A(i,:));
      n = length(neighbor);
      F = zeros(node.N,1);   
      for j = 1:n
             F = F + (node.Adouble(i,neighbor(j)))*node.F_inter(:,neighbor(j)); 
      end
      F_comb(:,i) = F + (node.Adouble(i,i))*node.F_inter(:,i);
 end
%%%%%%%%%%%%%% node K filter Design for Global Filter

 Pre_F = node.S*node.Xtrue;
 D = eye(node.N);
 
for j = 1:node.N
        h   = zeros(node.N); 
        
     if (norm(F_comb(:,j)) > norm(Pre_F))
        for i =2:node.N
                    h(i,i) =2*(F_comb(i,j)^2-Pre_F(i)^2)/(D(j,j)*F_comb(i,j)^2);
                    if h(i,i) <= 0 
                        h(i,i) = 0;
                    end
                    if h(i,i) >2
                        h(i,i) =2;
                    end
        end
     end
     h(1) = node.D(1,1);
     LGF =h*D;
     node.F(:,j) = (eye(node.N)-LGF)*F_comb(:,j);
end
%%% Update %X% 
X_inter = node.X;
for i = 1:node.N
       X_inter(i,i) = node.S(:,i)'*node.F(:,i);
end
RMS = zeros(node.N,1);
 for  i = 1:node.N
      neighbor = find(node.A(i,:));
      n = length(neighbor);
      X = zeros(node.N,1);
    
      for j = 1:n
       X = X +(node.Adouble(i,neighbor(j)))* X_inter(:,neighbor(j));
       
      end
    node.X(:,i) = X +(node.Adouble(i,i))* X_inter(:,i);
    RMS(i) = (1/node.N)*norm(node.Xtrue - node.X(:,i))^2;
end 
node.Localmsd = mean(RMS);

end
end




% % %% Distributed LMS
%   for i = 1:node.N
%       e = zeros(node.N,1);
%       e(i) = 1;
%     node.X_inter(:,i) = node.X(:,i) - node.Stepsize*node.X(:,i).*e....
%                         + node.Stepsize*node.Y(i).*e;
%      e(i) = 0;
%   end
%   for  i = 1:node.N
%       neighbor = find(node.A(i,:));
%       n = length(neighbor);
%       X = zeros(node.N,1);
%       for j = 1:n
%         X = X + (node.Adouble(i,neighbor(j)))*node.X_inter(:,neighbor(j));
%       end
%       node.X(:,i) = X + (node.Adouble(i,i))*node.X_inter(:,i);
%      RMS(i) = (1/node.N)*norm(node.Y_Lap - node.X(:,i))^2;
%   end
%   node.Mean_Rms = mean(RMS);
% end