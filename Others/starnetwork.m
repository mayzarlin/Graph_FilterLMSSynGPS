function [H] = starnetwork(node)

% L12 = rand;
% L13 = rand;
% L14 = rand;
% L = [L12+L13+L14 -L12 -L13 -L14;-L12 L12 0 0;-L13 0 L13 0;-L14 0 0 L14];
% [V,D] = eig(L);
% [V L]
SortL = zeros(node.N,node.N,node.N);
        for i = 1:node.N
             CurrL = node.Lk(:,:,i);
             k= 2;
             E = zeros(node.N);
                for j = 1:node.N
                   E(i,1) = 1;
                    SortL(1,:,i) = CurrL(i,:);
                    if( i ~= j && any(CurrL(j,:)) == 1)
                        SortL(k,:,i) = CurrL(j,:);
                        E(j,k) = 1;
                        k=k+1;
                    end
                    
                    
                end
                
                SortL(:,:,i) =  SortL(:,:,i)*E;
                Nsize(i)= length(find(abs(SortL(1,:,i))>0));
%                 SmallL = SortL(1:Nsize,1:Nsize,i);
%                 [V,D] = eig(SmallL);
%                 [~,List] = sort(diag(D));
%                 V_sort = zeros(size(V));
%                    for j =1:length(List)
%                         V_sort(:,j) = V(:,List(j));
%                    end
%                 S(1:Nsize,1:Nsize,i) = V_sort;  %% LocalS has the same eigenvector as original Local Laplacian
                
        end
        
        H = zeros(max(Nsize),2*max(Nsize),node.N);
        for i = 1:node.N
            clear SmallL V D List V_sort L S;
            
            S = zeros(max(Nsize),max(Nsize));
            SmallL = SortL(1:Nsize(i),1:Nsize(i),i);
                [V,D] = eig(SmallL);
                [~,List] = sort(diag(D));
                V_sort = zeros(size(V));
                   for j =1:length(List)
                        V_sort(:,j) = V(:,List(j));
                   end
                S(1:Nsize(i),1:Nsize(i)) = V_sort;
                L =  SortL(1:max(Nsize),1:max(Nsize),i);
                H(:,:,i) = [L S];
        end



end