function plotting(Lsym,Pos,Y_true,figNum)
n = size(Lsym,2);
figure(figNum);
clf
hold on;
        for i = 1:n
              neighbor = find(Lsym(i,:));
              n = length(neighbor);
              for j =1:n
                  P = [Pos(i,:);Pos(neighbor(j),:)];
               line(P(:,1),P(:,2));
              end
        end
 scatter(Pos(:,1),Pos(:,2),50,Y_true,'filled')%v,node.Y_true,'filled');
 T = num2str(Y_true);
 text(Pos(:,1),Pos(:,2)+0.05,T)
end