function plottingTVGraph(Graph,Pos,FigNum)
Ngraph = size(Graph,3);
n      = size(Graph,2);
pg     = 4;
TotalSubplot = round(Ngraph/pg)+1;
figure(FigNum)
nsubp = 1;
clf
hold on
    for j = 1:pg:Ngraph
        subplot(1,TotalSubplot,nsubp) 
        nsubp = nsubp + 1;
        hold on
        for i = 1:n
        neighbor = find(Graph(i,:,j));
              nNgh = length(neighbor);
              for k =1:nNgh
                  P = [Pos(i,:);Pos(neighbor(k),:)];
               line(P(:,1),P(:,2));
              end
        end
        scatter(Pos(:,1),Pos(:,2),50,'filled')%v,node.Y_true,'filled')
    end
end