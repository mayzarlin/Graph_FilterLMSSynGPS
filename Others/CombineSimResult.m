function CombineSimResult
close all
[filename,pathname] = uigetfile('figure/figure','MultiSelect','on');
 n = length(filename);
    for i = 1:n
        fig_toopen = char(fullfile(pathname,filename(i)));
        h = open(fig_toopen);
        H = findobj(h,'type','line');
        X = get(H,'ydata');
        X = cell2mat(X);
        Data(i).X = X';
        close all
    end
    CData = [Data(1).X;Data(2).X(100:end,:)];
    plot(CData)
    newfigname  = strcat('Combined',filename(1));
    newfilename = fullfile(pathname,newfigname);
    savefig(newfilename)
    
end