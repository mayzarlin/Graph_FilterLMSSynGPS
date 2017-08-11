function CovEst_Sim
global  Graph GlobalMSD LocalMSD Algnum EVEstError node


Algnum             = 3;         %% Number of Algorithm.
%%%% Inalize node and graph Class
%%% Time Varying Graph

PlotFig             = 6;
Mvar                = 15^2;     %% Variance of Measurement Noise
Simmax              = 1;        %% Number of Simulation
Mark                = 1;        %Imax/100;
%%%    savefile     = sprintf('%dNet%dNghV%dStepsize0_%dTVGraph',N_number,N_neighbor,Mvar,StepS);
NeworCont           = 1;        %% 1 Start New Simulation 0 continue from previous simulation
NewSimorCont(NeworCont)
Imax                = 100;   
SNR                 = zeros(Imax,1);
splot               = 1;
for a = 1:Algnum
        GlobalMSD(a).msd     = zeros(Imax,1);
        LocalMSD(a).msd      = zeros(Imax,1);
       
end
 EVEstError.local          = zeros(Imax,1);
 EVEstError.global          = zeros(Imax,1);
 EVEstError.plot            = 0;
i = 1; Run = 1;
            while(Run<=Simmax)
 %% Measurement Data
                Graph.UpdateXtrue;
                SNR(i)  = SNR(i) + (1/Simmax)*norm(Graph.Xtrue)^2/(Graph.N*Mvar) ; 
                Noise           =  sqrt(Mvar)*randn(Graph.N,1); 
                Graph.Y         = Graph.Xtrue + Noise;
                for a = 1:Algnum          
                        Algrun(node(a),a);
                        if(a == 2) %% Graph Filter with EigenSpace Estimation 
                        Error = EigenSpaceError;
                        EVEstError.global(i)  =  EVEstError.global(i) ....
                                                + (1/Simmax)*(Error.global);
                        EVEstError.local(i)  =  EVEstError.local(i) ....
                                                 + (1/Simmax)*(Error.local);
                        end
                        %% Update Covariance Matrix and Estimate the EigenVector 
                        GlobalMSD(a).msd(i)      = GlobalMSD(a).msd(i) + (1/Simmax)*node(a).Globalmsd;
                        LocalMSD(a).msd(i)       = LocalMSD(a).msd(i)  + (1/Simmax)*node(a).Localmsd;
                 end
                
                        if(i==Imax)
                            i=1; Run = Run +1;
                            if(Simmax > 1 && Run < Simmax)
                            NewSimorCont(NeworCont);
                            end
                        else
                            i=i+1;
                        end
            end
 AvgSNR = 10*log10(mean(SNR));
 title_text = sprintf('AvgSNR = %f',AvgSNR);
 MSDPlot(PlotFig,Mark,title_text)
 EVEstErrorPlot(EVEstError,PlotFig+1,title_text)
 savingplot(splot,AvgSNR,NeworCont)
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Algrun(node,a)
if(a == 1)
       Distributed_LMS(node);
         
elseif(a==2)
        GFEVtracking_LMS(node);
       
else
        GraphFilter_LMS(node);   
end
end
function ChangeGraph(UsePrevGraph,C)
global Graph TVL TVXtrue Graphcount r
    if(UsePrevGraph == 1)
        load TVLap.mat
        Graph.Lsym = TVL(:,:,Graphcount);
        [EigenV,Eigen]              = eig(Graph.Lsym);
        [Graph.D,EigenV]            = Graph.SortingEV(Eigen,EigenV);
        Graph.S                     = EigenV';
        Graph.Xtrue                 = TVXtrue(:,Graphcount);
    else
        Graph.UpdateL(C)
        Graph.UpdateXtrue(r)
         TVL(:,:,Graphcount) = Graph.Lsym;
         TVXtrue(:,Graphcount) = Graph.Xtrue;
         
    end
Graphcount       = Graphcount+1;
Graph.findLocalL;
Graph.UpdateAdouble;
end
function MSDPlot(PlotFig,Mark,title_text)
global GlobalMSD LocalMSD Algnum 
Plotcolorglobal    = ['r' 'b' 'g'];
GlobalPlotName            = {'LMS' 'GGF-EV' 'GGF'};
LocalPlotName             = {'DGGF' 'LGF-EV' 'LGF'};
if (Mark > 1)
Plotcolorlocal     =  ['+-r'; 'o-g'; 'd-m' ;'s-k']; %% Color for Plots
else
    Plotcolorlocal     =  ['k'; 'm' ;'c' ;'k']; %% Color for Plots
end
l =1;
figure(PlotFig)
            clf
            hold on
            grid on
            for a = 1:Algnum
                if(max(GlobalMSD(a).msd)>0)
                    h = 10*log10(GlobalMSD(a).msd);
                    plot(h,Plotcolorglobal(a),'LineWidth',1.5)
                    LegendList(l) = GlobalPlotName(a);
                    l = l+1;
                end
                if(max(LocalMSD(a).msd)>0)
                    h = 10*log10(LocalMSD(a).msd);
                    plot(1:Mark:length(h),h(1:Mark:end),Plotcolorlocal(a,:),'LineWidth',1.5)
                    LegendList(l) = LocalPlotName(a);
                    l = l+1;
                end
            end
            legend(LegendList);
            title(title_text);
end
function EVEstErrorPlot(EVEstError,PlotFig,title_text)   
            if(EVEstError.plot == 1)
            figure(PlotFig)
            clf
            hold on
            grid on
               hglobal = 10*log10(EVEstError.global);
                hlocal  = 10*log10(EVEstError.local);
                plot(hglobal,'b','LineWidth',1.5);
                plot(hlocal,'m','LineWidth',1.5);
                     
            legend('Global Estimate Error of GFT matrix','Local Estimated Error of LGFT matrix');
            title(title_text);
            end
end
function Graph = initializing(SelectedNode,Nneighbor,SS)
loadfile        = sprintf('30RandNetwork');
load (loadfile); 
Nodepos            = Nodepos(1:SelectedNode,:);
[~,Lsym, A] = adjacency2(Nodepos, Nneighbor(1), Nneighbor(2));
%L = L.*(abs(L)>eps(20^10));
Graph = Network(Lsym,A,SS,SelectedNode);
Graph.UpdateXtrue;
plotting(Graph.Lsym,Nodepos,Graph.Xtrue,1)
end
function Error = EigenSpaceError
global Graph
Ftrue = abs(Graph.S*Graph.Xtrue);
F     = abs(Graph.EstS*Graph.Xtrue);
Error.global = (1/Graph.N)*norm(Ftrue - F)^2; %% Global GFT transform Error
LError = zeros(Graph.N,1);
    for i = 1:Graph.N %% Local LGFT Transform Error
        LocalFtrue  = Graph.Local(i).S*Graph.Local(i).E'*Graph.Xtrue;
        LocalF      = Graph.Local(i).EstS*Graph.Local(i).E'*Graph.Xtrue;
        LError(i)   = norm(abs(LocalFtrue) - abs(LocalF))^2;
    end
Error.local  = mean(LError);
end
function savingplot(splot,AvgSNR,NeworCont)
    if(splot ==1)
     h1 = figure(PlotFig);
     h2 = figure(PlotFig+1);
     index = 1; contindex = 1;
     while(1)
         s1 = sprintf('figure/30MSDALLSNR%d_%d.fig',round(AvgSNR),index);
         s2 = sprintf('figure/30EvErrorSNR%d_%d.fig',round(AvgSNR),index);
         if( exist(s1,2))
             if(NeworCont == 0) % if the simulation is cont
                 while(1)
                s1 = sprintf('figure/30MSDALLSNR%d_%d_%d.fig',round(AvgSNR),index,contindex);
                s2 = sprintf('figure/30EvErrorSNR%d_%d_%d.fig',round(AvgSNR),index,contindex);
                    if(exist(s1,2))
                        contindex = contindex+1;
                    else
                        break;
                    end
                 end
             end
             index = index + 1;
         else
         savefig(h1,s1)
         savefig(h2,s2)
         break;
         end
     end
    end
end
function NewSimorCont(NeworCont)
global Graph node Algnum
    if (NeworCont== 1) %% Staring the New Simuation
        SelectedNode       = 10;
        Nneighbor          = [4 3];     %% maximum number of neighbors
        SS                 = 0.1;       %% Stepsize
        Graph = initializing(SelectedNode,Nneighbor,SS);
        for a = 1:Algnum 
            node(a)  = Node(SelectedNode);
        end
    else %% Continume Prev Simulation 
        load('ContSim.mat')
    end
end








% if(TrackSignal.plot == 1)
% TrackSignal.LMS             = zeros(SelectedNode,Imax);
% TrackSignal.GGF             = zeros(SelectedNode,Imax);
% TrackSignal.TrueSig         = zeros(SelectedNode,Imax);
% TrackSignal.GlobalEst       = zeros(SelectedNode,Imax);
% end

% function TrackingSignalPlot(PlotFig)
% global TrackSignal
% if(TrackSignal.plot ==1)
%             figure(PlotFig)
%             clf
%             hold on 
%             grid on
%             plot(TrackSignal.LMS','k')
%             plot(TrackSignal.GGF','b')
%             plot(TrackSignal.TrueSig','r')
%             plot(TrackSignal.LGF','g')
% end
% end
 %TrackingSignalPlot(PlotFig+2)
% if(a == 1 && Graph.LMS == 1 && TrackSignal.plot == 1)
%     TrackSignal(1).LMSEst(:,i) = node(a).GlobalX;
% end
% if(a == 3 && Graph.GGF == 1 && TrackSignal.plot == 1)
%     TrackSignal(1).GlobalEst(:,i) = node(a).GlobalX;
%     TrackSignal(1).TrueSig(:,i)   = Graph.Xtrue;
% end
% if(a == 4 && Graph.GGFEV == 1 && TrackSignal.plot == 1)
%     TrackSignal(2).GlobalEst(:,i)   = node(a).GlobalX;
% end