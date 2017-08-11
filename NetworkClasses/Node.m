% Node class contains properties and functions relating to each node in the network.
classdef Node < handle
    properties
        %% Network properties
        X              % Local Estimate of X 
        GlobalX     % Global Estimate of X
        F
        Localmsd  = 0;
        Globalmsd = 0;
    end
    methods 
    function node = Node(n)
    node.X        = zeros(n);
    node.GlobalX  = zeros(n,1);
    node.F        = zeros(n,1);
     end
    end
end


