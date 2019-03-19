%% Clear workspace
    close all; clc; format compact;

    
%% NOTE: Solution step figures are for 1st case (can easily be generalized)
%
%

    %%
    % 
    % <<Fig1.PNG>>
    % 
    %%
    % 
    % <<Fig2.PNG>>
    % 
%% Given constants 
     element.Area = 150.0;
     prob = 1; % 1 = HW2.P1, 2 = Ch. 2 Ex, 3 = Ch. 1 P. 8
     
 %% Define No. of elements/nodes    
     switch prob
         case {1} % HW2 Problem 1 
            element.Ufactor = [5.88, 2.27, 10.0, 0.581, 0.781, 2.22, 1.47];   
            element.N = 7; % No. of elements (Make sure element.N == length(element.Ufactor))
         case {2} % Ch 1. Ex.    
           element.Ufactor = [5.88, 1.23, 0.76, 0.091, 2.22, 1.47];  
           element.N = 6;
         case {3}% Values for Ch 1. P. 8
           element.Ufactor = [5.88, 1.23, 0.76, 0.053, 2.22, 1.47]; 
           element.N = 6;
     end
     node.N = element.N + 1; % No. of nodes
     element.k = element.Ufactor.*element.Area; % K-Matrix coefficients (coupling coeff.'s)
	node.loads = zeros(node.N,1);
	node.loads(1) = 10;
	node.loads(end) = 68;
     % 
     % <<Fig3.PNG>>
     % 
     %%
     % 
     % <<Fig4.PNG>>
     % 
     %%
     % 
     % <<Fig5.PNG>>
     % 
     %%
     % 
     % <<Fig6.PNG>>
     % 
     %%
     % 
     % <<Fig7.PNG>>
     % 
     %%
     % 
     % <<Fig8.PNG>>
     % 
     %%
     % 
     % <<Fig9.PNG>>
     % 
     %%
     % 
     % <<Fig10.PNG>>
     % 
     %%
     % 
     % <<Fig11.PNG>>
     % 
     
%% Build stiffness matrix [K]
    K = diag(-element.k, 1) + ...
        diag([element.k(1), element.k(1:end-1) + element.k(2:end), element.k(end)], 0) + ...
        diag(-element.k, -1);

%% Apply boundary conditions (Convection/conductions conditions)
    K(1,1) = 1; K(1,2) = 0;
    K(end, end) = 1; K(end, end-1) = 0;
     %%
     % 
     % <<Fig12.PNG>>
     % 
%% Apply loads
    node.F_ext = zeros(node.N, 1);
    node.F_ext = node.loads; % Index denotes node
    
%% Calculate nodal solution 
    node.soln = K\node.F_ext;
     %%
     % 
     % <<Fig13.PNG>>
     % 
%% Calculate elemental solution
    element.soln = element.Ufactor'.*diff(node.soln).*element.Area;
    
%% Print nodal results
    sep = [repmat('-', 1, 40), '\n']; % Decorative separator
    fprintf(sep); disp('Nodal Displacements');
    fprintf(sep); fprintf(' Node\t\t Disp. (in)');
    node.tbl = zeros(1, 2*node.N); 
    node.tbl(1:2:end) = 1:node.N; node.tbl(2:2:end) = node.soln;
    fprintf(repmat('\n (%d)\t\t %2.4f\n', 1, node.N), node.tbl); fprintf('\n\n');
   
%% Plot nodal results
    figure, % subplot(2,1,1),
    plotCMAP([(1:node.N)', node.soln, node.soln]); grid on;
    xlabel('Node');
    ylabel('Nodal Temp.(^oF)');set(gca, 'fontweight', 'bold', 'fontsize', 11);
    
%% Print elemental results
    fprintf(sep); disp('Element Stress mean(A) (pg. 18-19)');
    fprintf(sep); fprintf(' Element\t Stress (lb*in^-2)');
    element.tbl = zeros(1, 2*element.N); 
    element.tbl(1:2:end) = 1:element.N; element.tbl(2:2:end) = element.soln;
    fprintf(repmat('\n (%d)\t\t %2.4f\n', 1, element.N), element.tbl); fprintf('\n\n');
    
%% Plot elemental results
    figure, % subplot(2,1,2), 
%    element.soln = round(element.soln, 3);
    plotCMAP([linspace(0,element.N, element.N)', element.soln, element.soln]); 
    grid on; xlabel('element'); ylabel('Heat loss');set(gca, 'fontweight', 'bold', 'fontsize', 11);
     %%
     % 
     % <<Fig14.PNG>>
     %     
     
     %%
     % Wikipedia Ref: 
     % <https://en.wikipedia.org/wiki/Dew_point>
%% Calculate dew point (HW2_P2)
    Td = @(degF, RH) degF - (9/25)*((100 - RH)); % General equation to find dew pt. temp.
    Tdp = Td(68, 40);
    disp('Dew Point Temperature: ');
    disp(Tdp); disp('Check nodal temp.`s for condensation location');
    
    %%
    % *Additional notes:*
    %%
    % 
    % * Can use a linear, quadratic, etc. fit to model the coupling coefficients
    % * Will allow for more elements/finer meshing (see tutorial from prev.
    % lab
    % * Can futher improve code to separate generalizable steps into
    % functions
    % 
    
    