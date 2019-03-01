function FEM_1D_TUT
    close all; clc;
    
    %% Define No. of elements/nodes
    element.N = 4; % No. of elements
    node.N = element.N + 1; % No. of nodes

    %% Given constants
    LMAX = 10; % Total length in y-dir (in)
    w = [2, 1]; % Element widths @ x = [0, LMAX] (in)
    E = 10.4e6; % Young's modulus (lb*in^-2)
    t = 0.125; % Thickness (in)

    gbl.y = linspace(0, LMAX, node.N); % Global position 0<= gbl.y <= L
    element.L = max(gbl.y)/element.N; % Element lengths
    Func.Shape = polyfit([0, LMAX], w, 1); % Linear fit of cross-sec area (coefficients)
    element.Area = polyval(Func.Shape, gbl.y).*t; % A(y) -- (in^2)
    element.Abar = (element.Area(1:end-1) + element.Area(2:end))./2; % Mean cross-sec. area (in^2)
    element.Modulus = ones(1,element.N).*E; % Element Modulus of Elasticity (lb*in^-2)
    element.k = element.Modulus.*element.Abar./element.L; % Element effective stiffness (lb*in^-1)

    %% Build stiffness matrix [K]
    K = diag(-element.k, 1) + ...
        diag([element.k(1), element.k(1:end-1) + element.k(2:end), element.k(end)], 0) + ...
        diag(-element.k, -1);

    %% Apply boundary conditions and forces
    % boundary conditions
    K(1,1) = 1; K(1,2) = 0;
    % forces
    node.F_ext = zeros(node.N, 1);
    node.F_ext(end) = 1000; % Force on last node (lb)
    %node.F_ext(floor(node.N/2)) = -10000; % Can add extra external forces
    
    %% Calculate nodal displacement and elemental stresses
    node.disp = K\node.F_ext;
    element.stress = element.Modulus'.*diff(node.disp)./element.L;
    
    %% Print/plot nodal results
    sep = [repmat('-', 1, 40), '\n']; % Decorative separator
    fprintf(sep); disp('Nodal Displacements');
    fprintf(sep); fprintf(' Node\t\t Disp. (in)');
    node.tbl = zeros(1, 2*node.N); 
    node.tbl(1:2:end) = 1:node.N; node.tbl(2:2:end) = node.disp;
    fprintf(repmat('\n (%d)\t\t %2.6e\n', 1, node.N), node.tbl); fprintf('\n\n');
   
    figure, 
    subplot(2,1,1), plotCMAP([gbl.y', node.disp, node.disp]); grid on;
    %set(gca, 'xticklabel', ''); % Turn off x-axis tick labels
    ylabel('Nodal Disp. (in)');
    
    %% Print/plot elemental stresses
    sep = [repmat('-', 1, 40), '\n']; % Decorative separator
    fprintf(sep); disp('Element Stress mean(A) (pg. 18-19)');
    fprintf(sep); fprintf(' Element\t Stress (lb*in^-2)');
    element.tbl = zeros(1, 2*element.N); 
    element.tbl(1:2:end) = 1:element.N; element.tbl(2:2:end) = element.stress;
    fprintf(repmat('\n (%d)\t\t %2.6e\n', 1, element.N), element.tbl); fprintf('\n\n');
    
    subplot(2,1,2), 
    plotCMAP([linspace(0,element.L*element.N, element.N)', element.stress, element.stress]); 
    grid on; xlabel('y-position (in)'); ylabel('Stress (lb*in^2)');
    
    function plotCMAP(data)

        x = data(:,1) ; %// extract "X" column
        y = data(:,2) ; %// same for "Y"
        c = data(:,3) ; %// extract color index for the custom colormap

        %% // Prepare matrix data
        xx = [x x];           %// create a 2D matrix based on "X" column
        yy = [y y];           %// same for Y
        zz = zeros(size(xx)); %// everything in the Z=0 plane
        cc = [c c] ;         %// matrix for "CData"

        %// draw the surface (actually a line)
        hs = surf(xx,yy,zz,cc,'EdgeColor','interp','FaceColor','none','Marker','o') ;

        colormap('jet') ;     %// assign the colormap
        shading flat                    %// so each line segment has a plain color
        view(2) %// view(0,90)          %// set view in X-Y plane
        colorbar