function pltSurf1(x, y, c)
%     x = data(:,1) ; %// extract "X" column
%     y = data(:,2) ; %// same for "Y"
%     c = data(:,3) ; %// extract color index for the custom colormap
    
    %// draw the surface (actually a line)
    surf([x, x], [y, y], zeros(size([x, x])), [c, c], ...
        'EdgeColor', 'interp', ...
        'FaceColor', 'none', ...
        'Marker','o'); view(2); 
    colormap('jet'); 
    shading flat; 
    colorbar
end