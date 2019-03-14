function  plot_temperature( x, y, data, mesh,M)
%fun_fem2_IGA_Postprocessing_plotDispStress 
figure;
title(['µÚ', num2str(M), ' Ìì'])
hold on;
colormap jet;
axis off
axis equal
axis tight
set(gcf,'color','white')
colorbar;
colorbar('FontSize',12);
patch('Faces', mesh, 'Vertices', [x ,y], 'CData', data, 'FaceColor', 'interp', 'EdgeColor', 'none');
end

