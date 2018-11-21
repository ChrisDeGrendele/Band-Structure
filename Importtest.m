uiimport('Plotdata.mat') 
figure;
ax1=subplot('Position',[.4 .3 .3 .3]); %Position = [Left bottom width height];
scatter(abs(Real_k), Real_E, '.')
title('Real k v. E')
ylim([a,b])

hold on
ax2=subplot('Position',[.1 .3 .3 .3]);
set(ax2,'YTick',[],'XTick',[]);
scatter( -abs(Imag_k), Imag_E, '.')
title('Imaginary k v. E')
ylim([a,b])