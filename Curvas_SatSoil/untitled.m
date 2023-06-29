figure()
tcl = tiledlayout(2,2);
nexttile(tcl)
line1 = plot(1:10,rand(1,10),'b','DisplayName','Data Axes 1');
title('Axes 1'); 
nexttile(tcl)
line2 = plot(1:10,rand(1,10),'g','DisplayName','Data Axes 2');
title('Axes 2');
nexttile(tcl)
line3 = plot(1:10,rand(1,10),'r','DisplayName','Data Axes 3');
title('Axes 3');
nexttile(tcl)
line4 = plot(1:10,rand(1,10),'c','DisplayName','Data Axes 4');
title('Axes 4');
% Construct a Legend with the data from the sub-plots
hL = legend([line1,line2,line3,line4]); 
% Move the legend to the right side of the figure
hL.Layout.Tile = 'East';