function im = plotLandscape(DEM,drca,mapping,threshold_area,num)

%%%%%%%%%%% LIRAN GOREN, gorenl@bgu.ac.il, 20/12/2023 %%%%%%%%%%%%%%%%

figure(num)
axis equal
hold on

max_a = max(drca(:,4));
for i = 1:length(drca(:,1))
    if drca(i,4)>threshold_area
        [xd,yd] = ind2coord(DEM,drca(i,1));
        rec_line = mapping(drca(i,2)); %finding the line of the receiver
        [xr,yr] = ind2coord(DEM,drca(rec_line,1));
        a = drca(i,4);
        plot([xd xr],[yd yr],'b','LineWidth',(a/max_a)*5);

    end
end
set(gcf, 'color', 'white');   
axis off
im = gcf;

