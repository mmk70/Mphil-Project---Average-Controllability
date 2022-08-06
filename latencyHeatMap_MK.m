function latencyHeatMap_MK(metric,channels,fileID) 
%Input = first_spike, channels


pltIndex=[find(channels==21),find(channels==31),find(channels==41),... %count across columns (subplot index plots across columns, e.g. sublot (4,4,2) the plots in column 2 not row 2
        find(channels==51),find(channels==61),find(channels==71),find(channels==12),...
        find(channels==22),find(channels==32),find(channels==42),find(channels==52),...
        find(channels==62),find(channels==72),find(channels==82),find(channels==13),...
        find(channels==23),find(channels==33),find(channels==43),find(channels==53),...
        find(channels==63),find(channels==73),find(channels==83),find(channels==14),...
        find(channels==24),find(channels==34),find(channels==44),find(channels==54),...
        find(channels==64),find(channels==74),find(channels==84),find(channels==15),...
        find(channels==25),find(channels==35),find(channels==45),find(channels==55),...
        find(channels==65),find(channels==75),find(channels==85),find(channels==16),...
        find(channels==26),find(channels==36),find(channels==46),find(channels==56),...
        find(channels==66),find(channels==76),find(channels==86),find(channels==17),...
        find(channels==27),find(channels==37),find(channels==47),find(channels==57),...
        find(channels==67),find(channels==77),find(channels==87),find(channels==28),...
        find(channels==38),find(channels==48),find(channels==58),find(channels==68),...
        find(channels==78)];
metric = metric(pltIndex)
heatMatrix = zeros(8, 8); 
    heatMatrix(1, 1) = NaN; 
    heatMatrix(1, 8) = NaN; 
    heatMatrix(8, 1) = NaN;
    heatMatrix(8, 8) = NaN; 
heatMatrix(2:7) = metric(1:6);
    heatMatrix(9:56) = metric(7:54); 
    heatMatrix(58:63) = metric(55:60);
    heatMatrix(33) = NaN %channel 15 
heatMatrix = heatMatrix' 
    figure
    h =  imagesc(heatMatrix); 
    set(h, 'AlphaData', ~isnan(heatMatrix))
    text(0.9,33,'R','FontName','Arial','FontSize',18)
   % set(gcf,'Visible', 'off')
   
    c = spring;
c = flipud(c);
colormap(c);


cb = colorbar;
        ylabel(cb, 'AP latency post stimulation (ms)')
   
   cb.TickDirection = 'out';
    cb.Location = 'Southoutside';
    cb.Box = 'off';
  
 

     % make it square-ish
    set(gcf, 'Position', [100, 100, 800, 700])
    
    % set font size 
    set(gca, 'FontSize', 14)

    
        title({strcat(fileID(1:end-4),' Latency MEA grid'),' '});
   
end 
%make a seperate one for latency (opposite heatmap)