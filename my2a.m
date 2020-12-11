% Q: compare the histograms after remove the background

chdata = {};

for i=1:8
    for k=1:4
        ch = ColorHistogramLAB(ct.get_filename_lab(i,k));
        chdata{i,k} = ch;
    end
end


% mk = '---:';
% for i=1:8
%     clf
%     
%     %subplot(2,4,i)
%     hold on
%     
%     for k=1:4
%         y = datasorted{i,k};
%         x = [1:length(y)]/length(y);
%         plot(y(1:10),mk(k));
%     end
%     
%     legend('hamamatsu','leica','zeiss','truth');
%     title(sprintf('%d',i))
%     
%     saveas(gcf,sprintf('d%d.png',i))
% end
% 
% return

mk = '---:';
clf
for i=1:8
    
    subplot(2,4,i)
    hold on
    
    for k=1:4
        ch = chdata{i,k};
        y = ch.mcdfnormal;
        x = [1:ch.n_nonwhite]/ch.n_nonwhite;
        plot(x,y,mk(k));
    end
    
    axis equal
    axis([0 1 0 1])
    
    legend('hamamatsu','leica','zeiss','truth');
    legend('Location','southeast')
    title(sprintf('%d',i))
    
%    saveas(gcf,sprintf('%d.png',i))
end
saveas(gcf,'1_8.png')
