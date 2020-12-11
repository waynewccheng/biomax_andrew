% Q: compare the histograms

% data = {};
% 
% datadsorted = {};
% 
% for i=1:8
%     for k=1:4
%         ch = ColorHistogramLAB(ct.get_filename_lab(i,k));
%         d = ch.m;
%         d1d = d(:);
%         dnonzero = d1d(d1d~=0);
%         dsorted = sort(dnonzero,'descend');
%         dcdf = cumsum(dsorted);
%         dcdfnormal = dcdf / sum(dnonzero);
%         
%         data{i,k} = dcdfnormal;
%         datasorted{i,k} = dsorted;
%     end
% end

% plot only the top 10 bins
mk = '---:';
for i=1:8
    clf
    
    %subplot(2,4,i)
    hold on
    
    for k=1:4
        y = datasorted{i,k};
        x = [1:length(y)]/length(y);
        plot(y(1:10),mk(k));
    end
    
    legend('hamamatsu','leica','zeiss','truth');
    title(sprintf('%d',i))
    
    saveas(gcf,sprintf('d%d.png',i))
end

return

% plot the bins
mk = '---:';
for i=1:8
    clf
    
    %subplot(2,4,i)
    hold on
    
    for k=1:4
        y = data{i,k};
        x = [1:length(y)]/length(y);
        plot(x,y,mk(k));
    end
    
    axis equal
    axis([0 1 0 1])
    
    legend('hamamatsu','leica','zeiss','truth');
    title(sprintf('%d',i))
    
    saveas(gcf,sprintf('%d.png',i))
end
