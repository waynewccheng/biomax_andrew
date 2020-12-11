% Q: which color in CIELAB is the most populated?
% A: generate the 3D bubble plot

for i=1:8
    clf
    for k=1:4
        subplot(1,4,k)
        ch = ColorHistogramLAB(ct.get_filename_lab(i,k));
        
        % the 3D histogram matrix
        d = ch.m;
        [L a b] = ch.labmat_create;
        
        % number of pixels
        dsum = sum(d,'all');
        
        % normalization
        dnormal = d/dsum;
        
        % link histogram with a, b, and L
        d1 = dnormal(:);
        a1 = a(:);
        b1 = b(:);
        L1 = L(:);
        dLab = [d1 L1 a1 b1];
        
        % sort by #
        dLabsorted = sortrows(dLab,'descend');
        
        % non-zero only
        dLabnonzero = dLabsorted(dLabsorted(:,1)>0,:);
        
        % only the first ones
        n = size(dLabnonzero,1)
        n = 300
        n = nnz(dLabsorted(:,1)>0.001)
        
        dLabfront = dLabnonzero(1:n,:);
        
        % calculate size
        sz = dLabfront(:,1);
        sz2 = log10(log10(sz)+10);
        sz2min = min(sz2);
        sz2max = max(sz2);
        sz3 = 1 + 80/(sz2max-sz2min)*(sz2-sz2min);
        
        % add size
        LabS = [dLabfront(:,2:4) sz3];
        
        % pass only L, a, b, and S
        ct.LABlist_in_CIELAB_size(LabS)
        
        xlabel('{\it a}*')
        ylabel('{\it b}*')
        zlabel('{\it L}*')
        
        axis equal
        %axis([-100 100 -100 100 0 100])
        view(30,18)
        
    end
    saveas(gcf,sprintf('hist%d.png',i))
end
