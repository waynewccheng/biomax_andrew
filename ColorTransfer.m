

classdef ColorTransfer
    % Analysis for Andrew's 3-scanner comparison study
    % WCC 12/4/2020
    
    properties
        dname = {'Bladder_M13','Brain_H10','Breast_A1','Colon_H6','Kidney_H7','Liver_H9','Lung_J7','UterineCervix_B10'};
        sname = {'hamamatsu','leica','zeiss','truth'};
    end
    
    methods
        
        function obj = ColorTransfer
            mkdir findings
        end
        
        function [labscan labtruth] = get_lab_data (obj,slide_no,scanner_no)
            
            % input is a .mat file
            if scanner_no < 4                  % one of the 3 scanners
                
                % construct filename
                tn = ['BiomaxOrgan10_' obj.dname{slide_no}];
                sn = obj.sname{scanner_no};
                fn = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\' sn '\400 sRGB\final_images.mat'];
                
                load(fn,'LAB_scan');
                labscan = LAB_scan;
                
                fn = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\' sn '\400 sRGB\final_images.mat'];
                load(fn,'LAB_truth_reg_trimmed');
                labtruth = LAB_truth_reg_trimmed;
                
            end
            
        end
        
        function lab = get_filename_lab (obj,slide_no,scanner_no)
            
            % input is a .mat file
            if scanner_no < 4                  % one of the 3 scanners
                
                % construct filename
                tn = ['BiomaxOrgan10_' obj.dname{slide_no}];
                sn = obj.sname{scanner_no};
                fn = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\' sn '\400 sRGB\final_images.mat'];
                
                load(fn,'LAB_scan');
                lab = LAB_scan;
                
            else                               % the truth
                
                % construct filename
                tn = ['BiomaxOrgan10_' obj.dname{slide_no}];
                
                % use one of the scanners!!
                sn = obj.sname{1};
                
                fn = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\' sn '\400 sRGB\final_images.mat'];
                load(fn,'LAB_truth_reg_trimmed');
                lab = LAB_truth_reg_trimmed;
                
            end
            
        end
        
        function pixel_in_CIELAB_slide_scanner (obj,slide_no,scanner_no)
            
            lab = obj.get_filename_lab (slide_no,scanner_no);
            
            obj.LABpixel_in_CIELAB(lab);
            
        end
        
        
        function LABlist_in_CIELAB_size (obj, labs1)
            
            lab1 = labs1(:,1:3);
            sz1 = labs1(:,4);
            
            % convert to sRGB
            rgb1 = lab2rgb(lab1);
            
            % show
            scatter3(lab1(:,2),lab1(:,3),lab1(:,1),sz1,double(rgb1(:,:)),'o','filled','MarkerEdgeColor','k')
            
            axis equal
            
            xlabel('{\it a}*')
            ylabel('{\it b}*')
            zlabel('{\it L}*')
            
            view(30,18)
            axis equal
        end
        
        function LABlist_in_CIELAB (obj, lab1)
            
            % convert to sRGB
            rgb1 = lab2rgb(lab1);
            
            % show
            scatter3(lab1(:,2),lab1(:,3),lab1(:,1),50,double(rgb1(:,:)),'o','filled')
            
            axis equal
            
            axis([-100 100 -100 100 0 100])
            
        end
        
        function LABpixel_in_CIELAB (obj, lab)
            
            % linearize
            lab1 = reshape(lab,size(lab,1)*size(lab,2),3);
            
            % convert to sRGB
            rgb1 = lab2rgb(lab1);
            
            % show
            scatter3(lab1(:,2),lab1(:,3),lab1(:,1),1,double(rgb1(:,:)),'o')
            
            axis equal
            
            axis([-40 100 -50 50 0 100])
            
        end
        
        
        function dE_quiver (obj)
            % Q: what is the color shift?
            % A: use quiver3 to show the difference
            % WCC 12/4/2020
            
            for i = 1:8
                for k = 1:3
                    tn = ['BiomaxOrgan10_' obj.dname{i}];
                    sn = obj.sname{k};
                    
                    %         if i==8 & k==3
                    %             continue
                    %         end
                    
                    fn_scan = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\' sn '\400 sRGB\final_images.mat']
                    fn_truth = ['C:\Users\wcc\Desktop\biomax_andrew\' tn '\truth\900 sRGB\truth.tif']
                    
                    
                    load(fn_dE)
                    
                    step = 25;
                    
                    LAB_scan_small = LAB_scan(1:step:end,1:step:end,:);
                    LAB_truth_reg_trimmed_small = LAB_truth_reg_trimmed(1:step:end,1:step:end,:);
                    
                    im_diff = LAB_scan - LAB_truth_reg_trimmed;
                    im_diff_small = LAB_scan_small - LAB_truth_reg_trimmed_small;
                    v = im_diff_small;
                    
                    im_dE76 = sum(im_diff.^2,3).^0.5;
                    
                    clf
                    subplot(4,4,1)
                    im = imread(fn_truth);
                    image(im)
                    axis image
                    axis off
                    title('truth')
                    
                    subplot(4,4,2)
                    load(fn_scan,'unregistered_trimmed');
                    image(unregistered_trimmed)
                    axis image
                    axis off
                    title(sn)
                    
                    subplot(4,4,3)
                    imagesc(im_dE76)
                    axis image
                    axis off
                    colorbar
                    title('CIE {\Delta}{\itE}^*_{76}')
                    
                    subplot(4,4,5)
                    imagesc(im_diff(:,:,1))
                    axis image
                    axis off
                    colorbar
                    title('CIE {\Delta}{\itL}^*_{76}')
                    
                    subplot(4,4,6)
                    imagesc(im_diff(:,:,2))
                    axis image
                    axis off
                    colorbar
                    title('CIE {\Delta}{\ita}^*_{76}')
                    
                    subplot(4,4,7)
                    imagesc(im_diff(:,:,3))
                    axis image
                    axis off
                    colorbar
                    title('CIE {\Delta}{\itb}^*_{76}')
                    
                    subplot(4,4,[9 10 13 14])
                    quiver3(LAB_truth_reg_trimmed_small(:,:,2),LAB_truth_reg_trimmed_small(:,:,3),LAB_truth_reg_trimmed_small(:,:,1),v(:,:,2),v(:,:,3),v(:,:,1))
                    
                    xlabel('CIE {\Delta}{\ita}^*_{76}')
                    ylabel('CIE {\Delta}{\itb}^*_{76}')
                    zlabel('CIE {\Delta}{\itL}^*_{76}')
                    axis equal
                    title([tn ' : ' sn],'Interpreter','None')
                    view(34,25)
                    
                    subplot(4,4,[11 12 15 16])
                    quiver3(ones(size(v,1),size(v,2)),v(:,:,2),v(:,:,3),v(:,:,1))
                    
                    xlabel('CIE {\Delta}{\ita}^*_{76}')
                    ylabel('CIE {\Delta}{\itb}^*_{76}')
                    zlabel('CIE {\Delta}{\itL}^*_{76}')
                    axis equal
                    title([tn ' : ' sn],'Interpreter','None')
                    view(12,40)
                    
                    
                    saveas(gcf,sprintf('findings/%s_%s.png',obj.dname{i},sn))
                end
            end
        end
        
        function cdfplot_8 (obj)
            % Q: Analysis of Andrew's data from 3 scanners
            % A: use the PV20 style
            % WCC 10/30/2020
            
            clf
            for k = 1:8
                
                dn = obj.dname{k};
                fn = ['BiomaxOrgan10_' dn '/dE/'];
                
                % load the data
                load([fn 'dE_hamamatsu.mat'])
                dE1_hamamatsu = dE00(:);
                load([fn 'dE_leica.mat'])
                dE1_leica = dE00(:);
                load([fn 'dE_zeiss.mat'])
                dE1_zeiss = dE00(:);
                load([fn 'dE.mat'])
                dE1_mono = dE00(:);
                
                subplot(2,4,k)
                % retrieve the histogram data
                %clf
                h = histogram(dE1_hamamatsu,[0:100],'Normalization','cdf');
                H{1} = h.Values;
                h = histogram(dE1_leica,[0:100],'Normalization','cdf');
                H{2} = h.Values;
                h = histogram(dE1_zeiss,[0:100],'Normalization','cdf');
                H{3} = h.Values;
                h = histogram(dE1_mono,[0:100],'Normalization','cdf');
                H{4} = h.Values;
                
                % plot the CDF
                plot(H{1},'b-')
                hold on
                plot(H{2},'g:')
                plot(H{3},'r--')
                plot(H{4},'k:')
                xlabel('{\Delta}E_{00}')
                ylabel('Cumulative Density')
                legend('Hamamatsu','Leica','Zeiss','Mono')
                legend('Location','southeast')
                axis square
                axis([0 30 0 1])
                title(dn,'Interpreter','None')
                
            end
            
            saveas(gcf,sprintf('findings/cdfplot_8.png'))
            
        end
        
        function roi_8x4 (obj)
            % Q: Analysis of Andrew's data from 3 scanners
            % A: use the PV20 style
            % WCC 10/30/2020
            
            clf
            for i = 1:8
                
                dn = obj.dname{i};
                
                for k = 1:3
                    j = 4*(i-1) + k;
                    subplot(8,4,j)
                    sn = obj.sname{k};
                    fn = ['BiomaxOrgan10_' dn '/' sn '/400 sRGB/final_images.mat'];
                    load(fn,'unregistered_trimmed');
                    image(unregistered_trimmed)
                    axis image
                    axis off
                    title(sn)
                end
                
                for k=4
                    j = 4*(i-1) + k;
                    subplot(8,4,j)
                    sn = obj.sname{k};
                    fn = ['BiomaxOrgan10_' dn '/' sn '/900 sRGB/truth.tif'];
                    im = imread(fn);
                    image(im)
                    axis image
                    axis off
                    title(sn)
                end
                
            end
            
            saveas(gcf,sprintf('findings/roi_8x4.png'))
            
        end
        
    end
    
    
    methods (Static)
        
        function plot_color_gamut_distribution

            ct = ColorTransfer;
            
            for i = 1:8
                for k = 1:4
                    subplot(8,4,k+4*(i-1))
                    ct.pixel_in_CIELAB_slide_scanner(i,k)
                    title([ct.dname{i} ' : ' ct.sname{k}],'Interpreter','none')
                end
            end
        end
        
        
        function print_color_gamut_histogram

            ct = ColorTransfer;
            
            for i = 1:1
                
                clf
                
                for k = 1:4
                    
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
                    n_nonzero = size(dLabnonzero,1)
                    
                    %
                    % remove white
                    %
                    % threshold from white
                    chroma_th = 10;
                    
                    % calculate chroma
                    chroma = (dLabnonzero(:,3).^2 + dLabnonzero(:,4).^2) .^ 0.5;
                    
                    % filter
                    mask = chroma > chroma_th;
                    
                    dLabnonwhite = dLabnonzero(mask,:);
                    n_nonwhite = size(dLabnonwhite,1)
                    
                    % only the first ones
                    n = nnz(dLabnonwhite(:,1)>0.001)
                    
                    dLabfront = dLabnonwhite(1:n,:);
                    
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
                    title([ct.dname{i} ' : ' ct.sname{k}],'Interpreter','none')
                    
                    axis equal
                    %axis([-100 100 -100 100 0 100])
                    view(30,18)
                    
                end
                
                saveas(gcf,sprintf('hist%d.png',i))

            end
        end
        
    end
    
end

