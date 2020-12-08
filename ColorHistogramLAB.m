% new version!! 12/5/2020
%
% color normalization paper in SPIE 2021
% color normality

% migrate to 3Scanner project


classdef ColorHistogramLAB < handle
    
    properties (Constant)
        OFFSET_AB_RANGE = 110
        OFFSET_L = 1
        OFFSET_A = ColorHistogram.OFFSET_AB_RANGE
        OFFSET_B = ColorHistogram.OFFSET_AB_RANGE
        
        SIZE_L = 100 + 1
        SIZE_A = 1 + ColorHistogram.OFFSET_AB_RANGE * 2
        SIZE_B = 1 + ColorHistogram.OFFSET_AB_RANGE * 2
    end
    
    properties
        n_pixel  % total pixel count in the image
        n_present % count of occupied bins
        n_nonzero % count of nonzero bins
        n_nonwhite % count of nonzero bins that are not white
        
        m         % 3D bins
        L
        a
        b
        
        m1        % 1D bins
        L1
        a1
        b1
        
        mLab
        mLabSorted
        mLabNonwhite
        
        mask_nonzero    % mask
        mask_nonwhite
        
        mcdf
        mcdfnormal
    end
    
    methods
        
        function [labmatL labmata labmatb] = labmat_create (obj)
            
            % create the Lab values for each bin
            % meshgrid created confusion
            
            % L*
            labmatL = repmat([0:100]',1,ColorHistogram.SIZE_A,ColorHistogram.SIZE_B);
            
            % a*
            labmata = repmat([-ColorHistogram.OFFSET_AB_RANGE:ColorHistogram.OFFSET_AB_RANGE],ColorHistogram.SIZE_L,1,ColorHistogram.SIZE_B);
            
            % b*
            vecB = zeros(1,1,ColorHistogram.SIZE_B);
            vecB(1,1,:) = [-ColorHistogram.OFFSET_AB_RANGE:ColorHistogram.OFFSET_AB_RANGE];
            labmatb = repmat(vecB,ColorHistogram.SIZE_L,ColorHistogram.SIZE_A,1);
           
            obj.L = labmatL;
            obj.a = labmata;
            obj.b = labmatb;
            
            obj.L1 = labmatL(:);
            obj.a1 = labmata(:);
            obj.b1 = labmatb(:);

        end
        
        function [normality m1_ratio m2_ratio] = color_normality (obj1, obj2, threshold)
            % Calculate the "color normality" of two images after color normalization
            % inputs: two ColorHistogram m1 and m2, minimum pixel count required to be considered
            % outputs: color_normality, ratio of intersection to m1_gamut, ratio of intersection to m2_gamut
            
            m1 = obj1.m;
            m2 = obj2.m;
            
            m_union = (m1>threshold) | (m2>threshold);
            m_intersect = (m1>threshold) & (m2>threshold);
            n_union = nnz(m_union);
            n_intersect = nnz(m_intersect);
            
            m1_ratio = n_intersect / obj1.n_present;
            m2_ratio = n_intersect / obj2.n_present;
            normality = m1_ratio * m2_ratio;
            
            %[obj1.n_present obj2.n_present n_union n_intersect]
            %[m1_ratio m2_ratio]
        end
        
        function m1_ratio = diff (obj1, obj2, threshold)
            % similar to color_normality, but returns only m1 data
            m1 = obj1.m;
            m2 = obj2.m;
            m_union = (m1>threshold) | (m2>threshold);
            m_intersect = (m1>threshold) & (m2>threshold);
            n_union = nnz(m_union);
            n_intersect = nnz(m_intersect);
            m1_ratio = n_intersect / obj1.n_present;
            m2_ratio = n_intersect / obj2.n_present;
            
            %[obj1.n_present obj2.n_present n_union n_intersect]
            %[m1_ratio m2_ratio]
        end
        
        function obj = ColorHistogramLAB (lab)

            obj.labmat_create;
            
            %
            % input is LAB
            %
            
            % 3D to 1D
            lab1 = reshape(lab,size(lab,1)*size(lab,2),3);
            
            % assign total pixel count in the image
            obj.n_pixel = size(lab1,1);
            
            % double to integer
            lab1 = floor(lab1);
            
            % choose one plane of L*
            for l = 0:100       % l: layer in L*
                
                % mask for each layer in L*
                lab_masked = lab1(lab1(:,1)==l,:);
                
                % convert a* to matrix index
                lab_a = lab_masked(:,2) + obj.OFFSET_A + 1;
                
                % convert b* to matrix index
                lab_b = lab_masked(:,3) + obj.OFFSET_B + 1;
                
                % use an ONE matrix to do the counting
                many_ones = ones(size(lab_masked,1),1);
                
                % make sure evreything is in range
                assert(nnz(lab_a < 1)==0);
                assert(nnz(lab_a > obj.SIZE_A)==0);
                assert(nnz(lab_b < 1)==0);
                assert(nnz(lab_b > obj.SIZE_B)==0);
                
                %
                % make a sparse matrix with a* and b* and ones
                %
                % Accumulate Values into Sparse Matrix
                % https://www.mathworks.com/help/matlab/ref/sparse.html
                s{l+1} = sparse(lab_a,lab_b,many_ones,obj.SIZE_A,obj.SIZE_B);
                
            end
            
            %
            % convert back to full 3D matrix
            %
            obj.m = zeros(obj.SIZE_L,obj.SIZE_A,obj.SIZE_B);
            for l = 0:100
                layer = l+1;
                
                %
                % convert sparse matrix to full matrix
                %
                obj.m(layer,:,:) = full(s{layer});
            end
            obj.m1 = obj.m(:);
            
            %
            % assign count of occupied bins
            %
            obj.n_present = nnz(obj.m);
            
            obj.mLab = [obj.m1 obj.L1 obj.a1 obj.b1];
            
            obj.mLabSorted = sortrows(obj.mLab,'descend');

            obj.mask_nonzero = obj.mLabSorted(:,1)>0;
            
            mask_nonwhiteLAB = obj.remove_white(obj.mLabSorted(:,2:4));
            
            obj.mask_nonwhite = mask_nonwhiteLAB & obj.mask_nonzero;

            obj.n_nonzero = nnz(obj.mask_nonzero);
            obj.n_nonwhite = nnz(obj.mask_nonwhite);
            
            obj.mLabNonwhite = obj.mLabSorted(obj.mask_nonwhite,:);
            
            obj.mcdf = cumsum(obj.mLabNonwhite(:,1));
            obj.mcdfnormal = obj.mcdf/sum(obj.mLabNonwhite(:,1));
            
        end
        
        function mask = remove_white (obj,lab)
            %
            % remove white
            %
            % threshold from white
            chroma_th = 10;
            
            % calculate chroma
            chroma = (lab(:,2).^2 + lab(:,3).^2) .^ 0.5;
            
            % filter
            mask = chroma > chroma_th;
        end
        
        function obj = ColorHistogramRGB (fn)
            %
            % input is an image file
            %
            
            % 3D to 1D
            rgb = imread(fn);
            rgb1 = reshape(rgb,size(rgb,1)*size(rgb,2),3);
            
            % assign total pixel count in the image
            obj.n_pixel = size(rgb1,1);
            
            % RGB to LAB
            lab1 = rgb2lab(rgb1);
            
            % double to integer
            lab1 = floor(lab1);
            
            % choose one plane of L*
            for l = 0:100       % l: layer in L*
                
                % mask for each layer in L*
                lab_masked = lab1(lab1(:,1)==l,:);
                
                % convert a* to matrix index
                lab_a = lab_masked(:,2) + obj.OFFSET_A;
                
                % convert b* to matrix index
                lab_b = lab_masked(:,3) + obj.OFFSET_B;
                
                % use an ONE matrix to do the counting
                many_ones = ones(size(lab_masked,1),1);
                
                % make sure evreything is in range
                assert(nnz(lab_a < 1)==0);
                assert(nnz(lab_a > obj.SIZE_A)==0);
                assert(nnz(lab_b < 1)==0);
                assert(nnz(lab_b > obj.SIZE_B)==0);
                
                %
                % make a sparse matrix with a* and b* and ones
                %
                % Accumulate Values into Sparse Matrix
                % https://www.mathworks.com/help/matlab/ref/sparse.html
                s{l+1} = sparse(lab_a,lab_b,many_ones,obj.SIZE_A,obj.SIZE_B);
                
            end
            
            % convert back to full 3D matrix
            obj.m = zeros(obj.SIZE_L,obj.SIZE_A,obj.SIZE_B);
            for l = 0:100
                layer = l+1;
                
                %
                % convert sparse matrix to full matrix
                %
                obj.m(layer,:,:) = full(s{layer});
                
            end
            
            %
            % assign count of occupied bins
            %
            obj.n_present = nnz(obj.m);
            
        end
        
    end
    
end

