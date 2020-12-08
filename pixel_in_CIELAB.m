        function pixel_in_CIELAB (fn)
            
            % linearize
            rgb = imread(fn);
            rgb1 = reshape(rgb,size(rgb,1)*size(rgb,2),3);
            
            % convert to LAB
            lab1 = rgb2lab(rgb1);
            
            % remove masked green pixels
            if 1
                mask = (rgb1(:,1) == 0) & (rgb1(:,2) == 255) & (rgb1(:,3) == 0);
                lab1_masked = lab1(mask,:);
                rgb1_masked = rgb1(mask,:);
            end
            
            % show
            scatter3(lab1(:,2),lab1(:,3),lab1(:,1),1,double(rgb1(:,:))/255,'o')
            
           
            axis equal
            xlabel('{\it a}*')
            ylabel('{\it b}*')
            zlabel('{\it L}*')
            
        end