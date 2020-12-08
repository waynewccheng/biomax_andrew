i = 8;
k = 2;
[labscan labtruth] = ct.get_lab_data(i,k);
labscan1 = reshape(labscan,size(labscan,1)*size(labscan,2),3);
labtruth1 = reshape(labtruth,size(labtruth,1)*size(labtruth,2),3);

        OFFSET_AB_RANGE = 110
        OFFSET_L = 1
        OFFSET_A = OFFSET_AB_RANGE
        OFFSET_B = OFFSET_AB_RANGE
        
labindex = floor(labtruth1);
% labindex(:,1) = labindex(:,1) + OFFSET_L;
% labindex(:,2) = labindex(:,2) + OFFSET_A;
% labindex(:,3) = labindex(:,3) + OFFSET_B;

q = [labtruth1 labscan1 labindex];

targetindex = [99 0 -1]
targetindex = chdata{i,4}.mLabNonwhite(100,2:4)

mask = q(:,7)==targetindex(1) & q(:,8)==targetindex(2) & q(:,9)==targetindex(3);
nnz(mask)

q2 = q(mask,:);
qdiff = q2(:,4:6) - q2(:,1:3);
qdE = sum(qdiff.^2,2).^0.5;

clf
hold on
plot(q2(:,1),q2(:,4),'.')
plot(q2(:,2),q2(:,5),'.')
plot(q2(:,3),q2(:,6),'.')

clf
%plot3(qdiff(:,1),qdiff(:,2),qdiff(:,3),'.')

%n = size(qdiff,1);
%quiver3(zeros(n,1),zeros(n,1),zeros(n,1),qdiff(:,1),qdiff(:,2),qdiff(:,3))

quiver3(q2(:,1),q2(:,2),q2(:,3),qdiff(:,1),qdiff(:,2),qdiff(:,3))

grid on



