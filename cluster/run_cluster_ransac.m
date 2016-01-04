function run_cluster_ransac(id_start)
    % cd /n/fs/sun3d/DDD/TripleD/cluster
    % /n/fs/vision/ionicNew/starter.sh run_cluster_ransac 1000mb 1:00:00 1 100 1
    sequenceName = 'hotel_umd/maryland_hotel3/';
    write2path = fullfile('/n/fs/sun3d/DDD/matchdata/',sequenceName);
     
    load(fullfile(write2path,'matchpair.mat'));
    for id = id_start:100:length(frag_i)
           
        filename = [sequenceName sprintf('%d_%d-%d_%d',...
                            fragments(frag_i(id),:)-1,fragments(frag_j(id),:)-1)];
        %hotel_umd/maryland_hotel3/800_849-850_899_scores2.txt
        system(horzcat(['unset LD_LIBRARY_PATH; ','./cluster ',filename]));
    end
end