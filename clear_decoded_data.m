function decodedDataCell = clear_decoded_data(decodedDataCell)
    % calculate speed based on position/time
    count = 1;
    index_to_remove = [];
    for i = 1:size(decodedDataCell, 1)-1
       if (decodedDataCell{i,1}.plane == decodedDataCell{i+1,1}.plane)
            airplane_loc1 = lla2ecef([decodedDataCell{i,1}.lat, ...
                   decodedDataCell{i,1}.lon, decodedDataCell{i,1}.alt]);
            airplane_loc2 = lla2ecef([decodedDataCell{i+1,1}.lat, ...
                   decodedDataCell{i+1,1}.lon, decodedDataCell{i+1,1}.alt]);

            distance = sqrt(sum((airplane_loc2-airplane_loc1).^2));
            start_time = max(decodedDataCell{i,1}.even_cprtime, decodedDataCell{i,1}.odd_cprtime);
            end_time = max(decodedDataCell{i+1,1}.even_cprtime, decodedDataCell{i+1,1}.odd_cprtime);
            speed = distance/(end_time-start_time)*1e7;
            if (speed > 250)
                index_to_remove(count) = i+1;
                count = count+1;
            end
       end    
    end

    % remove data that gives wrong speed
    if (~isempty(index_to_remove))
            decodedDataCell([index_to_remove], :) = [];
    end

end

