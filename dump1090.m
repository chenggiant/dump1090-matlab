% Main function for decoding ADSB, save decoded info in decodedDataCell

function [decodedDataCell] = dump1090(seconds,dataPath,filename,sampRate)       

%     seconds = 2;
%     dataPath = 'D:\USRP_data\ADSB';
%     filename = 'adsb_feb13_id5_rx1_set3.dat';
%     sampRate = 0.4e7;
    
    step = 0.5;  % processing 0.5 second of data every step

    packetDataCell = [];
    planeIndex = 1;
    count = 1;
    planeList = [''];
    decodedDataCell = [];

    for skipSeconds = 0:step:seconds-step
        packetDataCell = [packetDataCell; decode_adsb_signal(skipSeconds,filename,dataPath,step,sampRate)];
    end

    % Saving plane ID in planeList
    for i = 1:size(packetDataCell,1)
        if ~ismember(packetDataCell{i,1}.planeId, planeList,'rows')
            planeList(planeIndex, :) = packetDataCell{i,1}.planeId;
            planeIndex = planeIndex + 1;
        end
    end

    % Decoding packets

    % a is to store packet data with position info
    a = struct;
    a.odd_cprlat = [];
    a.odd_cprlon = [];
    a.odd_cprtime = [];
    a.even_cprlat = [];
    a.even_cprlon = [];
    a.even_cprtime = [];
    a.alt = [];

    % b is to store packet data with velocity info
    b = struct;
    b.ew_dir = [];
    b.ew_velocity = [];
    b.ns_dir = [];
    b.ns_velocity = [];
    b.vert_rate_source = [];
    b.vert_rate_sign = [];
    b.vert_rate = [];
    b.cprtime = [];
    b.accuracy=[];

    % c is to store packet data with velocity info
    c = struct;
    c.hori_velocity = [];
    c.heading = [];
    c.vert_rate = [];
    c.vert_rate_sign = [];
    c.cprtime = [];

    % Go through all the DF 17 packets, decode Lat, Lon, Alt infomation 
    % and save it to decodedDataCell for further process
    for planeIndex = 1:size(planeList,1)
        for packetIndex = 1:size(packetDataCell,1)
            if ((packetDataCell{packetIndex,:}.DF == 17) && ...
                    isequal(packetDataCell{packetIndex,:}.planeId, planeList(planeIndex,:)))
                bits = packetDataCell{packetIndex,:}.bits;
                startIndex = packetDataCell{packetIndex,:}.startIndex;
                message = bits(33:33+56-1);
                if ((vec2dec(message(1:5)) >= 9) && (vec2dec(message(1:5)) <= 18))
                    fflag = message(22);  % cpr format bit, odd or even
                    raw_latitude = vec2dec(message(23:23+17-1));
                    raw_longitude = vec2dec(message(40:40+17-1));
                    if message(20)
                        raw_altitude = vec2dec(message(9:9+11-1));
                        a(planeIndex).alt = 25*raw_altitude - 1000;
                        a(planeIndex).alt = a(planeIndex).alt/3.2808; % convert feet to meter
                    else
                        a(planeIndex).alt = 0;
                    end
                    if fflag
                        a(planeIndex).odd_cprlat = raw_latitude;
                        a(planeIndex).odd_cprlon = raw_longitude;
                        a(planeIndex).odd_cprtime = startIndex;
                    else
                        a(planeIndex).even_cprlat = raw_latitude;
                        a(planeIndex).even_cprlon = raw_longitude;
                        a(planeIndex).even_cprtime = startIndex;
                    end
                    if (~isempty(a(planeIndex).odd_cprtime) && ~isempty(a(planeIndex).even_cprtime))
                        res = decode_cpr(a(planeIndex));
                        res.plane = planeList(planeIndex,:);
                        if (a(planeIndex).alt ~= 0)
                            decodedDataCell{count,1} = res;
                            count = count + 1; 
                        % for a specific plane, clear data 
    %                     if a(planeIndex).odd_cprtime < a(planeIndex).even_cprtime
    %                         a(planeIndex).odd_cprlat = [];
    %                         a(planeIndex).odd_cprlon = [];
    %                         a(planeIndex).odd_cprtime = [];
    %                     else
    %                         a(planeIndex).even_cprlat = [];
    %                         a(planeIndex).even_cprlon = [];
    %                         a(planeIndex).even_cprtime = [];
    %                     end
                            res
                        end
                    end
                elseif (vec2dec(message(1:5)) == 19)
                    if (vec2dec(message(6:8)) == 1) || (vec2dec(message(6:8)) == 2)
                        b(planeIndex).ew_dir = message(14);
                        b(planeIndex).ew_velocity = vec2dec(message(15:15+10-1));
                        b(planeIndex).ns_dir = message(25);
                        b(planeIndex).ns_velocity = vec2dec(message(26:26+10-1));
                        b(planeIndex).vert_rate_source = message(36);
                        b(planeIndex).vert_rate_sign = message(37);
                        b(planeIndex).vert_rate = vec2dec(message(38:38+9-1));
                        b(planeIndex).cprtime = startIndex;
                        b(planeIndex).accuracy = vec2dec(message(48:48+3-1));
                        res = decode_ground_ref_velocity(b(planeIndex));
                        res.plane = planeList(planeIndex,:);
                        res
                    elseif (vec2dec(message(6:8)) == 3)
                        if message(14)
                            c(planeIndex).heading = vec2dec(message(15:15+10-1));
                            c(planeIndex).hori_velocity = vec2dec(message(26:26+10-1));
                            c(planeIndex).vert_rate = vec2dec(message(38:38+9-1));
                            c(planeIndex).vert_rate_sign = message(37);
                            c(planeIndex).cprtime = startIndex;
                            res2 = decode_air_ref_velocity(c(planeIndex));
                            res2.plane = planeList(planeIndex,:);
                            res2
                        end
                    end
                end
            end
        end
    end

    decodedDataCell = clear_decoded_data(decodedDataCell);
end

function [b] = decode_ground_ref_velocity(b)
    ewv = (b.ew_velocity-1)*1.852/3.6;  % convert unit from knot to m/s
    nsv = (b.ns_velocity-1)*1.852/3.6;  % convert unit from knot to m/s
    vv = (b.vert_rate-1)*64*0.00508;  % convert unit from feet/min to m/s

    if b.ew_dir   % 0 means east direction, 1 means west
        ewv = -1*ewv;
    end

    if b.ns_dir   % 0 means north direction, 1 means south direction 
        nsv = -1*nsv;
    end

    if b.vert_rate_sign  % 0 means up, 1 means down
        vv = -1*vv;
    end

    b.velocity = [ewv nsv vv]; % positive for East, North, Up direction, unit is m/s

    if (b.accuracy==0)
        b.accuracy = 'Error > 10m/s';
    elseif (b.accuracy==1)
        b.accuracy = 'Error < 10m/s';
    elseif (b.accuracy==2)
        b.accuracy = 'Error < 3m/s';
    elseif (b.accuracy==3)
        b.accuracy = 'Error < 1m/s';
    elseif (b.accuracy==4)
        b.accuracy = 'Error < 0.3m/s';
    end
end

function [c] = decode_air_ref_velocity(c)
    heading  = c.heading*360/1024; % value is in degree

    hv = (c.hori_velocity-1)*1.852/3.6;  % convert unit from knot to m/s
    vv = (c.vert_rate-1)*64*0.00508;  % convert unit from feet/min to m/s

    if c.vert_rate_sign  % 0 means up, 1 means down
        vv = -1*vv;
    end

    ewv = hv*sin(heading*pi/180);
    nsv = hv*cos(heading*pi/180);

    c.velocity = [ewv nsv vv]; % positive for East, North, Up direction, unit is m/s
end

function [dec] = vec2dec(data)
    dec = bin2dec(regexprep(num2str(data),'[^\w'']',''));
end