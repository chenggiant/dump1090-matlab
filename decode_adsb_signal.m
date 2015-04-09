function [packetDataCell] = decode_adsb_signal(skipSeconds, filename, dataPath,step,sampRate)

    nbrSamples = step * sampRate;
    q = ceil(sampRate/1e7*10);
 
    fc = 1090e6;
    Ts = 1/sampRate;

    % Loading data
    fprintf(sprintf(' -- Loading data ... \n'));

    rxSignal = read_binary(sprintf('%s/%s',dataPath,filename),round(skipSeconds*sampRate),nbrSamples);

    rxSignal = resample(rxSignal, 10, q);
    sampRate = 10000000;
    nbrSamples = step * sampRate;

    % Low-pass filtering
    fprintf(sprintf(' -- Filtering data ... \n'));
    [b,a] = butter(5,0.25,'low'); % cutoff freq = 0.02*sampRate/2, filter order = 5
    rxSignal = filter(b,a,rxSignal);

    % Finding the start of packets
    fprintf(sprintf(' -- Correlating with packet headers ... \n'));

    preamble = [1 0 1 0 0 0 0 1 0 1];
    preamble = kron(preamble, ones(1,5));
    N = length(preamble);
    for cpt = 1:length(rxSignal)-length(preamble)
        corrFct(cpt) = (preamble) * abs(rxSignal(cpt:cpt+N-1));
    end

    treshold_corrFct = 200;
    treshold_totAmpl = 2500;
    r = abs(rxSignal);
    cpt2 = 0;
    for cpt = 1:length(corrFct) - 560-80-20 %20  !!!! be careful about this part
        % Verify if correlation function has peak
        if corrFct(cpt) > treshold_corrFct
            % Checking if header shape is detected
            if ( (r(cpt)>r(cpt+1*5)) && ...
                    (r(cpt+1*5)<r(cpt+2*5)) && ...
                    (r(cpt+2*5)>r(cpt+3*5)) && ...
                    (r(cpt+3*5)<r(cpt+0*5)) && ...
                    (r(cpt+4*5)<r(cpt+0*5)) && ...
                    (r(cpt+5*5)<r(cpt+0*5)) && ...
                    (r(cpt+6*5)<r(cpt+0*5)) && ...
                    (r(cpt+7*5)>r(cpt+8*5)) && ...
                    (r(cpt+8*5)<r(cpt+9*5)) && ...
                    (r(cpt+9*5)>r(cpt+6*5)) )
                % Compute power of peaks
                high = (r(cpt)+r(cpt+2*5)+r(cpt+7*5)+r(cpt+9*5))/4;
                % Verify if low between peaks is lower than average of peaks
                if ( (r(cpt+4*5)<high) && ...
                        (r(cpt+5*5)<high) )
                    % Verify if low after peaks is lower than average of peaks
                    if ( (r(cpt+11*5)<high) && ...
                            (r(cpt+12*5)<high) && ...
                            (r(cpt+13*5)<high) && ...
                            (r(cpt+14*5)<high) )
                        % Verify if the packet is long enough to be a data
                        % packet
                        L = 80+560;
                        amplMaxPacket = max(abs(r(cpt:cpt+L)));
                        amplMinPacket = min(abs(r(cpt:cpt+L)));
                        meanAmpl = (amplMaxPacket-amplMinPacket)/2;
                        I = find(abs(r(cpt:cpt+L))>meanAmpl);
                        if (length(I) > 0.5*L/2)
                            cpt2 = cpt2+1;
                            headerIndex(cpt2) = cpt;
                        end
                    end
                end
            end
        end
    end

    % Cleaning up headerIndex vector: when multiple indices are following each
    % other, we select the one with highest amplitude
    I = diff(headerIndex);
    cpt  = 1;
    cpt2 = 0;
    while(cpt<length(headerIndex)-4)
        if I(cpt)~=1
            % next sample far away
            cpt2 = cpt2+1;
            headerIndex_final(cpt2) = headerIndex(cpt);
            cpt = cpt+1;
        elseif ((I(cpt)==1) && (I(cpt+1)~=1))
            cpt2 = cpt2+1;
            [~,ind] = max(abs(r(headerIndex(cpt):headerIndex(cpt)+1)));
            headerIndex_final(cpt2) = headerIndex(cpt+ind-1);
            cpt = cpt+2;
        elseif ((I(cpt)==1) && (I(cpt+1)==1) && (I(cpt+2)~=1))
            cpt2 = cpt2+1;
            [~,ind] = max(abs(r(headerIndex(cpt):headerIndex(cpt)+2)));
            headerIndex_final(cpt2) = headerIndex(cpt+ind-1);
            cpt = cpt+3;
        elseif ((I(cpt)==1) && (I(cpt+1)==1) && (I(cpt+2)==1) && (I(cpt+3)~=1))
            cpt2 = cpt2+1;
            [~,ind] = max(abs(r(headerIndex(cpt):headerIndex(cpt)+3)));
            headerIndex_final(cpt2) = headerIndex(cpt+ind-1);
            cpt = cpt+4;
        elseif ((I(cpt)==1) && (I(cpt+1)==1) && (I(cpt+2)==1) && (I(cpt+3)==1) && (I(cpt+4)~=1))
            cpt2 = cpt2+1;
            [~,ind] = max(abs(r(headerIndex(cpt):headerIndex(cpt)+4)));
            headerIndex_final(cpt2) = headerIndex(cpt+ind-1);
            cpt = cpt+5;
        elseif ((I(cpt)==1) && (I(cpt+1)==1) && (I(cpt+2)==1) && (I(cpt+3)==1) && (I(cpt+4)==1) && (I(cpt+5)~=1))
            cpt2 = cpt2+1;
            [~,ind] = max(abs(r(headerIndex(cpt):headerIndex(cpt)+5)));
            headerIndex_final(cpt2) = headerIndex(cpt+ind-1);
            cpt = cpt+6;
        end
    end

    % Clean up multiple index for long packet
    j = 1;
    headerIndex_remove = [];
    for i = 1:length(headerIndex_final)-1
        if headerIndex_final(i+1)-headerIndex_final(i) < 56*10
            headerIndex_remove(j) = headerIndex_final(i+1);
            j = j+1;
        end
    end

    headerIndex_final = setdiff(headerIndex_final, headerIndex_remove);

    % Decode the packet
    fprintf(sprintf(' -- Decoding packets ... \n \n'));

    planeIndex = 1;
    packetIndex = 1;

    for indexStartPacket = headerIndex_final(1:end-1)   % the last packet may cause index exceeds matrix dimension

        L = (8+112)*10;

        % Undersample signal
        rxSymbols = abs(rxSignal(indexStartPacket:5:indexStartPacket+L-1));
        header  = rxSymbols(1:16);
        payload = rxSymbols(17:end);

        cptBits = 0;
        for cptSymb = 1:2:length(payload)
            cptBits = cptBits+1;
            if (payload(cptSymb) > payload(cptSymb+1))
                bits(cptBits) = 1;
            else
                bits(cptBits) = 0;
            end
        end

        if ( isequal(bits(1:5), [1 0 0 0 0]) || ...
                isequal(bits(1:5), [1 0 0 0 1]) || ...
                isequal(bits(1:5), [1 0 0 1 1]) || ...
                isequal(bits(1:5), [1 0 1 0 0]) || ...
                isequal(bits(1:5), [1 0 1 0 1]) )
            M = 112;
        else
            M = 56;
        end

        checksum_table_hex =  [...
            '3935ea' '1c9af5' 'f1b77e' '78dbbf' 'c397db' '9e31e9' 'b0e2f0' '587178' ...
            '2c38bc' '161c5e' '0b0e2f' 'fa7d13' '82c48d' 'be9842' '5f4c21' 'd05c14' ...
            '682e0a' '341705' 'e5f186' '72f8c3' 'c68665' '9cb936' '4e5c9b' 'd8d449' ...
            '939020' '49c810' '24e408' '127204' '093902' '049c81' 'fdb444' '7eda22' ...
            '3f6d11' 'e04c8c' '702646' '381323' 'e3f395' '8e03ce' '4701e7' 'dc7af7' ...
            '91c77f' 'b719bb' 'a476d9' 'adc168' '56e0b4' '2b705a' '15b82d' 'f52612' ...
            '7a9309' 'c2b380' '6159c0' '30ace0' '185670' '0c2b38' '06159c' '030ace' ...
            '018567' 'ff38b7' '80665f' 'bfc92b' 'a01e91' 'aff54c' '57faa6' '2bfd53' ...
            'ea04ad' '8af852' '457c29' 'dd4410' '6ea208' '375104' '1ba882' '0dd441' ...
            'f91024' '7c8812' '3e4409' 'e0d800' '706c00' '383600' '1c1b00' '0e0d80' ...
            '0706c0' '038360' '01c1b0' '00e0d8' '00706c' '003836' '001c1b' 'fff409' ...
            '000000' '000000' '000000' '000000' '000000' '000000' '000000' '000000' ...
            '000000' '000000' '000000' '000000' '000000' '000000' '000000' '000000' ...
            '000000' '000000' '000000' '000000' '000000' '000000' '000000' '000000' ...
            ];
        for cpt = 1:112
            checksum_table_dec(cpt,1) = hex2dec(checksum_table_hex((cpt-1)*6+1:cpt*6));
        end
        for cpt = 1:112
            checksum_table_bin(cpt,:) = dec2bin(checksum_table_dec(cpt,1),6*4);
        end
        % Reduce checksum table when dealing with short packets of 56 bits
        if M == 112
            checksum_table_bin_reduced = checksum_table_bin;
        elseif M == 56
            checksum_table_bin_reduced = checksum_table_bin(57:end,:);
            bits(57:end) = [];
        end

        % Do CRC computation, only analyze DF 11 and DF 17, since their CRC are
        % not XORed with Address

        if (isequal(bits(1:5), [0 1 0 1 1]) || ...
                isequal(bits(1:5), [1 0 0 0 1]) )

            crc_temp = '000000000000000000000000';

            for cptBits = 1:M
                if bits(cptBits)==1
                    crc_temp = xor_crc(crc_temp,checksum_table_bin_reduced(cptBits,:));
                end
            end
            crc = crc_temp;

            % Verify if CRC is correct
            crc_ok = verify_crc(bits(end-23:end),crc);
            if crc_ok
                packet = decode_packet(bits, sampRate,indexStartPacket,skipSeconds);
                packetDataCell{packetIndex, :} = packet;
                packetIndex = packetIndex + 1;

            else    % if crc check fails, flip one bit and run crc check again.
                for flipBit = 1:M
                    bits_temp = bits;

                    if bits(flipBit)==1
                        bits_temp(flipBit)=0;
                    elseif bits(flipBit)==0
                        bits_temp(flipBit)=1;
                    end

                    crc_temp = '000000000000000000000000';

                    for cptBits = 1:M
                        if bits_temp(cptBits)==1
                            crc_temp = xor_crc(crc_temp,checksum_table_bin_reduced(cptBits,:));
                        end
                    end
                    crc = crc_temp;

                    crc_ok = verify_crc(bits_temp(end-23:end),crc);

                    if crc_ok
                        packet = decode_packet(bits_temp, sampRate,indexStartPacket,skipSeconds);
                        packetDataCell{packetIndex, :} = packet;
                        packetIndex = packetIndex + 1;
                    end
                end
            end
        end
    end

    if ~exist('packetDataCell')
        packetDataCell = [];
    end
end


% Decode plane ID and packet data type
function [packet] = decode_packet(bits, sampRate,indexStartPacket,skipSeconds)

    doPrint = 1;

    planeId = dec2hex(bin2dec(regexprep(num2str(bits(9:9+23)),'[^\w'']','')));
    padId = regexprep(num2str(zeros(1, 6-length(planeId))),'[^\w'']','');
    planeId = strcat(padId, planeId);

    if doPrint
        fprintf(sprintf(' -- We find a packet ... \n'));
        fprintf(sprintf(' -- Start index is %d ... \n', ...
            indexStartPacket+round(skipSeconds*sampRate)));
        fprintf(sprintf(' -- DF type is %d ... \n', vec2dec(bits(1:5))));
        fprintf(sprintf(' -- Plane ID is %s ... \n \n', planeId));
    end

    % save packet to packetDataCell
    packet.planeId = planeId;
    packet.bits = bits;
    packet.startIndex = indexStartPacket+round(skipSeconds*sampRate);
    packet.DF = vec2dec(bits(1:5));
end

function crc_ok = verify_crc(bits,crc)
    crc_ok = 1;
    for cpt = 1:length(bits)
        if crc_ok
            if ( (bits(cpt)==1)&&(crc(cpt)=='1') )
                crc_ok = 1;
            elseif ( (bits(cpt)==0)&&(crc(cpt)=='0') )
                crc_ok = 1;
            else
                crc_ok = 0;
            end
        end
    end
end

function result = xor_crc(sig1,sig2)
    for cpt = 1:length(sig1)
        if sig1(cpt) == sig2(cpt)
            result(cpt) = '0';
        else
            result(cpt) = '1';
        end
    end
end

function [dec] = vec2dec(data)
    dec = bin2dec(regexprep(num2str(data),'[^\w'']',''));
end

% Read binary data, and output baseband complex signal

% jumper_length: sets the file position indicator offset bytes from 
% origin in the specified file

function [baseband_sig] = read_binary(file_name, jumper_length, sample_length)
    fid = fopen(file_name,'r');
    fseek(fid,jumper_length*4,'bof');
    IQdata = fread(fid, 2*sample_length, 'short');
    fclose(fid);
    baseband_sig = IQdata(1:2:2*sample_length) + 1i*IQdata(2:2:2*sample_length);
end