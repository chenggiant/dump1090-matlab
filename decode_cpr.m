function [a] = decode_cpr(a)

    AirDlat0 = 360.0 / 60;
    AirDlat1 = 360.0 / 59;
    lat0 = a.even_cprlat;
    lat1 = a.odd_cprlat;
    lon0 = a.even_cprlon;
    lon1 = a.odd_cprlon;

    % Compute the Latitude Index j
    j = floor(((59*lat0 - 60*lat1) / 131072) + 0.5);
    rlat0 = AirDlat0 * (cprModFunction(j,60) + lat0 / 131072);
    rlat1 = AirDlat1 * (cprModFunction(j,59) + lat1 / 131072);

    if (rlat0 >= 270)
        rlat0 = rlat0 - 360;
    end

    if (rlat1 >= 270)
        rlat1 = rlat1 - 360;
    end

    % Check that both are in the same latitude zone, or abort.
    if (cprNLFunction(rlat0) ~= cprNLFunction(rlat1))
        return;
    end

    % Compute ni and the longitude index m
    if (a.even_cprtime > a.odd_cprtime)
        % Use even packet
        ni = cprNFunction(rlat0,0);
        m = floor((((lon0 * (cprNLFunction(rlat0)-1))-(lon1 * cprNLFunction(rlat0))) / 131072) + 0.5);
        a.lat = rlat0;
        a.lon = cprDlonFunction(rlat0,0) * (cprModFunction(m,ni)+lon0/131072);
    else
        % Use odd packet
        ni = cprNFunction(rlat1,1);
        m = floor((((lon0 * (cprNLFunction(rlat1)-1))-(lon1 * cprNLFunction(rlat1))) / 131072) + 0.5);
        a.lat = rlat1;
        a.lon = cprDlonFunction(rlat1,1) * (cprModFunction(m,ni)+lon1/131072);
    end

    if (a.lon > 180)
        a.lon = a.lon - 360;
    end
end

function [n] = cprNLFunction(lat)
    if (lat < 0)
        lat = -lat;  % Table is simmetric about the equator
    end

    if (lat < 10.47047130)
        n = 59;
    elseif (lat < 14.82817437)
        n = 58;
    elseif (lat < 18.18626357)
        n = 57;
    elseif (lat < 21.02939493)
        n = 56;
    elseif (lat < 23.54504487)
        n = 55;
    elseif (lat < 25.82924707)
        n = 54;
    elseif (lat < 27.93898710)
        n = 53;
    elseif (lat < 29.91135686)
        n = 52;
    elseif (lat < 31.77209708)
        n = 51;
    elseif (lat < 33.53993436)
        n = 50;
    elseif (lat < 35.22899598)
        n = 49;
    elseif (lat < 36.85025108)
        n = 48;
    elseif (lat < 38.41241892)
        n = 47;
    elseif (lat < 39.92256684)
        n = 46;
    elseif (lat < 41.38651832)
        n = 45;
    elseif (lat < 42.80914012)
        n = 44;
    elseif (lat < 44.19454951)
        n = 43;
    elseif (lat < 45.54626723)
        n = 42;
    elseif (lat < 46.86733252)
        n = 41;
    elseif (lat < 48.16039128)
        n = 40;
    elseif (lat < 49.42776439)
        n = 39;
    elseif (lat < 50.67150166)
        n = 38;
    elseif (lat < 51.89342469)
        n = 37;
    elseif (lat < 53.09516153)
        n = 36;
    elseif (lat < 54.27817472)
        n = 35;
    elseif (lat < 55.44378444)
        n = 34;
    elseif (lat < 56.59318756)
        n = 33;
    elseif (lat < 57.72747354)
        n = 32;
    elseif (lat < 58.84763776)
        n = 31;
    elseif (lat < 59.95459277)
        n = 30;
    elseif (lat < 61.04917774)
        n = 29;
    elseif (lat < 62.13216659)
        n = 28;
    elseif (lat < 63.20427479)
        n = 27;
    elseif (lat < 64.26616523)
        n = 26;
    elseif (lat < 65.31845310)
        n = 25;
    elseif (lat < 66.36171008)
        n = 24;
    elseif (lat < 67.39646774)
        n = 23;
    elseif (lat < 68.42322022)
        n = 22;
    elseif (lat < 69.44242631)
        n = 21;
    elseif (lat < 70.45451075)
        n = 20;
    elseif (lat < 71.45986473)
        n = 19;
    elseif (lat < 72.45884545)
        n = 18;
    elseif (lat < 73.45177442)
        n = 17;
    elseif (lat < 74.43893416)
        n = 16;
    elseif (lat < 75.42056257)
        n = 15;
    elseif (lat < 76.39684391)
        n = 14;
    elseif (lat < 77.36789461)
        n = 13;
    elseif (lat < 78.33374083)
        n = 12;
    elseif (lat < 79.29428225)
        n = 11;
    elseif (lat < 80.24923213)
        n = 10;
    elseif (lat < 81.19801349)
        n = 9;
    elseif (lat < 82.13956981)
        n = 8;
    elseif (lat < 83.07199445)
        n = 7;
    elseif (lat < 83.99173563)
        n = 6;
    elseif (lat < 84.89166191)
        n = 5;
    elseif (lat < 85.75541621)
        n = 4;
    elseif (lat < 86.53536998)
        n = 3;
    elseif (lat < 87.00000000)
        n = 2;
    else
        n = 1;
    end
end

function [value]= cprDlonFunction(lat, isodd)
    value =  360.0 / cprNFunction(lat, isodd);
end

function [res] = cprModFunction(a, b)
    res = mod(a,b);
    
    if (res < 0)
        res = res + b;
    end
end

function [nl] = cprNFunction(lat, isodd)
    nl = cprNLFunction(lat) - isodd;

    if (nl < 1)
        nl = 1;
    end
end

