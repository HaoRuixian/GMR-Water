function Lcar = get_wave_length(gnss_system, band, sat)

load("glonasswlen.mat")

if gnss_system == "GPS"
    if band(2) == '1'
        Lcar = (299792458/1575.42e06); % for GPS L1       % 波长
    elseif band(2) == '2'
        Lcar = (299792458/1227.6e06);% L2
    elseif band(2) == '5'
        Lcar = (299792458/1176.45e06);% L5
    end
elseif gnss_system == "GLONASS"
    if sat>=24
        Lcar = 0;
        return
    end
    if band(2) == '1'
        Lcar = glonasswlen(sat); % for GLONASS L1       % 波长
    elseif band(2) == '2'
        Lcar = glonasswlen(sat);
        Lcar = 9/7 * Lcar;% L2
    else
        Lcar = 0;
    end
elseif gnss_system == "GALILEO"
    if band(2) == '1'
        Lcar = (299792458/1575.42e06); % E1 for galileo      % 波长
    elseif band(2) == '5'
        Lcar = (299792458/1176.45e06);% E5a
    elseif band(2) == '7'
        Lcar = (299792458/1207.14e06);% E5b
    elseif band(2) == '8'
        Lcar = (299792458/1191.795e06);% E5
    elseif band(2) == '6'
        Lcar = (299792458/1278.75e06);% E6
    end
elseif gnss_system == "BDS"
    if band(2) == '2'
        Lcar = (299792458/1561.098e06); % B1I for BDS      % 波长
    elseif band(2) == '7'
        Lcar = (299792458/1207.14e06);% B2I
    elseif band(2) == '6'
        Lcar = (299792458/1268.52e06);% B3I
    elseif band(2) == '1'
        Lcar = (299792458/1575.42e06);% B1C
    elseif band(2) == '5'
        Lcar = (299792458/1176.45e06);% B2a
    elseif band(2) == '9'
        Lcar = (299792458/1207.14e06);% B2b
    elseif band(2) == '8'
        Lcar = (299792458/1191.795e06);% B2a+b
    end

end