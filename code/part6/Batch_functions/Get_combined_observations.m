function [M,a,b,cur_band] = Get_combined_observations(Meth_id, cur_sys, inverse_data, indx, cur_sat)


if Meth_id == 3                              % Carrier phase dual/three frequency combination
    if string(cur_sys) == "GPS"
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        L2 = inverse_data{indx,get_col_name(inverse_data,indx,'L2')};
        L3 = inverse_data{indx,get_col_name(inverse_data,indx,'L5')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L2',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L5',cur_sat);
    elseif string(cur_sys) == "GALILEO"
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        L2 = inverse_data{indx,get_col_name(inverse_data,indx,'L8')};
        L3 = inverse_data{indx,get_col_name(inverse_data,indx,'L7')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L8',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L7',cur_sat);
    elseif string(cur_sys) == "BDS"
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        L2 = inverse_data{indx,get_col_name(inverse_data,indx,'L6')};
        L3 = inverse_data{indx,get_col_name(inverse_data,indx,'L5')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L6',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L5',cur_sat);
    elseif string(cur_sys) == "GLONASS"
        lamda = get_wave_length(cur_sys,'L7',cur_sat);
    end
    yita1 = lamda3^2 - lamda2^2;
    yita2 = lamda1^2 - lamda3^2;
    yita3 = lamda2^2 - lamda1^2;

    kama1 = yita1 * lamda1;
    kama2 = yita2 * lamda2;
    kama3 = yita3 * lamda3;

    M = kama1*L1 + kama2*L2 + kama3*L3;
    cur_band = 'Carrier';
    if string(cur_sys)     == "GPS"
        a = 0.1248;
        b = -0.024;
    elseif string(cur_sys) == "GALILEO"
        a = 0.1257;
        b = -0.0548;
    elseif string(cur_sys) == "BDS"
        a = 0.1207;
        b = -0.25;
    end
    return

elseif Meth_id == 4                         % pseudorange dual/three frequency combination
    if string(cur_sys) == "GPS"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C1')};
        C2 = inverse_data{indx,get_col_name(inverse_data,indx,'C2')};
        C3 = inverse_data{indx,get_col_name(inverse_data,indx,'C5')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L2',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L5',cur_sat);
    elseif string(cur_sys) == "GALILEO"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C1')};
        C2 = inverse_data{indx,get_col_name(inverse_data,indx,'C8')};
        C3 = inverse_data{indx,get_col_name(inverse_data,indx,'C7')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L8',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L7',cur_sat);
    elseif string(cur_sys) == "BDS"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C2')};
        C2 = inverse_data{indx,get_col_name(inverse_data,indx,'C7')};
        C3 = inverse_data{indx,get_col_name(inverse_data,indx,'C6')};
        lamda1 = get_wave_length(cur_sys,'L2',cur_sat);
        lamda2 = get_wave_length(cur_sys,'L7',cur_sat);
        lamda3 = get_wave_length(cur_sys,'L6',cur_sat);
    end
    yita1 = lamda3^2 - lamda2^2;
    yita2 = lamda1^2 - lamda3^2;
    yita3 = lamda2^2 - lamda1^2;

    M = yita1*C1 + yita2*C2 + yita3*C3;
    cur_band = 'Pseudorange';
    if string(cur_sys) == "GPS"
        a = 0.1248;
        b = -0.024;
    elseif string(cur_sys) == "GALILEO"
        a = 0.1257;
        b = -0.0548;
    elseif string(cur_sys) == "BDS"
        a = 0.1207;
        b = -0.0013;
    end
    return

elseif Meth_id == 5                         % Single frequency carrier phase pseudo range combination
    if string(cur_sys) == "GPS"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C1')};
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
    elseif string(cur_sys) == "GALILEO"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C1')};
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
    elseif string(cur_sys) == "BDS"
        C1 = inverse_data{indx,get_col_name(inverse_data,indx,'C1')};
        L1 = inverse_data{indx,get_col_name(inverse_data,indx,'L1')};
        lamda1 = get_wave_length(cur_sys,'L1',cur_sat);
    end
    M = C1 - L1*lamda1;
    cur_band = 'CP';
    if string(cur_sys) == "GPS"
        a = 0.0951;
        b = 0.0016;
    elseif string(cur_sys) == "GALILEO"
        a = 0.0951;
        b = 0.0016;
    elseif string(cur_sys) == "BDS"
        a = 0.0951;
        b = 0.0016;
    end
    return

end
end

function name = get_col_name(inverse_data,indx,point)
name = inverse_data.Properties.VariableNames(startsWith(inverse_data.Properties.VariableNames,point));
nan_counts = sum(ismissing(inverse_data(indx, name)), 1);
[~, idx] = min(nan_counts);
name = name{idx};
end