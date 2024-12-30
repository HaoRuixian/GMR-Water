function [GPS_azi, GPS_elv] = compute_azi_elv(secs, GPS_data_time, GPS_data_prn_num, xyzt, staxyz,plat,plon,d, error_ind,pace,l)
num = numel(GPS_data_time);
GPS_azi = nan(num,1);
GPS_elv = nan(num,1);
ind = 1;
bool = 0;
[all_ok,~] = size(GPS_data_time);
for j = 1:numel(secs)
    time_indx = GPS_data_time==secs(j);
    satnum = GPS_data_prn_num(time_indx); % the prn number in this time
    for satindx = 1:numel(satnum)
        satindx_num = satnum(satindx);
        if error_ind(ind) ~= 0
            ind = ind+1;
        end

        if (secs(j) == GPS_data_time(ind) && satindx_num == GPS_data_prn_num(ind))
            [azi,elv] = gnss2azelv(staxyz,xyzt(j,:,satindx_num),plat,plon);
            GPS_azi(ind,1) = azi;
            GPS_elv(ind,1) = elv;
            ind = ind + 1;
            d.Value = (ind/all_ok)*(1/4)*pace + l;
            if ind > numel(GPS_data_time)
                bool = 1;
                break
            end
        end

    end
    if bool == 1
        return
    end
end
end