function [snr_data] = load_SNRdata(startdate, enddate, sp3option, station,d, all_days)
addpath('functions')
global station_l station_xyz elv azimask azi station_belong station_name obs_path sp3_path dt
tdatenum = startdate - 1;
datastr = strcat(obs_path,'\');

% LOAD SNR DATA
while tdatenum < enddate
    tdatenum = tdatenum + 1;
    curdt = datetime(tdatenum,'convertfrom','datenum');%转换为 datetime 数组
    curjd = juliandate(curdt);%将 MATLAB 日期时间转换为儒略日
    [gpsw,sow,~] = jd2gps(curjd);%gps周 ；seconds of week since 0 hr, Sun
    dow = sow/86400;%天

    strday=char(datetime(curdt,'format','DDD'));
    stryr=char(datetime(curdt,'format','yy'));
    stryrl=char(datetime(curdt,'format','yyyy'));

    if numel(dt) == 1
        dts = strcat('0',dt);
    else
        dts = dt;
    end
    if exist([datastr,upper(station),'00',station_belong,'_R_',stryrl,strday,'0000_01D_',num2str(dts),'S_MO.rnx'],'file') == 2
        obsstr=[datastr,upper(station),'00',station_belong,'_R_',stryrl,strday,'0000_01D_',num2str(dts),'S_MO.rnx'];
    elseif exist([datastr,'/',station,strday,'0.',stryr,'o'],'file')
        obsstr=[datastr,'/',station,strday,'0.',stryr,'o'];
    else
        disp('snr data does not exist!')
        continue
    end

    % the different precision ephemeris mechanism
    if sp3option == 'com'
        % OPTION 1
        sp3str = [sp3_path,'com',num2str(gpsw),num2str(round(dow)),'.sp3'];
    elseif sp3option == 'igs'
        % OPTION 2
        sp3str=[sp3_path,'\igs',num2str(gpsw),num2str(round(dow)),'.sp3'];
    elseif sp3option == 'COD'
        % OPTION 3
        sp3str=[sp3_path,'COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
    elseif sp3option == 'GFZ'
        % OPTION 4
        sp3str=[sp3_path,'GFZ0MGXRAP_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
    end


    elv_lims = elv;
    azi_lims = azi;
    azi_mask = azimask;
    plat = station_l(2);
    plon = station_l(1);
    staxyz = station_xyz;

    [snr_azi, snr_all ,snr_data] = rinex2snrfile_2(app,obsstr,sp3str,elv_lims,azi_lims,azi_mask,staxyz,plat,plon,d,all_days);
    save([app.path_SNR.Value, '\', station_name, num2str(tdatenum),'.mat'],'snr_data');
    save([app.path_SNR.Value, '\', station_name, num2str(tdatenum),'unselected.mat'],'snr_all');
    save([app.path_SNR.Value, '\', station_name, num2str(tdatenum),'azi.mat'],'snr_azi');


    clear slvlr lspy


end
end