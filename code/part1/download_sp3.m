function download_sp3(gpsweek,sow,dow,sp3_option,save_path,d,day_num,numday)
%This script is to download eph files from igs.ign.fr
%
%--------------------------------------------------------------------------
%         By: Ruixian Hao
%           Contact: Vitamin_N@outlook.com
%           China University of Mining and Technology-Beijing
%
%           Dec 2024
%--------------------------------------------------------------------------

jd = gps2jd(gpsweek,sow,0);
time_gps = gpsweek*10 + dow;
gpsweek = num2str(gpsweek);
curdt=datetime(jd,'convertfrom','juliandate');
strday=char(datetime(curdt,'format','DDD'));
stryrl=char(datetime(curdt,'format','yyyy'));

if sp3_option == 'COD'
    ftpobj = ftp('igs.ign.fr');
    cd(ftpobj,strcat('pub/igs/products/mgex/',gpsweek));

    filename = strcat('COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3.gz');
    mget(ftpobj,filename,save_path);
    gunzip(strcat(save_path,'\',filename),strcat(save_path,'\',sp3_option));
    delete(strcat(save_path,'\',filename))

    % progress bar
    d.Value = day_num / numday;
    progress = strcat('(',num2str(day_num),'/',num2str(numday),')'); % update progress
    d.Title = strcat('Download progress ',progress);
    d.Message = strcat('File：','COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3','Saved to: ',save_path,'\',sp3_option);

elseif sp3_option == 'GFZ'
    ftpobj = ftp('igs.ign.fr');
    cd(ftpobj,strcat('pub/igs/products/mgex/',gpsweek));

    filename = strcat('GFZ0MGXRAP_',stryrl,strday,'0000_01D_05M_ORB.SP3.gz');
    mget(ftpobj, filename, save_path);
    gunzip(strcat(save_path,'\',filename),strcat(save_path,'\',sp3_option));
    delete(strcat(save_path,'\',filename))

    % progress bar
    d.Value = day_num / numday;
    progress = strcat('(',num2str(day_num),'/',num2str(numday),')');
    d.Title = strcat('Download progress ',progress);
    d.Message = strcat('File：','GFZ0MGXRAP_',stryrl,strday,'0000_01D_05M_ORB.SP3','Saved to: ',save_path,'\',sp3_option);

elseif sp3_option == 'com'
    filename = strcat('com',num2str(time_gps),'.sp3.Z');
    mget(ftpobj,filename,save_path);
    
    % progress bar
    d.Value = day_num / numday;
    progress = strcat('(',num2str(day_num),'/',num2str(numday),')');
    d.Title = strcat('Download progress ',progress);
    d.Message = strcat('File：','com',num2str(time_gps),'.sp3','Saved to：',save_path,'\',sp3_option);
elseif sp3_option == 'WUM'
    ftpobj = ftp('igs.ign.fr');
    cd(ftpobj,strcat('pub/igs/products/mgex/',gpsweek));

    filename = strcat('WUM0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3.gz');
    mget(ftpobj, filename, save_path);
    gunzip(strcat(save_path,'\',filename),strcat(save_path,'\',sp3_option));
    delete(strcat(save_path,'\',filename))

    % progress bar
    d.Value = day_num / numday;
    progress = strcat('(',num2str(day_num),'/',num2str(numday),')');
    d.Title = strcat('Download progress ',progress);
    d.Message = strcat('File：','WUM0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3','Saved to: ',save_path,'\',sp3_option);
end
end