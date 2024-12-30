%--------------------------------------------------------------------------
% Inverse the RH based on multiple methods
%         By: Ruixian Hao
%           Contact: Vitamin_N@outlook.com
%           China University of Mining and Technology-Beijing
%
%           Dec 2024
%--------------------------------------------------------------------------

%% Data path
snr_savepath         = [settings.Out_path,'/SNR_file'];
carrier_savepath     = [settings.Out_path,'/Carrier_file'];
pseudorange_savepath = [settings.Out_path,'/Pseudorange_file'];

%% Inversion is carried out by five methods
if string(settings.par) == "None"
    for tdatenum = start_date: end_date
        inverse_code
    end
else
    par = parpool(str2double(settings.par));
    parfor tdatenum = start_date: end_date
        varNames = {'Time','System','BAND','PRN','ROC','MIN_elv','MAX_elv','MEAN_AZI','RH','trop_c'};
        varTypes = {'datetime','string','string','double','double','double','double','double','double','double'};
        RH_info = table('Size',[0,length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
        rid   = 0;

        for Meth_id = 1:5
            % presuppose
            curdt = datetime(tdatenum,'convertfrom','datenum'); 
            curjd = juliandate(curdt);
            curyr = datetime(curdt,'format','yyyy');
            curyr = double(int16(convertTo(curyr,'yyyymmdd') / 10000));
            curday = datetime(curdt,'format','DDD');
            curday = day(curday,'dayofyear');

            % Trop pre settings
            lla       = ecef2lla(staxyz);
            hell      = lla(3);
            gpt3_grid = gpt3_5_fast_readGrid;
            curjd2    = doy2jd(curyr,curday)-2400000.5;
            hgtlim    = nan(2,1);
            hgtlim(1) = hell - sta_asl+tide_range(1);
            hgtlim(2) = hell - sta_asl+tide_range(2);
            dlat      = sta_lat*pi/180.d0;
            dlon      = sta_lon*pi/180.d0;
            it        = 0;
            plim = nan(2,1);
            tlim = nan(2,1);
            tmlim = nan(2,1);
            elim = nan(2,1);

            [pant,~,~,tmant,eant,ah,aw,lambda,~] = gpt3_5_fast (curjd2, dlat,dlon, hell,it, gpt3_grid);
            [plim(1),tlim(1),~,tmlim(1),elim(1),~,~,~,~] = gpt3_5_fast (curjd2, dlat,dlon, hgtlim(1),it, gpt3_grid);
            [plim(2),tlim(2),~,tmlim(2),elim(2),~,~,~,~] = gpt3_5_fast (curjd2, dlat,dlon, hgtlim(2),it, gpt3_grid);

            % load Data
            if settings.methods(Meth_id)
                if Meth_id < 3      % inverse base on SNR
                    Data     = load([snr_savepath, '\', settings.station_name, num2str(tdatenum),'.mat']);
                    All_data = Data.snr_data;
                elseif Meth_id == 3 % base on carrier
                    Data     = load([carrier_savepath, '\', settings.station_name, num2str(tdatenum),'.mat']);
                    All_data = Data.carrier_data;
                elseif Meth_id == 4 % base on psedorange
                    Data     = load([pseudorange_savepath, '\', settings.station_name, num2str(tdatenum),'.mat']);
                    All_data = Data.pseudorange_data;
                elseif Meth_id == 5
                    Data      = load([carrier_savepath, '\', settings.station_name, num2str(tdatenum),'.mat']);
                    All_data1 = Data.carrier_data;
                    Data      = load([pseudorange_savepath, '\', settings.station_name, num2str(tdatenum),'.mat']);
                    All_data2 = Data.pseudorange_data;
                    All_data  = mergeStructTables(All_data1, All_data2);
                end
            else
                continue
            end

            %% Start inverse using All_data
            system     = fieldnames(All_data);
            system_num = numel(system);

            for s = 1:system_num                                    % GNSS
                cur_sys = system{s};
                cur_data= All_data.(cur_sys);
                prn     = unique(cur_data{:,1});
                sat_num = numel(prn);
                for sat = 1:sat_num                                 % Satellite
                    cur_sat = prn(sat);
                    inverse_indx = cur_data{:,1}==cur_sat;
                    inverse_data = cur_data(inverse_indx,:);
                    time_a   = inverse_data{:,2};
                    headers  = cur_data.Properties.VariableNames;
                    bands    = headers(5:end);
                    band_num = numel(bands);

                    %% For SNR invers
                    % spectral inverse
                    if Meth_id == 1
                        for band_id = 1:band_num          % band
                            cur_band = bands{band_id};
                            if cur_band(1) ~= 'S'   % only for snr
                                continue
                            end
                            time_gap = diff(time_a);
                            gap_indx = find((time_gap>10*settings.Rinex_dt) == 1);
                            for g = 1:numel(gap_indx)+1
                                if numel(gap_indx) == 0
                                    indx = 1:size(inverse_data,1);
                                else
                                    if g == 1
                                        indx = 1:gap_indx(g);
                                    elseif g == numel(gap_indx)+1
                                        indx = (gap_indx(g-1)+1):size(inverse_data,1);
                                    else
                                        indx = (gap_indx(g-1)+1) : gap_indx(g);
                                    end
                                end
                                time_o= inverse_data{indx,2};
                                azi_o = inverse_data{indx,3};
                                snr_o = inverse_data{indx,cur_band};
                                elv_o = inverse_data{indx,4};

                                elv_one = all(diff(elv_o)<0) || all(diff(elv_o)>0);
                                if elv_one
                                    arc_num = 1;
                                else
                                    arc_num = 2;
                                    try
                                        gap_elv = find(elv_o == findpeaks(elv_o));
                                    catch
                                        elv_o = -elv_o;
                                        gap_elv = find(elv_o == findpeaks(elv_o));
                                    end
                                end
                                for arc_id = 1:arc_num
                                    if arc_num ~= 1
                                        if arc_id == 1
                                            time = time_o(1:gap_elv);
                                            azi = azi_o(1:gap_elv);
                                            snr = snr_o(1:gap_elv);
                                            elv = elv_o(1:gap_elv);
                                        else
                                            time = time_o(gap_elv+1,:);
                                            azi = azi_o(gap_elv+1,:);
                                            snr = snr_o(gap_elv+1,:);
                                            elv = elv_o(gap_elv+1,:);
                                        end
                                    else
                                        time = time_o;
                                        azi = azi_o;
                                        snr = snr_o;
                                        elv = elv_o;
                                    end

                                    % get wave length
                                    wave_length = get_wave_length(cur_sys,cur_band,cur_sat);
                                    if wave_length == 0
                                        continue
                                    end

                                    % check nan
                                    nan_indx = isnan(snr);
                                    snr(nan_indx) = [];
                                    elv(nan_indx) = [];
                                    azi(nan_indx) = [];
                                    time(nan_indx) = [];
                                    sinelv = sind(elv);

                                    % inverse
                                    if numel(sinelv) < 5
                                        continue
                                    end

                                    [refl_h, id, psd, pks] = snr2RH_lsp(sinelv, snr, wave_length, hell, hgtlim);
                                    if isnan(pks)
                                        continue
                                    end

                                    % Tropospheric correction
                                    pre_h = refl_h(id);
                                    hsfc = hell - pre_h;
                                    if hsfc > hgtlim(1) && hsfc < hgtlim(2)
                                        psfc = interp1(hgtlim,plim,hsfc,'linear');
                                        tsfc = interp1(hgtlim,tlim,hsfc,'linear');
                                        tmsfc = interp1(hgtlim,tmlim,hsfc,'linear');
                                        esfc = interp1(hgtlim,elim,hsfc,'linear');

                                        dele = (1/60) * 510 * psfc / ((9/5*tsfc+492) * 1010.16) .*cotd(elv+7.31./(elv+4.4));
                                        sinelv = sind(elv + dele);
                                        elv = elv+dele;

                                        theta = elv;
                                        curdt = datetime(tdatenum,'convertfrom','datenum');
                                        curjd = juliandate(curdt);
                                        thetarefr = [];
                                        for jj = 1:numel(theta)
                                            tau = trop_delay_tp(curjd,sta_lat,hell,hsfc,theta(jj),pant,tmant,eant,ah,aw,lambda,psfc,tmsfc,esfc);
                                            thetarefr(jj) = asind( sind(theta(jj)) + 0.5*tau/pre_h );
                                        end
                                        sinelv = sind(thetarefr).';

                                        [refl_h, id, psd, pks] = snr2RH_lsp(sinelv, snr, wave_length, hell, hgtlim);
                                        if isnan(pks)
                                            continue
                                        end

                                    end
                                    cur_rh = refl_h(id);  % Final rh
                                    trop_c = cur_rh - pre_h; % Tropospheric correction
                                    valid = 1;
                                    if cur_rh > sta_asl-tide_range(1) || cur_rh < sta_asl-tide_range(2) || max(psd)<5*mean(pks(1:end-1))
                                        valid = 0;
                                    end

                                    if valid == 1
                                        rid = rid + 1; % row id
                                        Time = datetime(tdatenum + mean(time)/86400, 'ConvertFrom', 'datenum');       % datenum
                                        System = string(cur_sys);
                                        BAND = string(cur_band);
                                        PRN = cur_sat;                               % Sat prn
                                        ROC = tand(mean(elv)) / (( (pi/180) * (elv(end)-elv(1)) )...
                                            /(time(end)-time(1)));                   % tan(th)/dth/dt :rate of change
                                        MIN_elv = min(elv);                          % THETA MIN
                                        MAX_elv = max(elv);                          % THETA MAX
                                        MEAN_AZI = nanmean(azi);                     % MEAN AZI
                                        RH = cur_rh;                                 % slvl
                                        RH_info(rid,:) = table(Time,System,BAND,PRN,ROC,MIN_elv,MAX_elv,MEAN_AZI,RH,trop_c);
                                    end

                                end
                            end
                        end
                    end

                    % inverse modeling
                    if Meth_id == 2

                    end

                    %% For carrier and pseudorange inverse
                    if Meth_id > 2
                        time_gap = diff(time_a);
                        gap_indx = find((time_gap>10*settings.Rinex_dt) == 1);
                        for g = 1:numel(gap_indx)+1
                            if numel(gap_indx) == 0
                                indx = 1:size(inverse_data,1);
                            else
                                if g == 1
                                    indx = 1:gap_indx(g);
                                elseif g == numel(gap_indx)+1
                                    indx = (gap_indx(g-1)+1):size(inverse_data,1);
                                else
                                    indx = (gap_indx(g-1)+1) : gap_indx(g);
                                end
                            end
                            time_o = inverse_data{indx,2};
                            azi    = inverse_data{indx,3};
                            elv_o  = inverse_data{indx,4};

                            if string(cur_sys) == "GLONASS"
                                continue
                            end
                            [M,a,b,cur_band] = Get_combined_observations(Meth_id, cur_sys, inverse_data, indx, cur_sat);

                            % Start inversion based on combined observations
                            nan_indx = isnan(M);
                            info = [elv_o,M,time_o];
                            info(nan_indx,:) = [];
                            if numel(info) == 0
                                continue
                            end
                            info = sortrows(info);
                            elv = info(:,1);
                            M = info(:,2);
                            time = info(:,3);
                            if numel(M) < 10
                                continue
                            end
                            sinelv = sind(elv);
                            if Meth_id == 5 % Remove the ionosphere error
                                p = polyfit(sinelv, M, 10);
                                m_fit = polyval(p, sinelv);
                                M = M-m_fit;
                            end
                            [refl_h, id, psd, pks] = mp2RH_lsp(sinelv, M, hell, hgtlim, a, b);
                            if isnan(pks)
                                continue
                            end

                            % Tropospheric correction
                            pre_h = refl_h(id);
                            hsfc = hell - pre_h;
                            if hsfc > hgtlim(1) && hsfc < hgtlim(2)
                                psfc = interp1(hgtlim,plim,hsfc,'linear');
                                tsfc = interp1(hgtlim,tlim,hsfc,'linear');
                                tmsfc = interp1(hgtlim,tmlim,hsfc,'linear');
                                esfc = interp1(hgtlim,elim,hsfc,'linear');

                                dele = (1/60) * 510 * psfc / ((9/5*tsfc+492) * 1010.16) .*cotd(elv+7.31./(elv+4.4));
                                sinelv = sind(elv + dele);
                                elv = elv+dele;

                                theta = elv;
                                curdt = datetime(tdatenum,'convertfrom','datenum');
                                curjd = juliandate(curdt);
                                thetarefr = [];
                                for jj = 1:numel(theta)
                                    tau = trop_delay_tp(curjd,sta_lat,hell,hsfc,theta(jj),pant,tmant,eant,ah,aw,lambda,psfc,tmsfc,esfc);
                                    thetarefr(jj) = asind( sind(theta(jj)) + 0.5*tau/pre_h );
                                end
                                sinelv = sind(thetarefr).';

                                [refl_h, id, psd, pks] = mp2RH_lsp(sinelv, M, hell, hgtlim, a, b);

                            end
                            cur_rh = refl_h(id);  % Final rh
                            trop_c = cur_rh - pre_h; % Tropospheric correction

                            valid = 1;
                            if  cur_rh > sta_asl-tide_range(1) || cur_rh < sta_asl-tide_range(2) || max(psd)<5 *mean(pks(1:end-1))
                                valid = 0;
                            end

                            if valid == 1
                                rid = rid + 1; % row id
                                Time = datetime(tdatenum + mean(time)/86400, 'ConvertFrom', 'datenum');       % datenum
                                System = string(cur_sys);
                                BAND = string(cur_band);
                                PRN = cur_sat;                               % Sat prn
                                ROC = tand(mean(elv)) / (( (pi/180) * (elv(end)-elv(1)) )...
                                    /(time(end)-time(1)));                   % tan(th)/dth/dt :rate of change
                                MIN_elv = min(elv);                          % THETA MIN
                                MAX_elv = max(elv);                          % THETA MAX
                                MEAN_AZI = nanmean(azi);                     % MEAN AZI
                                RH = cur_rh;                                 % L1 slvl
                                RH_info(rid,:) = table(Time,System,BAND,PRN,ROC,MIN_elv,MAX_elv,MEAN_AZI,RH,trop_c);
                            end
                        end

                    end
                end
            end
        end
        % save the RH_file
        RH_file = [settings.Out_path,'/RH_file/',settings.station_name,num2str(tdatenum),'RH_info.mat'];
        if ~exist([settings.Out_path,'/RH_file'], "dir")
            mkdir([settings.Out_path,'/RH_file'])
        end
        parsave(RH_file,RH_info,'RH_info')
    end
    delete(par)
end