function Analysis_Callback(app)

%%
% set tides data
tide_file = app.TidefileEditField.Value;
Tide = load(tide_file);
tide = [Tide.xaxis, Tide.slvl];
[~, idx] = unique(tide(:,1), 'stable');
tide = tide(idx, :);
xaxis = tide(:,1);
slvl  = tide(:,2);

% set Final file
Final_path = app.FilepathEditField.Value;
Final_files= app.FinalFileListTextArea.Value;
% loop by files
for f = 1:app.FilenumberEditField.Value
    %% Load file & get information
    % get time from file name
    file_info = strsplit(Final_files{f},'_');
    start_date = datenum(file_info{2});
    end_date   = datenum(file_info{3});

    % get RH_info_all
    if string(file_info{4}) == "SNR-Inverse.mat"
        FF = load(fullfile(Final_path,Final_files{f}));
        RH_info_all = FF.Final_info;
        RH_info_all(isnan(RH_info_all{:,"RH_Inverse"}),:) = [];
    else
        FF = load(fullfile(Final_path,Final_files{f}));
        Final_info = FF.Final_info;
        fieldNames = fieldnames(Final_info);
        RH_info_all = [];
        for i = 1:length(fieldNames)
            currentTable = Final_info.(fieldNames{i});
            RH_info_all = [RH_info_all; currentTable];
        end
    end

    % get tide info
    tide_indx = xaxis>=start_date & xaxis<=end_date+1;
    time_tide_raw = datetime(xaxis(tide_indx),'ConvertFrom','datenum');
    sea_level_tide = slvl(tide_indx);

    % interp
    time_tide = datetime(start_date,'ConvertFrom','datenum'):hours(1):datetime(end_date+1,'ConvertFrom','datenum');
    tide = [datenum(time_tide_raw) sea_level_tide];
    [~, idx] = unique(tide(:,1), 'stable');
    tide = tide(idx, :);
    xaxis = tide(:,1);
    slvl  = tide(:,2);
    sea_level_tide = interp1(datetime(xaxis,'ConvertFrom','datenum'), slvl, time_tide);

    %% Figure base settings
    colors = [
        0.3, 0.6, 0.85;
        0.9, 0.5, 0.3;
        0.95, 0.75, 0.4;
        0.65, 0.4, 0.7;
        0.6, 0.8, 0.4;
        0.5, 0.8, 0.9;
        0.8, 0.4, 0.4;
        ];
    set(groot, 'DefaultLineLineWidth', 0.5);
    set(groot, 'DefaultAxesFontName', 'Times New Roman');
    set(groot, 'DefaultAxesFontSize', 14);
    set(groot, 'DefaultTextFontName', 'Times New Roman');

    %% Analysis of different Methdos
    % SNR-Spectral method
    if string(file_info{4}) == "SNR-Spectral.mat"
        % Excluded Bad frequency points
        for exf = 1:numel(app.ExcludedfrequencypointsTextArea.Value)
            out_ind = RH_info_all{:,3}==string(app.ExcludedfrequencypointsTextArea.Value{exf});
            RH_info_all(out_ind,:) = [];
        end

        Result_display_SNR_Spectral(app, RH_info_all,time_tide,sea_level_tide, colors,start_date,end_date)
        % save picture
        name = [app.AnalysisResultsEditField.Value,'/',...
            app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral'];
        saveas(gcf, [name,'.fig']);
        if app.CorrelationscatterCheckBox.Value
            figure
            ir_t = RH_info_all{:,"Time"};
            ir_h = app.AntennaheightmEditField.Value - RH_info_all{:,"RH"};
            cor_scatter(ir_h, datenum(ir_t), sea_level_tide, datenum(time_tide))
            ylabel('Tide gauge (m)')
            xlabel('GNSS-IR raw inversion value (m)')
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral_scatter'];
            saveas(gcf, [name,'.fig']);
        end

        if app.RobustregressionstrategyCheckBox.Value % Robust regression strategy
            fig_set.option  = 1;
            fig_set.sta_asl = app.AntennaheightmEditField.Value;
            fig_set.tidal_info = {time_tide, sea_level_tide};
            [t_rrs,rh_rrs, RMSE_rrs, Bias_rrs] = rrs(RH_info_all, start_date, end_date, fig_set);
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral_rrs'];
            saveas(gcf, [name,'.fig']);
            save([name,'.mat'],"t_rrs","rh_rrs", "RMSE_rrs", "Bias_rrs")
            if app.CorrelationscatterCheckBox.Value
                figure
                cor_scatter(rh_rrs', datenum(t_rrs'), sea_level_tide, datenum(time_tide))
                ylabel('Tide gauge (m)')
                xlabel('GNSS-IR Robust regression strategy (m)')
                name = [app.AnalysisResultsEditField.Value,'/',...
                    app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral_rrs_scatter'];
                saveas(gcf, [name,'.fig']);
            end
        end

        if app.BsplineCheckBox.Value % B-spline
            fig_set.option  = 1;
            fig_set.sta_asl = app.AntennaheightmEditField.Value;
            fig_set.tidal_info = {time_tide, sea_level_tide};
            [t_b_spline,rh_b_spline, RMSE_b, Bias_b] = B_spline(RH_info_all, start_date, end_date, fig_set);
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral_b_spline'];
            saveas(gcf, [name,'.fig']);
            save([name,'.mat'],"t_b_spline","rh_b_spline", "RMSE_b", "Bias_b")
            if app.CorrelationscatterCheckBox.Value
                figure
                cor_scatter(rh_b_spline', datenum(t_b_spline'), sea_level_tide, datenum(time_tide))
                ylabel('Tide gauge (m)')
                xlabel('GNSS-IR B-spline (m)')
                name = [app.AnalysisResultsEditField.Value,'/',...
                    app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_SNR-Spectral_b_spline_scatter'];
                saveas(gcf, [name,'.fig']);
            end
        end

    % SNR-Inverse modeling
    elseif string(file_info{4}) == "SNR-Inverse.mat"
        time = RH_info_all.Time;
        RH_inverse_modeling = -app.AntennaheightmEditField.Value + RH_info_all.RH_Inverse;
        RH_inverse_modeling = RH_inverse_modeling-nanmean(RH_inverse_modeling);
        RH_spectral_b = -app.AntennaheightmEditField.Value + RH_info_all.RH_Spectral_bsp;
        RH_spectral_b = RH_spectral_b-nanmean(RH_spectral_b);
        interp_tide = interp1(time_tide,sea_level_tide,time,"linear");
        figure
        tiledlayout(3,6,"TileSpacing","compact")
        nexttile([2,4])
        plot(time, RH_spectral_b, 'r--', DisplayName='Spectral Retrieval')
        hold on
        plot(time, RH_inverse_modeling, 'blue',DisplayName='Inverse Modeling')
        plot(time, interp_tide, 'black', DisplayName='Tide gauge')
        xlabel('Time','FontWeight','bold')
        ylabel('Sea level (m)','FontWeight','bold')
        xlim([min(time),max(time)+1])
        legend

        % cor-scatter
        nexttile(5,[1,2])
        cor_scatter(RH_inverse_modeling, datenum(time), sea_level_tide, datenum(time_tide))
        ylabel('Tide gauge (m)')
        xlabel('Inverse Modeling (m)')
        nexttile(11,[1,2])
        cor_scatter(RH_spectral_b, datenum(time), sea_level_tide, datenum(time_tide))
        ylabel('Tide gauge (m)')
        xlabel('Spectral Retrieval (m)')

        % RMSE
        nexttile(17,[1,2])
        RMSE(1) = sqrt(nanmean((interp_tide-RH_spectral_b).^2));
        RMSE(2) = sqrt(nanmean((interp_tide-RH_inverse_modeling).^2));
        b = barh(RMSE*100,0.5, 'FaceColor',[0.3010 0.7450 0.9330],'EdgeColor',[0.3010 0.7450 0.9330],'BarWidth',0.9, 'LineWidth', 2);
        set(gca, 'YTickLabel', {'Spectral Retrieval','Inverse Modeling'})
        xlabel('RMSE(cm)','FontWeight','bold')
        xtips1 = b(1).YEndPoints + 0.01;
        ytips1 = b(1).XEndPoints;
        rms_band = string(b(1).YData);
        text(xtips1,ytips1,rms_band,'HorizontalAlignment','right','Color',[1, 0.2, 0.2],'FontWeight','bold')
        xlim([0, max(RMSE)*100+5])

        % res-his
        for m = 1:2
            if m == 1
                data = interp_tide-RH_spectral_b;
                nexttile(13,[1,2])
            else
                data = interp_tide-RH_inverse_modeling;
                nexttile(15,[1,2])
            end
            h = histfit(data, 50, 'normal');
            h(1).FaceColor = [0.4, 0.6, 0.8];
            h(1).EdgeColor = [1,1,1];
            mu = mean(data);
            sigma = std(data);
            hold on;
            xline(mu, '--r', 'LineWidth', 2);
            text(mu, max(ylim)*0.9, sprintf('Mean:%.4f', mu),    'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'right');
            text(mu, max(ylim)*0.65, sprintf('STD :%.4f', sigma), 'Color', 'b', 'FontSize', 10, 'HorizontalAlignment', 'right');
            if m == 1
                ylabel('number')
                xlabel('Spectral Residual')
            else
                xlabel('Inverse Modeling Residual')
            end
        end

        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        name = [app.AnalysisResultsEditField.Value,'/',...
            app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4)];
        saveas(gcf, [name,'.fig']);

        % OBS-carrier / pseudo range
    else
        obs_method = unique(RH_info_all.BAND);
        Result_display_OBS(RH_info_all, app, time_tide, sea_level_tide, end_date, start_date, colors);
        % save picture
        name = [app.AnalysisResultsEditField.Value,'/',...
            app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4)];
        saveas(gcf, [name,'.fig']);
        if app.CorrelationscatterCheckBox.Value
            figure
            ir_t = RH_info_all{:,"Time"};
            ir_h = app.AntennaheightmEditField.Value - RH_info_all{:,"RH"};
            ir_h = ir_h-nanmean(ir_h);
            cor_scatter(ir_h, datenum(ir_t), sea_level_tide, datenum(time_tide))
            ylabel('Tide gauge (m)')
            xlabel('GNSS-IR raw inversion value (m)')
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4),'_scatter'];
            saveas(gcf, [name,'.fig']);
        end

        if app.RobustregressionstrategyCheckBox.Value % Robust regression strategy
            fig_set.option  = 1;
            fig_set.sta_asl = app.AntennaheightmEditField.Value;
            fig_set.tidal_info = {time_tide, sea_level_tide};
            [t_rrs, rh_rrs, RMSE_rrs, Bias_rrs] = rrs(RH_info_all, start_date, end_date, fig_set);
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4),'_rrs'];
            saveas(gcf, [name,'.fig']);
            save([name,'.mat'],"t_rrs","rh_rrs", "RMSE_rrs", "Bias_rrs")
            if app.CorrelationscatterCheckBox.Value
                figure
                cor_scatter(rh_rrs', datenum(t_rrs'), sea_level_tide, datenum(time_tide))
                ylabel('Tide gauge (m)')
                xlabel('GNSS-IR Robust regression strategy (m)')
                name = [app.AnalysisResultsEditField.Value,'/',...
                    app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4),'_rrs_scatter'];
                saveas(gcf, [name,'.fig']);
            end
        end

        if app.BsplineCheckBox.Value % B-spline
            fig_set.option  = 1;
            fig_set.sta_asl = app.AntennaheightmEditField.Value;
            fig_set.tidal_info = {time_tide, sea_level_tide};
            [t_b_spline,rh_b_spline, RMSE_b, Bias_b] = B_spline(RH_info_all, start_date, end_date, fig_set);
            name = [app.AnalysisResultsEditField.Value,'/',...
                app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4),'_b_spline'];
            saveas(gcf, [name,'.fig']);
            save([name,'.mat'],"t_b_spline","rh_b_spline", "RMSE_b", "Bias_b")
            if app.CorrelationscatterCheckBox.Value
                figure
                cor_scatter(rh_b_spline', datenum(t_b_spline'), sea_level_tide, datenum(time_tide))
                ylabel('Tide gauge (m)')
                xlabel('GNSS-IR B-spline (m)')
                name = [app.AnalysisResultsEditField.Value,'/',...
                    app.StationnameEditField.Value, '_',file_info{2},'_',file_info{3},'_',file_info{4}(1:end-4),'_b_spline_scatter'];
                saveas(gcf, [name,'.fig']);
            end
        end
    end


end
