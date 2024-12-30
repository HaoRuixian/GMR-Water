function main_batch(settings)

%% Settings
start_time   = settings.Time(1);
end_time     = settings.Time(2);
start_date   = datenum(start_time);
end_date     = datenum(end_time);

%% file preparation
load(settings.Station_info_file)

% extract SNR & carrier pseudorange from rinex
if settings.flow(1)
    RINEX_info_EXTRCATION
end

% Inverse
if settings.flow(2)
    Inverse_RH
end

% Tidal correction & save
if settings.flow(3)
    Tidal_correction_save
end
end
