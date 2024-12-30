function [t_rh,sea_level_ir, RMSE, Bias] = B_spline(RH_info_all, start_date, end_date, fig_set)

p = 2;
t = datenum(RH_info_all.Time);
rh= RH_info_all.RH;
xinit = t;
hinit = rh;
tanthter = RH_info_all.ROC;
kspac = 0.125;
tlen = end_date+1-start_date;
knots = [start_date*ones(1,p) ...
    start_date:kspac:start_date+tlen ...
    (start_date+tlen)*ones(1,p)]; % knot vector, multiplicity = 4
nsfac = tlen/kspac + p;
sfacs_0 = fig_set.sta_asl*ones(1,nsfac); % control points
tempfun_init = @(sfacs) bspline_spectral(sfacs, p, knots, tanthter, xinit,1)-hinit.';
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off'); % off for now but maybe should check sometimes
sfacs_init = lsqnonlin(tempfun_init,sfacs_0,[],[],options);

t_rh = start_date:(6/(24*60)):end_date+1;
rh_b_spline = bspline_deboor(p+1,knots,sfacs_init,t_rh);
sea_level_ir = fig_set.sta_asl - rh_b_spline - nanmean(fig_set.sta_asl - rh_b_spline);
t_rh = datetime(t_rh,'ConvertFrom','datenum');

time_tide      = fig_set.tidal_info{1};
sea_level_tide = fig_set.tidal_info{2};
if fig_set.option
    figure
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0.1 0.1 0.7 0.5]);
    tiledlayout(7,1,'TileSpacing','compact');
    nexttile([5,1])
    scatter(datetime(t,'ConvertFrom','datenum'), fig_set.sta_asl-rh-nanmean(fig_set.sta_asl-rh), ...
        30, [0.7, 0.7, 0.7],"filled", 'DisplayName', 'Raw retrievals')
    hold on
    plot(t_rh, sea_level_ir,'Color','r','LineWidth',1.5,...
        'DisplayName', 'Multi-GNSS combined retrievals')
    hold on 
    sea_level_tide = interp1(time_tide,sea_level_tide, t_rh);
    plot(t_rh,sea_level_tide,'Color','black','LineWidth',1.5,...
        'DisplayName', 'Tide gauge')
    %     xlabel('Time')
    set(gca,'XTickLabel',[])
    ylabel('Sea level (m)','FontWeight','bold')
    legend
    box on
    xlim([min(time_tide),max(time_tide)])
    title('B spline','FontWeight','bold','FontSize',18)

    nexttile([2,1])
    plot(t_rh,sea_level_ir- sea_level_tide)
    xlabel('Time','FontWeight','bold')
    datetick('x', 'dd-mmm-yyyy', 'keepticks', 'keeplimits')
    ylabel('Residual error (m)','FontWeight','bold')
    xlim([min(time_tide),max(time_tide)])
end

rh_b_spline = sea_level_ir;
Bias = nanmean(sea_level_ir - sea_level_tide);
RMSE = sqrt(nanmean((sea_level_ir - sea_level_tide).^2));
end