function = Build_Pseudo(Te,z)
  %Te is the transfer Matrix, z is the specific energy value that we are working with,
  %This is to name the plots

    resolution = 400;%Resolution
    output_size = [3000 1800];%Size in pixels
    save = "no"

    opts.npts=100;
    %opts.levels=-7.5:0.5:-2;
    %opts.ax=[-1 0.2 -1.05 0];
    %opts.proj_lev=Inf;
    %opts.colour=1;
    %opts.thick_lines=1;
    %opts.scale_equal=1;
    %opts.grid=0;
    %opts.dim=1;
    %opts.no_ews=0;
    %opts.no_psa=0;
    opts.fov=1;
    %opts.direct=1;
    %opts.print_plot=1;
    eigtool(Te,opts)
    if save == "yes"
        set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
        print(compose("T_%d.png",z),'-dpng',['-r' num2str(resolution)]);
        saveas(fig, compose("T_%",z), 'png')
    elseif save == "no"
        
end