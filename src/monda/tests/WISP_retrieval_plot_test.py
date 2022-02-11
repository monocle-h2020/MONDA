import monda
# monda.import_submodules(monda)

print('running WISP test')
instrument= 'WISPstation001'
Date = '2019-06-21';
start_time ='9:00:00';
end_time = '14:15:00';
plot1=monda.WISP.data_analysis.WISPplots.plot_reflectance_station(instrument, Date, start_time, end_time);
plot1.savefig('/tmp/WISP_test_plot1.png')
plot2=monda.WISP.data_analysis.WISPplots.plot_radiances_station(instrument, Date, start_time, end_time);
plot2.savefig('/tmp/WISP_test_plot2.png')

exit()
