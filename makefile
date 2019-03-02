#=======================================================================
# set lib and include directoris
#=======================================================================
netcdf_path=/opt/apps/netcdf/version_3.6.3/intel_15
init_info_path=/home/guoyk/Program/wrfversion3/EnKF_LDA/include
#=======================================================================
# set name of the EXE
#=======================================================================

init_ens_exe=init_ens
anal_result_exe=anal_result
#obs_sim_exe=obs_sim
#ensrf_sim_exe=ensrf_sim
mean_io_exe=mean_io
#plot_sim_radar_exe=plot_sim_radar
#pert_sounding_exe=pert_sounding
ensrf_real_exe=ensrf_real
#extract_latlon_exe=extract_latlon
#sim_obs_exe=sim_obs
#=======================================================================
# set name of the obj
#=======================================================================

init_ens_obj    =  init_ens.o init_ens_sub.o lib.o io_ensrf.o anal_sub.o
anal_result_obj =  anal_result.o anal_sub.o io_ensrf.o lib.o
obs_sim_obj     =  obs_sim.o io_ensrf.o ensrf_lib.o obs_sim_sub.o radar_lib.o lib.o
ensrf_sim_obj   =  ensrf_sim.o ensrf_sim_sub.o ensrf_lib.o anal_sub.o obs_sim_sub.o \
                   io_ensrf.o radar_lib.o read_obs_sim.o innovation_rmse.o lib.o
mean_io_obj     =  mean_io.o io_ensrf.o anal_sub.o #ensrf_lib.o radar_lib.o
plot_sim_radar_obj = plot_sim_radar.o io_ensrf.o ensrf_lib.o obs_sim_sub.o radar_lib.o
pert_sounding_obj = pert_sounding.o lib.o
ensrf_real_obj  =  ensrf_real.o ensrf_real_sub.o ensrf_lib.o io_ensrf.o anal_sub.o obs_real_sub.o \
                    radar_lib.o read_obs_real.o innovation_rmse_real.o lib.o
extract_latlon_obj     =  extract_latlon.o io_ensrf.o
#sim_obs_soilvar_boj      = sim_obs_soilvar.o
#=======================================================================
#compile options
#=======================================================================

FFLAGS	       =	-O3 -convert big_endian  #-openmp
NETCDFLIB      =        -L$(netcdf_path)/lib -lnetcdf #-lm
NETCDFINCLUDE  =        -I$(netcdf_path)/include

#=======================================================================
.c.o:
	icc -c $<
.f90.o:
	ifort -c $(FFLAGS) $<
#=======================================================================
$(init_ens_exe):$(init_ens_obj)
	ifort $(FFLAGS)  -o init_ens.exe  $(init_ens_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(anal_result_exe):$(anal_result_obj)
	ifort $(FFLAGS)  -o anal_result.exe $(anal_result_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(obs_sim_exe):${obs_sim_obj}
	ifort $(FFLAGS)  -o obs_sim.exe  $(obs_sim_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(ensrf_sim_exe):$(ensrf_sim_obj)
	ifort $(FFLAGS)  -o ensrf_sim.exe  $(ensrf_sim_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(mean_io_exe):$(mean_io_obj)
	ifort $(FFLAGS)  -o mean_io.exe  $(mean_io_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(plot_sim_radar_exe):$(plot_sim_radar_obj)
	ifort $(FFLAGS)  -o plot_sim_radar.exe  $(plot_sim_radar_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(pert_sounding_exe) : $(pert_sounding_obj)
	ifort $(FFLAGS)  -o pert_sounding.exe $(pert_sounding_obj)
$(ensrf_real_exe):$(ensrf_real_obj)
	ifort $(FFLAGS)  -o ensrf_real.exe  $(ensrf_real_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(extract_latlon_exe):$(extract_latlon_obj)
	ifort $(FFLAGS) -o extract_latlon.exe	 $(extract_latlon_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
$(sim_ob_exe):$(sim_ob_obj)
	ifort -o sim_obs_soilvar.exe	 $(sim_ob_obj) $(NETCDFLIB) $(NETCDFINCLUDE)
#=======================================================================
init_ens.o : init_ens.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
init_ens_sub.o : init_ens_sub.f90
	ifort  -c $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
io_ensrf.o : io_ensrf.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
lib.o : lib.f90
	ifort  -c $(FFLAGS) $<
anal_result.o : anal_result.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
anal_sub.o : anal_sub.f90
	ifort  -c $(FFLAGS) $<
obs_sim.o : obs_sim.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
ensrf_lib.o : ensrf_lib.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $<
obs_sim_sub.o : obs_sim_sub.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $<
ensrf_sim.o : ensrf_sim.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
ensrf_sim_sub.o : ensrf_sim_sub.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $<
mean_io.o : mean_io.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
radar_lib.o : radar_lib.f90
	ifort  -c $(FFLAGS) $<
plot_sim_radar.o : plot_sim_radar.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
pert_sounding.o : pert_sounding.f90 $(init_info_path)/pert_sounding.inc
	ifort  -c $(FFLAGS) -I$(init_info_path) $<
read_obs_sim.o : read_obs_sim.f90 $(init_info_path)/namelist.inc
	ifort  -c $(FFLAGS) -I$(init_info_path) $<
innovation_rmse.o : innovation_rmse.f90	$(init_info_path)/namelist.inc
	ifort  -c $(FFLAGS) -I$(init_info_path) $<
ensrf_real.o : ensrf_real.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
ensrf_real_sub.o : ensrf_real_sub.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $<
obs_real_sub.o : obs_real_sub.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $<
read_obs_real.o : read_obs_real.f90 $(init_info_path)/namelist.inc
	ifort  -c $(FFLAGS) -I$(init_info_path) $<
innovation_rmse_real.o : innovation_rmse_real.f90 $(init_info_path)/namelist.inc
	ifort  -c $(FFLAGS) -I$(init_info_path) $<
extract_latlon.o : extract_latlon.f90 $(init_info_path)/namelist.inc
	ifort  -c -I$(init_info_path) $(FFLAGS) $(NETCDFLIB) $(NETCDFINCLUDE) $<
#sim_obs_soilvar.o: sim_obs_soilvar.f90
#	ifort  -c  $(NETCDFLIB) $(NETCDFINCLUDE) $<
all:
	make init_ens ensrf_real mean_io
clean:
	rm *.exe *.o
