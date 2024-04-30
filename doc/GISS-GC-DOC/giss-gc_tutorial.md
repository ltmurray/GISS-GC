# Running the GISS ModelE2.1 GCM to drive the GEOS-Chem CTM 
Lee T. Murray
Last Updated: April 30, 2024

## Overview
[ModelE2.1](https://www.giss.nasa.gov/tools/modelE/) is a global general circulation model (GCM) developed primarily at the NASA Goddard Institute for Space Studies in New York, NY. Version "E2.1" of the model is the version that was frozen for the majority of NASA’s contributions to the Coupled Model Intercomparison Project Phase 6 (CMIP6) experiment done in support of the Sixth IPCC Assessment Report [(AR6)](https://www.ipcc.ch/assessment-report/ar6/). 

The emphasis of this tutorial is getting ModelE2.1 to compile and run on your own UNIX-based operating system (MacOS, Linux), and produce the diagnostics necessary to drive the [GEOS-Chem](http://www.geos-chem.org) Chemical Transport Model.

The tutorial was tested using the Ubuntu Desktop 22.04.4 LTS Linux operating system.

If you are new to Linux/UNIX systems, first read the following guides
- [UNIX tutorial for beginners](https://linuxclass.heinz.cmu.edu/doc/Unix-Tutorial-surrey/)
- [How to edit files with emacs](http://www.jesshamrick.com/2012/09/10/absolute-beginners-guide-to-emacs/)

## Table of Contents
1. [First-Time Setup](#first)
2. [Setup a ModelE2.1 run](#config)
3. [Perform a ModelE2.1 simulation](#run)
4. [Process ModelE2.1 output](#postprocess)

## 1. First-Time Setup<a name="first"></a>

### 1a. Download GISS-GC

This tutorial will assume that the source code is located in the following directory: `$GISS_HOME`

```console
git clone git@github.com:fetch4/GISS-GC.git $GISS_HOME
```

### 1b. Install Software Dependencies

The most important dependency is the [netcdf-fortran](https://github.com/Unidata/netcdf-fortran) library compiled with MPI parallelization. If using a package manager like [miniforge](https://github.com/conda-forge/miniforge) (_recommended_), this will come with all other necessary packages for running ModelE2.1, although you will also want cmake, git and a text editor.

Example environmental yaml files for miniforge are available in `$GISS_HOME/doc/GISS-GC-DOC/envs`, and may be installed using
```console
conda env create -n giss -f giss_ubuntu_22.04.4_lts.yml
```

To load the environment containing all the variables necessary to compile and run ModelE2.1
```console
conda activate giss
```

### 1c. Install GISS Post-Processing Tools<a name="scaleacc"></a>

The tools located in `$GISS_HOME/model/mk_diags` are necessary for post-processing the model output into a useful format. I recommend copying that directory to a common location to which you can point your `$PATH` variable to the binaries.

```console
mkdir $HOME/tools
cp -rp $GISS_HOME/model/mk_diags $HOME/tools/
cd $HOME/tools/mk_diags
```

and then compile the tools (be sure to have activated your conda environment first)
```console
/bin/bash compscr
```

You may point your $PATH to this directory by either executing the following command or adding it to your `$HOME/.bashrc` (Linux) or `$HOME/.zshrc` (macOS) file.
```console
export PATH=$HOME/tools/mk_diags:$PATH
```

### 1d. Setup GISS Environment Variables

The GCM uses several variables that are defined in the `$HOME/.modelErc` configuration file, which we will set up now. You will need to pass the script `$ModelE_Support`, which will be the path to a directory that will serve as the root location for where model input and output are stored.
```console
cd $GISS_HOME/decks
make config COMPILER=gfortran ModelE_Support=$ModelE_Support OVERWRITE=YES
```

This will generate the file `$HOME/.modelErc`. Replace the default options shown in red with the green changes (ignore the leading - and +). Note, `$CONDA_PREFIX` is set when you activate your conda environment, and points to where all the software dependencies were installed. You will need make sure that `$ModelE_Support` and `$CONDA_PREFIX` are given as actual paths in the file.
```diff
# This file contains global options for modelE. 
# By default it assumes that the directory structure for modelE runs
# is set under $ModelE_Support .

## Directory structure ##

# DECKS_REPOSITORY - a directory for permanenet storage of run info.
# All rundecks that you create will be copied to this directory. 
DECKS_REPOSITORY=$ModelE_Support/prod_decks

# CMRUNDIR - directory to which all run directories will be linked.
# This directory will be searched by most scripts for locations of 
# specific runs.
CMRUNDIR=$ModelE_Support/prod_runs

# GCMSEARCHPATH - directory to search for gcm input files.
# All necessary input files should be copied or linked to this directory.
GCMSEARCHPATH=$ModelE_Support/prod_input_files

# EXECDIR - path to directory with modelE scripts and with some
# executables. This directory should contain the scripts from modelE/exec.
EXECDIR=$ModelE_Support/exec

# SAVEDISK - a directory where all run directories (which will contain
# all output files such as rsf, acc etc.) will be created. This should
# be big enough to accomodate all model output.
SAVEDISK=$ModelE_Support/huge_space

## External libraries ##

# Some of these options can be provided by environment modules (if you 
# use them). Specify here only what is necessary. Options specified 
# here will overwrite options proviided by environment modules.

# NETCDFHOME - path to location of netcdf installation directory. 
-# NETCDFHOME=/opt/netcdf/3.6.3
+NETCDFHOME=$CONDA_PREFIX

# MPI - set to YES if you want to compile the model for parallel 
# execution on multiple CPU cores. Keep in mind, that functional 
# MPI library should be installed on your computer and its type 
# and location should be specified below.
# This option can be overwritten from the compile line.
-MPI=NO
+MPI=YES

# MPIDISTR - the MPI distribution you are using. Currently supported 
# distributions are: 'intel, 'openmpi', 'mpich2', 'mvapich2', 'SCALI',
# 'mpt' 
-# MPIDISTR=openmpi
+MPIDISTR=openmpi

# MPIDIR - path to the MPI installation directory. (Needs to be set
# only if compiler can't find correct MPI library and include files by
# default)
-# MPIDIR=/opt/openmpi
+MPIDIR=$CONDA_PREFIX

# MPILIBDIR - path to the location of MPI library. Set it only if 
# it is different from the default $MPIDIR/lib
# MPILIBDIR=/opt/openmpi/lib

# MPIINCLUDEDIR - path to location of MPI include files. Set it only
# if it is different from the default $MPIDIR/include
# MPIINCLUDEDIR=/opt/openmpi/include

# ESMF5_DIR - path to the installation directory of ESMF (version 5)
# library. (Required only for Cubed Sphere simulations)
# ESMF5_DIR=

# ESMF_BOPT - optimization level of ESMF library. (Should only be used
# togeteher with ESMF5_DIR)
# ESMF_BOPT=O
+ESMF=NO

## Architecture and compiler

# ABI - Application Binary Interfaces. This variable specifies the
# architecture you are using. The valid values are '64' and '32'. 
# On most modern systems you should use '64'. Use '32' if your
# hardware or compiler support only 32-bit binaries.
ABI=64

# COMPILER - specifies the Fortran compiler you are using. Currently
# only 'intel' and 'gfortran' are supported. ('nag' has partial
# support on development branch.) If you are using Modules for
# Environment Management, then this variable may already be set in the
# environment. In this case you don't need to set it here.
-# COMPILER=gfortran
COMPILER=gfortran

## General User Preferences ##

# MAILTO - email address of the user. When the program ends/crashes
# all notifications will be sent to this address. Leave empty 
# or unset if you don't want to receive these emails
MAILTO=

# UMASK - the value of 'umask' you want to use for model runs. The files
# inside the run directory will have permissions set according to this
# mask.
UMASK=002

# OVERWRITE - can "gmake rundeck" overwrite files already in repository?
# (i.e. in the directory DECKS_REPOSITORY)
OVERWRITE=NO

# OUTPUT_TO_FILES - if set to YES all errors and warnings will be sent
# to files with the names <source_name>.ERR
OUTPUT_TO_FILES=NO

# VERBOSE_OUTPUT - if set to YES gmake will show compilation commands
# and some other information. Otherwise most of the output will be
# suppressed
VERBOSE_OUTPUT=YES
```

### 1e. Other Helpful Tools

There are many useful tools that may be used to manipulate [NetCDF](https://docs.unidata.ucar.edu/nug/current/netcdf_introduction.html) files. These were installed with the software environment above.

- [Climate Data Operators (cdo)](https://code.mpimet.mpg.de/projects/cdo)
- [NetCDF Common Operators (nco)](https://nco.sourceforge.net)
- [ncview](https://cirrus.ucsd.edu/ncview/)

Also, NASA has developed the easy-to-use [Panoply](https://www.giss.nasa.gov/tools/panoply/) tool for easily exploring NetCDF files, which requires a Java Runtime Environment be installed on your system.

## Setup a ModelE2.1 run<a name="config"></a>

Now that your dependencies are installed and software environment is set up, you can prepare a model simulation.

### 2a. Prepare Rundeck

Go to the code directory, and enter the decks directory, which holds the “rundeck” files (*.R) that define a given simulation.

```console
cd $GISS_HOME/decks
```
Let us generate a rundeck based on the `E6F40` template
```console
make rundeck RUN=E6F40_TEST RUNSRC=E6F40 OVERWRITE=YES
```
This creates a new simulation with the name assigned to `RUN` (it should be a unique name for every simulation you perform), using as a template the `RUNSRC` template. Here "E6F40" means a CMIP6-class atmosphere-only simulation with fine (“F”; 2º x 2.5º) horizontal resolution and 40 vertical layers from the surface to the mesosphere, and the default configuration is for a preindustrial simulation ca. 1850s in which all composition is prescribed (i.e., no interactive chemistry).

Let’s modify the rundeck we generated (`$GISS_HOME/decks/E6F40_TEST.R`) to run for one month. The row with YEARI indicates when the model starts; the row in red shows when the model will end.  Note: these are dummy years in this setup as master_yr is set to 1850 earlier in the file, which makes all boundary conditions and climate forcers equal to 1850 climatology.

Setting KDIAG equal to 13*0 reduces the amount of text archived in the log file.

```diff
&INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
-YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
+YEARE=1950,MONTHE=01,DATEE=1,HOURE=0,     KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
```
### 2b. Download Input Files

No we need to download the input files that are specified in the rundeck file `E6F40_TEST.R`. In practice, it is easiest to try the next step, which will fail but provide a list of the missing files you need.

The input files are available for download from NASA at [https://portal.nccs.nasa.gov/GISS_modelE/](https://portal.nccs.nasa.gov/GISS_modelE/).

The following will download the necessary input files using the ``wget`` command line tool for the `E6F40_TEST` simulation and put them into the directory in which the model will be looking.

```bash
#/bin/bash
cd $ModelE_Support/prod_input_files
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/NCARIC.144x90.D7712010_ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/GIC.144X90.DEC01.1.ext_1.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/OST_144x90.1876-1885avg.CMIP6.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/SICE_144x90.1876-1885avg.CMIP6.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/ZSIfac_144x90.1876-1885avg.CMIP6.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/Z2HX2fromZ1QX1N.BS1.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/RD_Fd.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/RD_Fd.names.txt
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/CD144X90.ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/V144x90_EntMM16_height_trimmed_scaled_ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/V144x90_EntMM16_lai_trimmed_scaled_ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/Irrig144x90_1848to2100_FixedFuture_v3.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/S144X900098M.ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/top_index_144x90_a.ij.ext.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/ZVAR2X25A.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/soil_textures_top30cm_2x2.5
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/soilcarb_top30cm_2x2.5.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/GLMELT_144X90_gas.OCN.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/sgpgxg.table8
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/LWTables33k_lowH2O_CO2_O3_planck_1-800
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/LWCorrTables33k
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/H2Ocont_MT_CKD
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/miescatpar.abcdv2
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/oct2003.relhum.nr.Q633G633.table
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/STRATAER.VOL.1850-2014_CMIP6_hdr
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/cloud.epsilon4.72x46
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/solar.CMIP6official.ann1850-2299_with_E3_fastJ.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/topcld.trscat8
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/ISCCP.tautables
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/GHG.CMIP6.1-2014.txt
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/CO2profile.Jul16-2017.txt
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/dH2O_by_CH4_monthly
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/o3_2010_shindell_144x90x49_April1850.nc
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/MSU_SSU_RSS_weights.txt
wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/REG2X2.5
for spec in BCdalbsn DUST SUL SSA NIT OCA BCA BCB O3; do
  mkdir -p $ModelE_Support/prod_input_files/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/$spec
  cd $ModelE_Support/prod_input_files/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/$spec
  for year in 1854 1864 1874 1884 1894 1904 1914 1924 1934 1944 1954 1964 1974 1984 1994 2004; do
    wget https://portal.nccs.nasa.gov/GISS_modelE/modelE_input_data/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/$spec/$year.nc
  done
done
```

### 2c. Compile the Model and Setup the Run Directory

Now we can compile the model.
```console
cd $GISS_HOME/decks
make -j setup RUN=E6F40_TEST F90=mpif90
```

If the model successfully compiled and set up its run directory, then you should see a message similar to the following:
```text
===> linking
mpif90 main.o -m64  DIAG_RES_F.o FFT144.o IO_DRV.o ATMDYN.o MOMEN2ND.o QUS_DRV.o QUS3D.o STRATDYN.o STRAT_DIAG.o GEOM_B.o DIAG_ZONAL.o GCDIAGb.o DIAG_PRT.o POUT.o MODEL_COM.o MODELE_DRV.o MODELE.o ATM_COM.o ATM_DRV.o ATM_UTILS.o QUS_COM.o QUSDEF.o SURFACE.o SURFACE_LANDICE.o FLUXES.o GHY_COM.o GHY_DRV.o VEG_DRV.o ENT_DRV.o ENT_COM.o PBL_COM.o PBL_DRV.o PBL.o IRRIGMOD.o ATURB.o LAKES_COM.o LAKES.o SEAICE.o SEAICE_DRV.o LANDICE.o LANDICE_COM.o LANDICE_DRV.o ICEDYN_DRV.o ICEDYN.o RAD_COM.o RAD_DRV.o RADIATION.o RAD_UTILS.o ALBEDO.o READ_AERO.o ocalbedo.o DIAG_COM.o DIAG.o DEFACC.o OCN_DRV.o OCEAN.o OCNML.o  Atm144x90.o AtmLayering.o AtmRes.o ATMDYN_COM.o CLOUDS2.o CLOUDS2_DRV.o CLOUDS_COM.o   \
  giss_LSM/libgiss_LSM.a dd2d/libdd2d.a solvers/libsolvers.a Ent/libEnt.a MPI_Support/libMPI_Support.a shared/libshared.a profiler/libprofiler.a -L/home/ltmurray/miniforge3/envs/giss/lib -lnetcdf -L/home/ltmurray/miniforge3/envs/giss/lib -L/opt/local/lib -lnetcdff -lnetcdf -L/home/ltmurray/miniforge3/envs/giss/lib -lmpi_mpifh -lmpi -lrt -o E6F40_TEST.bin  
===> linking ok

make[2]: Leaving directory '$GISS_HOME/model'
make[1]: Leaving directory '$GISS_HOME/model'
mv $GISS_HOME/model/E6F40_TEST.bin $GISS_HOME/decks/E6F40_TEST_bin/E6F40_TEST.exe
----    Looks like compilation finished successfully     ----
make setup_script RUN=E6F40_TEST
make[1]: Entering directory '$GISS_HOME/decks'
----- Saving Rundeck and other info to global repository -----
---------        Starting setup for E6F40_TEST        ----------
--------------------------------------------------------------
Using settings from  
CMRUNDIR = $ModelE_Support/prod_runs
GCMSEARCHPATH = $ModelE_Support/prod_input_files
SAVEDISK = $ModelE_Support/huge_space
MAILTO = 
UMASK = 002
MPIRUN_COMMAND = 
setting up run E6F40_TEST
output files will be saved in $ModelE_Support/huge_space/E6F40_TEST
---------       Setup  finished successfully         ---------
make[1]: Leaving directory '$GISS_HOME/decks'
--- The directory E6F40_TEST was created and prepared for the run ---
---                                                        ---
--- To start the first hour run from the initial conditions --
--- you may want to execute:                               ---
---                                                        ---
   ../exec/runE E6F40_TEST -cold-restart -np <NPROC>          
---                                                        ---
--- or use your favorite script which does the same        ---
--- You can use "make setup-run ... " to combine setup     ---
--- and first hour run.                                    ---

```

### 2d. Run the First Hour of the Model

When a new simulation begins, we first like to test to make sure that the model can get through an entire time step successfully.

Following successful compilation, a symbolic link to the run directory should appear, which is located in the `$ModelE_Support/huge_space` directory that you defined in your `$HOME/.modelErc`. 

Change into that directory. Here and for the remainder of the tutorial, `$RUNID` is the unique name of your run (E6F40_TEST if you’ve been following from above). 

```console
cd $ModelE_Support/huge_space/$RUNID
```

Now create the links to all in the necessary input files for the simulation
```console
/bin/bash ${RUNID}ln
```

Now you can execute the first hour to make sure everything runs okay. This means run in parallel with 6 processors (`mpirun -np 6`) the program `${RUNID}.exe`, and feed it the input file (`-i`) named `I`, re-start from the very beginning of the simulation (`-cold-restart`), and save the output into log file `cold_restart.log`. 
```console
mpirun -np 6 ./${RUNID}.exe -i I -cold-restart &> cold_restart.log &
```
Here we are asking for 6 processors, but you can ask for any number that your system has up to 88 to make the simulation faster.

If all goes well, you should see the model run for the first hour nad then gracefully stop.
You can follow along as the model outputs `cold_restart.log` using the `tail` command:
```console
tail -f cold_restart.log
```

If the model completes successfully, the log file will say
```text
 PROGRAM TERMINATED NORMALLY - Internal clock time:   16034


 ************************************************************************************************************************************
 ************************************************************************************************************************************
 ************************************************************************************************************************************
 ************************************************************************************************************************************

-------------------------------------------------------------------------------------------------------
        Name         Inclusi  Inclusive   Exclusive  Trips   Maximum   MAX   Minimum   MIN   Average
                        %      min/day     min/day           seconds         seconds         seconds
-------------------------------------------------------------------------------------------------------
main                  100.00     6.88768     5.10859      6  17.219318    0  17.219035    1   2.128581
main                   25.83     1.77909     0.41510     12   4.448063    4   4.446606    0   0.086479
Main Loop               3.81     0.26238     0.26238     12   0.659168    5   0.649636    1   0.054663
 Atm. Dynamics          0.00     0.00008     0.00008     24   0.000239    4   0.000156    0   0.000008
 MELT_SI()              1.02     0.07058     0.07058     12   0.177643    1   0.175558    5   0.014703
 CONDSE()               9.07     0.62503     0.62503     12   1.722705    1   1.262316    5   0.130214
 RADIA()                1.91     0.13146     0.00007     12   0.329430    4   0.326927    0   0.000015
 Surface                0.00     0.00000     0.00000      0   0.000000    0   0.000000    0   0.000000
  Precip                1.91     0.13139     0.12241     12   0.329306    4   0.326612    0   0.025502
  SURFACE()             0.05     0.00361     0.00361     12   0.009411    0   0.008481    4   0.000752
 DYNSI()                0.00     0.00006     0.00006     24   0.000241    5   0.000100    2   0.000006
 UNDERICE()             0.00     0.00030     0.00030     24   0.001950    5   0.000161    2   0.000031
 GROUND_SI()            0.00     0.00002     0.00002     12   0.000172    0   0.000014    3   0.000004
 GROUND_LI()            0.00     0.00008     0.00008     12   0.000366    4   0.000080    0   0.000016
 GROUND_LK()            0.01     0.00036     0.00036     12   0.001267    1   0.000332    4   0.000074
 RIVERF()               0.07     0.00483     0.00091     12   0.012605    0   0.011457    4   0.000190
 OCEANS                 0.00     0.00000     0.00000      0   0.000000    0   0.000000    0   0.000000
 Daily                  3.92     0.26974     0.26974     36   0.682456    0   0.672545    2   0.018732
 Diagnostics            0.07     0.00504     0.00504     12   0.020429    0   0.002189    1   0.001050
PRECIP_LK()             0.05     0.00333     0.00333     24   0.041992    0   0.000079    2   0.000347
SURFACE_LANDICE()       0.00     0.00000     0.00000      0   0.000000    0   0.000000    0   0.000000
-------------------------------------------------------------------------------------------------------
 WARNING: defined in rundeck but not used: subdd =                                                                                                                                 


 ************************************************************************************************************************************
 ************************************************************************************************************************************

  Program terminated due to the following reason:
  >>  Terminated normally (reached maximum time)  <<

 ************************************************************************************************************************************
 ************************************************************************************************************************************
```

If the first hour runs successfully, then you are ready for much longer simulations!

## Perform a ModelE2.1 simulation<a name="run"></a>

Once we've confirmed the model can get through its first hour, we are good to try a longer run, which we refer to as a "production" run.

### 3a. Execute the Run Interactively (Development and Testing)

You can now continue the run from wherever it last ended by using

```console
mpirun -np 6 ./${RUNID}.exe -i I &> ${RUNID}.PRT &
```

Where this is the same as before, except for the `-cold-restart` being removed (which means continue from the last checkpoint file), and the  output is now sent to the `${RUNID}.PRT` file. 

As before, you can check on its contents using tail as the program is running. (You can break out of tail with “Ctrl-c”).

The model will continue until it crashes or completes its simulation, and inform you of its status either way.

### 3b. Stopping and Restarting the Model

To gracefully pause the model while it is running, simply change the contents of the `flagGoStop` text file in the run directory from
```text
___GO___
```
to
```text
__STOP__
```
The model will stop at the next chance it can and save its state for restarting later. 

To restart the model from where it stopped, just re-run the model without specifying `-cold-restart`
```console
mpirun -np 6 ./${RUNID}.exe -i I &> ${RUNID}.PRT &
```

### 3c. Execute a Simulation in a Queue

Coming Soon.

## Process ModelE2.1 output<a name="postprocess"></a>

ModelE2.1 aggregates many diagnostics every month in an “accumulation” file, which you may identify as having “.acc” in the name. A post-processing program called `scaleacc` that you compiled back in [Section 1c](#scaleacc) must be used to convert them to readable NetCDF files.

```console
scaleacc MONYYYY.acc${RUNID}.nc aij	  # Generate 2-D atmospheric data output for Month MON/YYYY
scaleacc MONYYYY.acc${RUNID}.nc aijl  # Generate 3-D atmospheric data output for Month MON/YYYY
```

You may also combine files using the `sumfiles` tool before converting using `scaleacc`. For example,
```console
sumfiles *${YYYY}.acc${RUNID}.nc
```
generates an `ANN${YYYY}.acc${RUNID}.nc` file that contains the annual average values for year `${YYYY}`. And
```console
sumfiles JAN*.acc${RUNID}.nc
```
generates a JANYYYY-YYYY file that contains the multi-year climatology. After aggregation, these acc files may be converted using scaleacc as above.

Once in NetCDF format, you can use whatever software and plotting tool that you prefer for analysis and interpretation (e.g., [R](https://cran.r-project.org/web/packages/ncdf4/index.html), [python](http://unidata.github.io/netcdf4-python/), [Matlab](https://www.mathworks.com/help/matlab/network-common-data-form.html), [NCL](https://www.ncl.ucar.edu/Applications/netcdf4.shtml)).