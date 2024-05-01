# Running the GISS ModelE2.1 GCM to drive the GEOS-Chem CTM 
Lee T. Murray (lee.murray@rochester.edu)
Last Updated: May 1, 2024

## Overview
[ModelE2.1](https://www.giss.nasa.gov/tools/modelE/) is a global general circulation model (GCM) developed primarily at the NASA Goddard Institute for Space Studies in New York, NY. Version "E2.1" of the model is the version that was frozen for the majority of NASA’s contributions to the Coupled Model Intercomparison Project Phase 6 (CMIP6) experiment done in support of the Sixth IPCC Assessment Report [(AR6)](https://www.ipcc.ch/assessment-report/ar6/). 

The emphasis of this tutorial is getting ModelE2.1 to compile and run on your own UNIX-based operating system (MacOS, Linux), and produce the diagnostics necessary to drive the [GEOS-Chem](http://www.geos-chem.org) Chemical Transport Model.

The tutorial was tested using the Ubuntu Desktop 22.04.4 LTS Linux operating system.

If you are new to Linux/UNIX systems, first read the following guides
- [UNIX tutorial for beginners](https://linuxclass.heinz.cmu.edu/doc/Unix-Tutorial-surrey/)
- [How to edit files with emacs](http://www.jesshamrick.com/2012/09/10/absolute-beginners-guide-to-emacs/)

## Table of Contents
1. [First-Time Setup](#first)
2. [Setup a ModelE2.1 Run](#config)
3. [Perform a ModelE2.1 Simulation](#run)
4. [Process ModelE2.1 output](#postprocess)
5. [Archive Diagnostics Necessary for Driving GEOS-Chem](#subdd)

## 1. First-Time Setup<a name="first"></a>

### 1a. Download GISS-GC

This tutorial will assume that the source code is located in the following directory: `$GISS_HOME`

```console
git clone git@github.com:fetch4/GISS-GC.git $GISS_HOME
```

### 1b. Install Software Dependencies

The most important dependency is the [netcdf-fortran](https://github.com/Unidata/netcdf-fortran) library compiled with MPI parallelization. If using a package manager like [miniforge](https://github.com/conda-forge/miniforge) (_recommended for most users_), this will come with all other necessary packages for running ModelE2.1, although you will also want cmake, git and a text editor.

Example environmental yaml files for miniforge are available in `$GISS_HOME/doc/GISS-GC-DOC/envs`, and may be installed using
```console
conda env create -n giss -f giss_ubuntu_22.04.4_lts.yml
```

To load the environment containing all the dependencies necessary to compile and run ModelE2.1
```console
conda activate giss
```

### 1c. Install GISS Post-Processing Tools<a name="scaleacc"></a>

The tools located in `$GISS_HOME/model/mk_diags` are necessary for post-processing the model output into a useful format, particularly the `scaleacc` program. I recommend copying that directory to a common location to which you can point your `$PATH` variable to the binaries.

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

Let’s modify the rundeck we generated (`$GISS_HOME/decks/E6F40_TEST.R`) to run for one month. The row with YEARI indicates when the model starts; the row with YEARE shows when the model will end.  Note: these are placeholder years in this setup as master_yr is set to 1850 earlier in the file, which makes all boundary conditions and climate forcers equal to 1850 climatology.

Setting KDIAG equal to 13*0 reduces the amount of text archived in the log file.

```diff
&INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
-YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
+YEARE=1950,MONTHE=01,DATEE=1,HOURE=0,     KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
```
### 2b. Download Input Files

Now we need to download the input files that are specified in the rundeck file `E6F40_TEST.R`. In practice, it is easiest to try the next step, which will fail, but output a list of the missing files you need.

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

Now you can execute the first hour to make sure everything runs okay. This means run in parallel with 6 processors (`mpirun -np 6`) the program `${RUNID}.exe`, and feed it the input file named `I` (`-i I`), restart from the very beginning of the simulation (`-cold-restart`), save the output into log file (`&> cold_restart.log`), and run in the background (`&`). 
```console
mpirun -np 6 ./${RUNID}.exe -i I -cold-restart &> cold_restart.log &
```
Here we are asking for 6 processors, but you can ask for any number that your system has up to 88 to make the simulation faster.

If all goes well, you should see the model run for the first hour and then gracefully stop. You can follow along as the model outputs `cold_restart.log` using the `tail` command:
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
scaleacc MONYYYY.acc${RUNID}.nc all   # Generate all model output for Month MON/YYYY
```

You may also combine files using the `sumfiles` tool before using `scaleacc`. For example,
```console
sumfiles *${YYYY}.acc${RUNID}.nc
```
generates an `ANN${YYYY}.acc${RUNID}.nc` file that contains the annual average values for year `${YYYY}`. And
```console
sumfiles JAN*.acc${RUNID}.nc
```
generates a JANYYYY-YYYY file that contains the multi-year climatology. After aggregation, these acc files may be converted using scaleacc as above.

Once in NetCDF format, you can use whatever software and plotting tool that you prefer for analysis and interpretation (e.g., [R](https://cran.r-project.org/web/packages/ncdf4/index.html), [python](http://unidata.github.io/netcdf4-python/), [Matlab](https://www.mathworks.com/help/matlab/network-common-data-form.html), [NCL](https://www.ncl.ucar.edu/Applications/netcdf4.shtml)).

## Archive Diagnostics Necessary for Driving GEOS-Chem<a name="subdd"></a>

The default diagnostics for ModelE2.1 are archived as monthly averages. GEOS-Chem needs input meteorology archived at higher temporal resolution (1-hr for 2-D fields, and 3-hr for 3-D fields). Therefore, we need to set up a simulation that includes sub-daily ("subdd") diagnostics. To do so, we are going to use the "GCAP2" rundeck template to create our run directory and compile the model.

```console
cd $GISS_HOME/decks
make rundeck RUN=GCAP2_TEST RUNSRC=GCAP2 OVERWRITE=YES
make -j setup RUN=GCAP2_TEST F90=mpif90
```

The relevant changes are the addition of the following preprocessor statements,
```diff
+#define NEW_IO_SUBDD
+#define CACHED_SUBDD
+#define GCAP
+#CALCULATE LIGHTNING
```
the inclusion of the following modules,
```diff
+lightning
+SUBDD
```
an increased frequency to calling the radiation code,
```diff
-NRAD=5
+NRAD=1
```
and the following diagnostics
```diff
-SUBDD=''
+SUBDD='PS:6i QV:6i T:6i'
+SUBDD1='ALBEDO:2 CLDTOT:2 EFLUX:2 FRSEAICE:2 FRSNO:2 GWETTOP:2'
+SUBDD2='GWETROOT:2 HFLUX:2 LAI:2 LWI:2 PARDF:2 PARDR:2 PBLH:2'
+SUBDD3='PRECANV:2 PRECCON:2 PRECLSC:2 PRECSNO:2 PRECTOT:2 SLP:2'
+SUBDD4='SNODP:2 SNOMAS:2 SWGDN:2 T2M:2 TO3:2 TROPPT:2 TS:2'
+SUBDD5='U10M:2 USTAR:2 V10M:2 Z0M:2 FLASH_DENS:2 CTH:2 QV2M:2'
+SUBDD6='DTRAIN:6 OMEGA:6 RH:6 U:6 V:6'
+SUBDD7='CLOUD:6 OPTDEPTH:6 QI:6 QL:6 TAUCLI:6 TAUCLW:6'
+SUBDD8='DQRCU:6 DQRLSAN:6 REEVAPCN:6 REEVAPLS:6'
+SUBDD9='CMFMC:6 PFICU:6 PFILSAN:6 PFLCU:6 PFLLSAN:6'

+NSUBDD=1         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
+DAYS_PER_FILE=1
```

Then, we will execute the first hour
```console
cd $ModelE_Support/huge_space/$RUNID
/bin/bash ${RUNID}ln
mpirun -np 6 ./${RUNID}.exe -i I -cold-restart &> cold_restart.log &
```
And finally, the full simulation (whether interactively or within a slurm script)
```console
mpirun -np 6 ./${RUNID}.exe -i I &> ${RUNID}.log &
```

In addition to the monthly accumulation (acc) files, the model will now archive daily sub-daily diagnostic (subdd) files named, `${YYYY}${MM}${DD}.subdd${RUNID}.nc`. As with the accumulation files, these must be passed through scaleacc to be converted to a useful form.
```console
scaleacc ${YYYY}${MM}${DD}.subdd${RUNID}.nc all
```
which will generate the following files
```console
${YYYY}${MM}${DD}.aijh2${RUNID}.nc
${YYYY}${MM}${DD}.aijlh6${RUNID}.nc
${YYYY}${MM}${DD}.aijleh6${RUNID}.nc
${YYYY}${MM}${DD}.aijh6i${RUNID}.nc
${YYYY}${MM}${DD}.aijlh6i${RUNID}.nc
```
The nomenclature for naming the files, which each contain diagnostics with common spatial and temporal characteristics are: `a` = atmosphere, `ij` = 2-D fields, `ijl` = 3-D fields at level midpoint, `ijle` = 3-D fields at level edges, `h2` = hourly averages (i.e., 2 time steps), `h6` = 3-hr averages, and `h6i` = 3-hr instantaneous values.

These files are very close to being able to drive GEOS-Chem directly in a GCAP2 simulation. However, the ModelE2.1 daily output files include times from 01:00 to 24:00, and GEOS-Chem requires the input files go from 00:00 to 23:00.

Here is a script that makes use of the extremely useful [cdo](https://code.mpimet.mpg.de/projects/cdo) and [nco](https://nco.sourceforge.net) tools (installed with the conda environment) that can do that conversion.
```bash
#!/bin/bash

conda activate giss

runid=${RUNID}
bd=`pwd`

# Remove any existing tmp files
rm tmp* *.tmp

# Process topographic file
./${runid}ln
if [ ! -d $bd/2x2.5 ]; then mkdir -p $bd/2x2.5; fi
cp TOPO 2x2.5/TOPO
ncap2 -O -s 'focean@units="1"; flake@units="1"; fgrnd@units="1"; fgice@units="1"; zatmo@units="m";' TOPO 2x2.5/TOPO

today=1949-12-01
while [ "$today" != 1950-01-01 ]; do

    # Get today's string
    year=${today:0:4}; month=${today:5:2}; day=${today:8:2}; ymd=${year}${month}${day}

    # Get yesterdays's string
    yest=$(date -I -d "$today - 1 day")
    yyear=${yest:0:4}; ymonth=${yest:5:2}; yday=${yest:8:2}

    # Get tomorrow's string
    tomm=$(date -I -d "$today + 1 day")

    # GISS uses a 365 day calendar
    if [[ $((10#$ymonth)) -eq 02 ]] && [[ $((10#$yday)) -eq 29 ]]; then
	yday=28
    fi
    if [[ $((10#$month)) -eq 02 ]] && [[ $((10#$day)) -eq 29 ]]; then 
	today=$tomm 
	continue
    fi

    if [ ! -f ${year}${month}${day}.subdd${runid}.nc ]; then 
	today=$tomm 
	continue 
    fi

    ########################
    # Special Case on Day 1
    ########################
    if [ "$today" == 1949-12-01 ]; then
	echo $year $month $day
	if [ ! -d $bd/2x2.5/$year/$/month ]; then mkdir -p $bd/2x2.5/$year/$month; fi
	scaleacc ${year}${month}${day}.subdd${runid}.nc all
	for ftype in aijh6i aijlh6i aijh2 aijleh6 aijlh6; do
	    
	    # Fix latitude so N and S pole are correctly 89N and 89S
	    ncap2 -O -s 'lat=array(-89e0,2e0,$lat);' ${year}${month}${day}.${ftype}${runid}.nc ${year}${month}${day}.${ftype}${runid}.nc

	    # Copy and compress the native model resolution 
	    echo 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4

	    # Special treatment for instantaneous files (GEOS includes 0; GISS includes 24 in a given day)
	    if [[ $ftype = "aijh6i" ]] || [[ $ftype = "aijlh6i" ]] ; then
		ncks -O -d time,0   ${year}${month}${day}.${ftype}${runid}.nc tmp1.nc
		ncap2 -O -s 'time=time-3;' tmp1.nc tmp1.nc
		ncks -O -d time,0,6 ${year}${month}${day}.${ftype}${runid}.nc tmp2.nc
		ncrcat -O tmp1.nc tmp2.nc ${year}${month}${day}.${ftype}${runid}.nc
		rm tmp*
	    fi
	    nccopy -k4 -d9 ${year}${month}${day}.${ftype}${runid}.nc 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4 && rm -f ${year}${month}${day}.${ftype}${runid}.nc

	    # Fix the metadata
	    f=2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4
	    f2=${year}${month}${day}.${ftype}${runid}.nc4

	    # Add some optional CF coordinate variables
	    ncatted -O -h -a scale_factor,,c,f,1 $f
	    ncatted -O -h -a add_offset,,c,f,0 $f
	    crDate=`date`
	    case "$ftype" in
		aijh2) title="GISS ModelE2.1 1-hour time-averaged parameters (aijh2), processed for GEOS-Chem input"; dt="010000" ;;
		aijh6i) title="GISS ModelE2.1 3-hour instantaneous parameters (aijh6i), processed for GEOS-Chem input"; dt="030000" ;;
		aijleh6) title="GISS ModelE2.1 3-hour time-averaged parameters on model edges (aijleh6), processed for GEOS-Chem input"; dt="030000" ;;
		aijlh6) title="GISS ModelE2.1 3-hour time-averaged parameters (aijlh6), processed for GEOS-Chem input"; dt="030000" ;;
		aijlh6i) title="GISS ModelE2.1 3-hour instantaneous parameters (aijlh6i), processed for GEOS-Chem input"; dt="030000" ;;
	    esac
	    echo $f $title

	    # Delete and replace global metadata
	    ncatted -O -h -a ,global,d,, $f
	    ncatted -O -h -a Title,global,c,c,"$title" $f
	    ncatted -O -h -a Contact,global,c,c,"Lee T. Murray (lee.murray@rochester.edu)" $f
	    ncatted -O -h -a References,global,c,c,"www.geos-chem.org; wiki.geos-chem.org; http://ees.rochester.edu/atmos/data" $f
	    ncatted -O -h -a Filename,global,c,c,"${f2}" $f
	    ncatted -O -h -a History,global,c,c,"File generated on:  ${crDate}" $f
	    ncatted -O -h -a ProductionDateTime,global,c,c,"File generated on: ${crDate}" $f
	    ncatted -O -h -a ModificationDateTime,global,c,c,"File generated on: ${crDate}" $f
	    ncatted -O -h -a Format,global,c,c,"NetCDF-4" $f
	    ncatted -O -h -a SpatialCoverage,global,c,c,"global" $f
	    ncatted -O -h -a Conventions,global,c,c,"COARDS" $f
	    ncatted -O -h -a Version,global,c,c,"GISS ModelE2.1" $f
	    ncatted -O -h -a VersionID,global,c,c,"$runid" $f
	    ncatted -O -h -a Nlayers,global,c,c,"40" $f
	    ncatted -O -h -a Start_Date,global,c,c,"$today" $f
	    ncatted -O -h -a Start_Time,global,c,c,"00:00:00.0" $f
	    ncatted -O -h -a End_Date,global,c,c,"$today" $f
	    ncatted -O -h -a End_Time,global,c,c,"23:59:59.99999" $f
	    ncatted -O -h -a Delta_Time,global,c,c,"$dt" $f
	    ncatted -O -h -a Delta_Lon,global,c,c,"2.5" $f
	    ncatted -O -h -a Delta_Lat,global,c,c,"2" $f
	done # ftype
    fi

    ########################
    # Normal Cases
    ########################

    if [ ! -f ${yyear}${ymonth}${yday}.subdd${runid}.nc ]; then today=$tomm; continue; fi
    if [ ! -d $bd/2x2.5/$year/$/month ]; then mkdir -p $bd/2x2.5/$year/$month; fi

    echo $year $month $day

    scaleacc ${year}${month}${day}.subdd${runid}.nc all
    scaleacc ${yyear}${ymonth}${yday}.subdd${runid}.nc aijh6i,aijlh6i

    for ftype in aijh6i aijlh6i aijh2 aijleh6 aijlh6; do

	# Fix latitude so N and S pole are correctly 89N and 89S
	ncap2 -O -s 'lat=array(-89e0,2e0,$lat);' ${year}${month}${day}.${ftype}${runid}.nc ${year}${month}${day}.${ftype}${runid}.nc

	# Copy and compress the native model resolution 
	echo 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4
	nccopy -k4 -d9 ${year}${month}${day}.${ftype}${runid}.nc 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4 && rm -f ${year}${month}${day}.${ftype}${runid}.nc

	# Special treatment for instantaneous files (GEOS includes 0; GISS includes 24 in a given day)
	if [[ $ftype = "aijh6i" ]] || [[ $ftype = "aijlh6i" ]] ; then
	    ncap2 -O -s 'lat=array(-89e0,2e0,$lat);' ${yyear}${ymonth}${yday}.${ftype}${runid}.nc ${yyear}${ymonth}${yday}.${ftype}${runid}.nc
	    nccopy -k4 -d9 ${yyear}${ymonth}${yday}.${ftype}${runid}.nc tmp${ymd}.nc4 && rm -f ${yyear}${ymonth}${yday}.${ftype}${runid}.nc
	    ncrcat -O tmp${ymd}.nc4 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4 tmp2${ymd}.nc4 && rm -f tmp${ymd}.nc4
	    cdo -O seldate,${year}-${month}-${day} tmp2${ymd}.nc4 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4 && rm -f tmp2${ymd}.nc4
	fi

	# Fix the metadata
	f=2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4
	f2=${year}${month}${day}.${ftype}${runid}.nc4

	# Add some optional CF coordinate variables                                                                                                                                                                                                                              
	ncatted -O -h -a scale_factor,,c,f,1 $f
	ncatted -O -h -a add_offset,,c,f,0 $f

	crDate=`date`

	case "$ftype" in
	    aijh2) title="GISS ModelE2.1 1-hour time-averaged parameters (aijh2), processed for GEOS-Chem input"; dt="010000" ;;
	    aijh6i) title="GISS ModelE2.1 3-hour instantaneous parameters (aijh6i), processed for GEOS-Chem input"; dt="010000" ;;
	    aijleh6) title="GISS ModelE2.1 3-hour time-averaged parameters on model edges (aijleh6), processed for GEOS-Chem input"; dt="010000" ;;
	    aijlh6) title="GISS ModelE2.1 3-hour time-averaged parameters (aijlh6), processed for GEOS-Chem input"; dt="010000" ;;
	    aijlh6i) title="GISS ModelE2.1 3-hour instantaneous parameters (aijlh6i), processed for GEOS-Chem input"; dt="010000" ;;
	esac
	echo $f $title

	# Delete and replace global metadata                                                                                                                                                                                                                                    
	ncatted -O -h -a ,global,d,, $f
	ncatted -O -h -a Title,global,c,c,"$title" $f
	ncatted -O -h -a Contact,global,c,c,"Your Name (your.email@domain)" $f
	ncatted -O -h -a References,global,c,c,"www.geos-chem.org; wiki.geos-chem.org; http://atmos.earth.rochester.edu/giss-gc" $f
	ncatted -O -h -a Filename,global,c,c,"${f2}" $f
	ncatted -O -h -a History,global,c,c,"File generated on:  ${crDate}" $f
	ncatted -O -h -a ProductionDateTime,global,c,c,"File generated on: ${crDate}" $f
	ncatted -O -h -a ModificationDateTime,global,c,c,"File generated on: ${crDate}" $f
	ncatted -O -h -a Format,global,c,c,"NetCDF-4" $f
	ncatted -O -h -a SpatialCoverage,global,c,c,"global" $f
	ncatted -O -h -a Conventions,global,c,c,"COARDS" $f
	ncatted -O -h -a Version,global,c,c,"GISS ModelE2.1" $f
	ncatted -O -h -a VersionID,global,c,c,"$runid" $f
	ncatted -O -h -a Nlayers,global,c,c,"40" $f
	ncatted -O -h -a Start_Date,global,c,c,"$today" $f
	ncatted -O -h -a Start_Time,global,c,c,"00:00:00.0" $f
	ncatted -O -h -a End_Date,global,c,c,"$today" $f
	ncatted -O -h -a End_Time,global,c,c,"23:59:59.99999" $f
	ncatted -O -h -a Delta_Time,global,c,c,"010000" $f
	ncatted -O -h -a Delta_Lon,global,c,c,"2.5" $f
	ncatted -O -h -a Delta_Lat,global,c,c,"2" $f

    done # ftype

    # Delete subdd file if successfully processed
    if [ -f 2x2.5/$year/$month/${year}${month}${day}.${ftype}${runid}.nc4 ]; then rm -f ${yyear}${ymonth}${yday}.subdd${runid}.nc; fi

    today=$tomm
done # Date

exit 0;
```

The contents of `${ModelE_Support}/huge_space/${RUNID}/2x2.5` may now be used to drive a GEOS-Chem simulation.