# GalevNB--------------------------------------
 ReadMe of software package "GalevNB"
--------------------------------------

 X. Pang, C. Olczak, Q. Shu, J. Li

 First edit: 26 February 2013
 Last  edit: 2 October   2022

--------------------------------------


1. Introduction:

   GalevNB is a software package based on the original GALEV code (Kotulla et al. 2009) developed to deliver (integrated or individual) observational
   magnitudes and spectra of modelled stellar systems using individual stellar masses, temperatures, luminosities, and metallicities.  The code is
   tuned to process data output from the state-of-the-art direct N-body integrator NBODY6 / NBODY6++ (Spurzem 1999, Nitadori & Aarseth 2012).



2. Motivation:

   The generation of (individual or integrated) observational magnitudes and spectra of model stellar systems is the key for a *direct* comparison of
   numerical simulations and observations that serves three fundamental scientific goals:
   a) determination of the initial conditions (insight into past evolution)
   b) understanding of present-day features and processes
   c) prediction of future evolution



3. Structure:

   GalevNB has two subfolders:
    ------galevnb_linux_bev
    ------galevnb_linux_sevsin

      Each subfolder contains a complete set of programs, and you can view these two programs as two different versions of galevnb:
      Program in galevnb_linux_bev applies to binaries
      Program in galevnb_linux_sevsin applies to single stars

      For both programs in galevnb_linux_bev and galevnb_linux_sevsin:
        
            1) The main program GalevNB.f90 parses snapshot files generated by NBODY6 / NBODY6++ and uses seven routines of the GALEV package to convert
                effective temperatures, stellar luminosities, metallicities, and masses into observational magnitudes and spectra.

            2) The GalevNB package contains four subfolders:
                a) "spectral_templates": All the spectral template files from the BaSeL library of model atmospheres (Lejeune, Cuisinier & Buser 1997, 1998).

                b) "standard_filters": A large set of filter response curves (FUV,NUV,U,B,V,R,I,J,H,K) that are used as standard filters for reference.

                c) "filter_response_curves": Filter response curves of the five magnitude systems HST, 2MASS, SDSS, Johnson, and Cousins are provided in
                    separate subfolders, intended for a user-specific choice of filter sets for generated magnitudes.

                Information about the entire set of filters available is included in the file "filterlist.dat". The file contains seven columns:
                - col 1: filter name
                - col 2: corresponding file path
                - col 3: id of selected zero point; please change this id number according to the selected zero point id!!!Very important!!
                - col 4: standard zero point in the Vega magnitude system (id=1)
                - col 5: standard zero point in the AB magnitude system   (id=2)
                - col 6: standard zero point in the ST magnitude system   (id=3)
                - col 7: optional user-defined zero point                 (id=4)

                    ******************************************************************************************************************
                    * NOTE: The file "filterlist.dat" must be present in the same directory as the NBODY6 / NBODY6++ snapshot files. *
                    ******************************************************************************************************************
            


4. Installation:  

        1) Unpack the GalevNB archive in your favorite location, subsequently called <GalevNB_installation_path>.

        Since the structure of both programs is the same, here we use galevnb_linux_bev as an example (The steps to install the program in galevnb_linux_bev
        are the same):

            2) cd <GalevNB_installation_path>/galevnb_linux_sevsin

            3) edit zmet in file GalevNB.f90 (line 230):
                    zmet = n 
                    [Here n = \log_{}{\frac{metallicity}{solar~abundance} }, metallicity is the the metallicity of target cluster]
                    For examplem, if the metallicity of target cluster equals to solar abundance, set zmet=0
                    If the metallicity of target cluster equals to 0.01*solar abundance, set zmet=-2
                    
            4) Type "make" to compile.
                If you need to recompile the program, just type "make" again.
                In case you need more debugging output during run time compile via "make debug".
                For cleaning the directory type "make clean".

        5) Add the GalevNB installation path to your personal shell configuration file (".bashrc",".tcshrc",...) in your home directory.
            
            a) for bash:
                - add this line to your .bashrc: "export PATH=${PATH}:<GalevNB_installation_path>galevnb_linux_bev:<GalevNB_installation_path>galevnb_linux_sevsin"
                -> NOTE: this will take effect only in a new terminal instance
                        for immediate effect in your current terminal type "source .bashrc"

            b) for tcsh:
                - add this line to your .tcshrc: "setenv PATH ${PATH}:<GalevNB_installation_path>galevnb_linux_bev:<GalevNB_installation_path>galevnb_linux_sevsin"
                -> NOTE: this will take effect only in a new terminal instance
                        for immediate effect in your current terminal type "source .tcshrc"



5. Execution:

   1) Make sure you have a set of snapshot files that were generated from a NBODY6 / NBODY6++ simulation.

   2) Copy the file "<GalevNB_installation_path>/galevnb_linux_sevsin(or galevnb_linux_bev)/filter_response_curves/filterlist.dat" to the location of the snapshot files.

   3) Edit the file "filterlist.dat": uncomment (i.e. remove the leading "#" of) your preferred filters.
      NOTE: The current maximum are 20 filters.

   4) Type "GalevNBsevs(or "GalevNBbev" for binaries)" to process the snapshot files.
      The following output files will be created:

      a) "stellar_magnitudes-*": snapshot files containing the stellar properties and observational magnitudes of individual stars for the chosen filters.
   
      b) "cluster_int_mag.dat": list of integrated magnitudes of the stellar system for each snapshot and the chosen filters.

      c) "spec_out-*": snapshots of the integrated spectrum of the stellar system covering a wavelength range of 90..10^6 Angstrom.

      [TODO:  d) Individual stellar spectra: this file can be produced upon user request.    !!(Christoph, we may add this choice for the user in the main file?) ]



6. Input parameters:
* Structure of the slibrary vector ******************************
 *
 * metall_bin:                  < -2.0    index 0                  Solar metallicity is [Fe/H]=0.0; 
 *               -2.0 <= [Fe/H] < -1.0          1
 *               -1.0 <=        <  0.0          2
 *                0.0 <=                        3
 *
 * ---------------------------------------------------------------
 *                
 * teff:         2000 <= Teff < 3000            0
 *               3000 <= Teff < 4000            1
 *               .
 *               .
 *               .
 *              49000 <= Teff < 50000           47
 *              50000 <= Teff ------------> Blackbody
 *
 * ----------------------------------------------------------------
 * 
 * logg                  logg < -2.0            0
 *               -2.0 <= logg < -1.0            1
 *               .
 *               .
 *                5.0 <= logg < 6.0             8
 *                6.0 <= logg                   9
 * 
 ***************************************************************** */
E.g., in the spectra file (under directory:spectral_templates/Lejeune), 
file name:09000_p2.0_m1.5_p00.spec means 9000K effective temperature, log(g)=+2.0, [Fe/H]=-1.5, [alpha/H]=+0.0 ... p/m are for plus/minus



7. Bug fixes:
   - C.O., 2013-08-13: Fixing two memory leaks (spectra.h, spectra_blank.cpp).