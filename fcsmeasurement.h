#ifndef FCSMEASUREMENT_H
#define FCSMEASUREMENT_H


//#include <boost/date_time/gregorian/greg_date.hpp>

#include <cmath>
#include <iostream>
#include <string>
#include <stdint.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

//#include <boost/thread/thread.hpp>

#include "browniandynamics.h"
#include "../lib/tools.h"
#include "../lib/jkINIParser2.h"
#include "../lib/multitau-correlator.h"
#include "fluorophordynamics.h"
#include "fluorescencemeasurement.h"
#include "janb/correlator_multitau.h"


/*! \brief Fluorescence Correlation Spectroscopy (FCS) class
    \ingroup diff4_measurement

   This class implements the FCS measurement at a single spot in the simulation volume.
   You can set the position of the laser focus and the detection volume independently, so you
   can experiment with shifts of them. The laser focus and the detection volume are modeled
   as gaussian curves. This is a strong simplification, but it is a simplification that is
   often met in literature. This class also has the possibility to parallelize the calculations
   on a given number of threads.

   Basically this class generates a time series of intensity measurements at a given position.
   Then it may (also on-line) calculate the autocorrelation function of this time series by use of
   a Multiple-Tau correlator (see MultiTauCorrelator calss). The average number of photons in every
   sample timestep \f$ \Delta t_{sample} \f$ (which may be any integer multiple of the simulation
   timestep \f$ \Delta t_{sim} \f$ ) is calculates as:
      \f[ n_{phot}(t_0)=\sum\limits_{t\in\{t_0, t_0+\Delta t_{sim},...,t_0+\Delta t_{sample} \}}\sum\limits_{i=1}^{N_{walker}}\frac{I_0\cdot \sigma_{abs,i}\cdot\Delta t_{sim}}{E_{photon}}\cdot q_{fluor,i}\cdot q_{detector}\cdot f_{spectral}(i,\lambda_{ex})\cdot f_{ex}(x_i,y_i,z_i,\vec{mu}_i)\cdot f_{det}(x_i,y_i,z_i,\vec{mu}_i) \f]
   Where:
      \f[ f_{ex}(x,y,z,\vec{\mu})=(1-F_{ex}+F_{ex}\cdot(\vec{\mu}_i\cdot\vec{p}_{ex})^2\cdot\exp\left(-2\cdot\frac{(x-x_{ex})^2+(y-y_{ex})^2}{\mbox{psf\_r0}^2}-2\cdot\frac{(z-z_{ex})^2}{\mbox{psf\_z0}^2}\right) \f]
      \f[ f_{det}(x,y,z,\vec{\mu})=(\vec{\mu}_i\cdot\vec{p}_{det})^2\cdot\exp\left(-2\cdot\frac{(x-x_{det})^2+(y-y_{det})^2}{\mbox{psf\_r0}^2}-2\cdot\frac{(z-z_{det})^2}{\mbox{psf\_z0}^2}\right) \f]
   Here:
     - \f$ I_0 \f$ is the laser beam intensity peak intensity
     - \f$ \sigma_{abs,i} \f$ is the absorbtion cross section of the i-th walker
     - \f$ \Delta t_{sim} \f$ is the simulation timestep
     - \f$ \Delta t_{sample} \f$ is the sample time of the correlator
     - \f$ E_{photon} \f$ is the energy of an excitation light photon
     - \f$ q_{fluor,i} \f$ is the quantum yeald of the fluorescence of the i-th walker
     - \f$ q_{detector} \f$ is the detection efficiency
     - \f$ (x_i, y_i, z_i)^t \f$ is the position of the i-th walker
     - \f$ (x_{ex}, y_{ex}, z_{ex})^t \f$ is the position of the excitation laser
     - \f$ (x_{det}, y_{det}, z_{det})^t \f$ is the position of the detection volume
     - \f$ \mbox{psf\_r0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the x- and y-direction
     - \f$ \mbox{psf\_z0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the z-direction
     - \f$ \vec{\mu}_i \f$ is the dipole moment direction ( \f$ |\vec{\mu}_i|=1 \f$ ) of the i-th walker
     - \f$ \vec{p}_{ex} \f$ is the (linear) polarisation vector ( \f$ |\vec{p}_{ex}|=1 \f$ ) of the excitation light
     - \f$ F_{pol} \f$ is the fraction of polarized light in the excitation light (0: no polarized light, 1: fully polarized light)
     - \f$ \vec{p}_{det} \f$ is the (linear) polarisation vector ( \f$ |\vec{p}_{det}|=1 \f$ ) of the detection
     - \f$ f_{spectral}(i,\lambda_{ex})\in[0..1] \f$ is a spectral factor that gives the value of the absorbtion spectrum for the i-th
           walker at the given wavelength \f$ \lambda_{ex} \f$. The spectrum is normalized so that its highest peak is set to 1
     - \f$ \lambda_{ex} \f$ is the wavelength of the excitation laser
   .

  The number of photons actually counted is calculated by choosing a random number from a
  poissonian distribution with average value \f$ n_{phot} \f$ . This introduces the photon
  counting statistics of the detector.

 */
class FCSMeasurement: public FluorescenceMeasurement {
    public:
        /** Default constructor */
        FCSMeasurement(FluorophorManager* fluorophors, std::string objectname=std::string(""));
        /** Default destructor */
        virtual ~FCSMeasurement();


        /** \brief initialize the simulation */
        virtual void init();

        /** \brief propagate the detection one step further */
        virtual void propagate();

        /** \brief report the object state */
        virtual std::string report();

        /** \brief save the created data (time series and correlation function) into file
         *         ts001.dat, ts002.dat, corr001.dat, corr002.dat ... as comma-separated values
         */
        virtual void save();

        GetSetMacro(double, psf_r0);
        GetSetMacro(double, psf_z0);
        GetSetMacro(double, lambda_ex);
        GetSetMacro(double, I0);
        GetSetMacro(double, q_det);
        GetSetMacro(double, img_z0);
        GetSetMacro(double, img_x0);
        GetSetMacro(double, img_y0);
        GetSetMacro(double, ex_z0);
        GetSetMacro(double, ex_x0);
        GetSetMacro(double, ex_y0);
        GetSetMacro(double, duration);
        GetSetMacro(double, corr_taumin);
        GetSetMacro(unsigned int, S);
        GetSetMacro(unsigned int, m);
        GetSetMacro(unsigned int, P);
        GetSetMacro(unsigned int, correlator_type);
        GetSetMacro(unsigned long long, timesteps);

    protected:
        /** \brief run the actual FCS simulation */
        void run_fcs_simulation();

        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);
        /** \brief clear all internal data structures */
        void clear();

        /** \brief photon energy [Joule] (calculated)*/
        double Ephoton;

        /** \brief radius (1/e^2 width) of the gaussian PSF in x-y-plane in [micrometers] */
        double psf_r0;

        /** \brief size (1/e^2 width) of the gaussian PSF in z-direction in [micrometers] */
        double psf_z0;

        /** \brief wavelength of excitation light [nanometers] */
        double lambda_ex;

        /** \brief intensity I_0 of excitation laser [uW/m^2] */
        double I0;

        /** \brief quantum efficiency of detection [0..1]*/
        double q_det;

        /** \brief z-position of the laser focus [microns] */
        double img_z0;
        /** \brief x-position of the laser focus [microns] */
        double img_x0;
        /** \brief y-position of the laser focus [microns] */
        double img_y0;

        /** \brief z-position of the detection volume [microns] */
        double ex_z0;
        /** \brief x-position of the detection volume [microns] */
        double ex_x0;
        /** \brief y-position of the detection volume [microns] */
        double ex_y0;

        /** \brief x-component of detector polarisation */
        double d_x;
        /** \brief y-component of detector polarisation */
        double d_y;
        /** \brief z-component of detector polarisation */
        double d_z;

        /** \brief indicates whether to use polarised detection */
        bool polarised_detection;

        /** \brief x-component of excitation polarisation */
        double e_x;
        /** \brief y-component of excitation polarisation */
        double e_y;
        /** \brief z-component of excitation polarisation */
        double e_z;

        /** \brief fraction of polarisation in the excitation light [0..1], 0=unpolarized, 1=completely polarized */
        double e_pol_fraction;

        /** \brief indicates whether to use polarised excitation */
        bool polarised_excitation;

        /** \brief all particle farther than psf_region_factor*psf_r0 or psf_region_factor*psf_z0
         *         away from the laser spot are assumed to not send any fluorescence photons */
        double psf_region_factor;

        /** \brief if this is \c true the object also saves a binned version of the timeseries with
         *         a binning time of save_binning_time */
        bool save_binning;
        /** \brief if save_binning is \c true the object also saves a binned version of the timeseries with
         *         a binning time of save_binning_time */
        double save_binning_time;

        /** \brief if this is \c true the intensity values are shifted into the correlator immediately
         *         and no complete timeseries is built up */
        bool online_correlation;

        /** \brief smallest time interval in correlation function [seconds] */
        double corr_taumin;
        /** \brief number of linear correlators (i.e. number of decades) */
        unsigned int S;
        /** \brief binning ratio between two subsequent linear correlators */
        unsigned int m;
        /** \brief items in each correlator */
        unsigned int P;

        /** \brief this array holds the generated time series */
        uint16_t* timeseries;

        /** \brief this array holds the generated binned time series */
        uint32_t* binned_timeseries;

        /** \brief the number of timesteps recorded in timeseries */
        unsigned long long timesteps;

        /** \brief this array holds the calculated correlation function */
        double* corr;

        /** \brief the times for the corr array */
        double* corr_tau;

        /** \brief number of items in \a corr and \a corr_tau */
        unsigned int slots;

        /** \brief summed runtime of correlation [seconds] */
        double correlation_runtime;

        /** \brief Multiple-Tau correlator class */
        MultiTauCorrelator<uint16_t, uint32_t>* correlator;

        /** \brief JanB's correlator */
        correlatorjb* corrjanb;

        /** \brief end of current integration period */
        double endCurrentStep;

        /** \brief number of integrated photons in the current step */
        double nphot_sum;

        /** \brief counter for the current timestep */
        uint64_t current_timestep;

        /** \brief counts the photons in one display period */
        double display_temp;
        /** \brief used to calculate the binned time series (sum over bins) */
        uint32_t bin_sum;
        /** \brief used to calculate the binned time series (how many bins have been summed up?) */
        uint32_t bin_counter;
        /** \brief index in the binned timeseries  */
        unsigned long long bin_i;
        /** \brief which correlator tu use: 0: JanK's multi-tau,  1:JanB's Multi-tau */
        unsigned int correlator_type;

    private:
};

#endif // FCSMEASUREMENT_H
