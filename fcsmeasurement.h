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
#include <gsl/gsl_sf_trig.h>

//#include <boost/thread/thread.hpp>

#include "browniandynamics.h"
#include "../../../LIB/trunk/tools.h"
#include "../../../LIB/trunk/jkiniparser2.h"
#include "../../../LIB/trunk/multitau-correlator.h"
#include "fluorophordynamics.h"
#include "fluorescencemeasurement.h"
#include "../../../LIB/trunk/correlator_multitau.h"
#include "../../../LIB/trunk/statistics_tools.h"
#include "../../../LIB/trunk/jkimage.h"
#include "../../../LIB/trunk/tinytiffwriter.h"

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
     - \f$ \mbox{detpsf\_r0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the x- and y-direction
     - \f$ \mbox{detpsf\_z0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the z-direction
     - \f$ \mbox{expsf\_r0} \f$ is the \f$ 1/e^2 \f$-width of the excitation profile in the x- and y-direction
     - \f$ \mbox{expsf\_z0} \f$ is the \f$ 1/e^2 \f$-width of the excitation profile in the z-direction
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

  The background photons are added to the photons from the fluorescence as a second poissonian process, so
    \f[ n_{phot+back}=\mathcal{P}_{Poisson}(n_{phot})+\mathcal{P}_{Poisson}(r_{back}\cdot\Delta t_{sample}) \f]
  where \f$ r_{back} \f$ is the background count rate in photons/second.

  In addition a number of offset photons may be added with a gaussian distribution of definable variance:
    \f[ n_{phot+back+offset}= n_{phot+back}+\mathcal{P}_{Gaussian}(n_{offset}, \sigma_{offset}^2) \f]


  The detection may be limited to a certain number of photons per simulation timestep.

  It is possible to select different PDPs for the detection and illumination. If gaussian functions are selected,
  the with is ALWAYS given as 1/e^2-width, so the gaussian takes the form
    \f[  g(x)=\exp(-2\cdot x^2/\sigma^2) \f]

  These are the illumintaion modes are documented \link ill_distribution here.\endlink

  The detection PSF modes are documented \link det_distribution here.\endlink

  It is possible to simulate different detectors with different properties:
    - <b>photon counting detectors:</b> Here the number of detected photons is simply determined by drawing a random number from a poissonian
      distribution (mean=variance!!!) with the given men photon number.
    - <b>linear detectors (camera FCS):</b> The average detected ADU value is calculated as
        \f[ n_{lin}(t)=\mbox{lindet\_gain}\cdot  n_{phot}(t) \f]
      Then the actual detected value is drawn from a gaussian distribution with average \f$ n_{lin}(t) \f$  and variance
      \f$ n_{lin}(t)\cdot\text{lindet\_var\_factor}+\text{lindet\_readnoise}^2 \f$. Then the value is cast to an integer (i.e. quantisation by an analog to digital converter).
      The range of these values is then limited to \f$ 0..2^{\text{lindet\_bits}}-1 \f$ to account for the finite resolution of the ADC.
  .

  Some detector PSFs are given as a radial symetric image \f$ \mbox{PSF}(r,z) \f$ only. The simulation is able to use such an image and
  simulate the detection efficiency a square pixel of length \f$ a \f$ (parameter \a pixel_size ) would have in an optical system with
  such a PSF. To do so it calculates this integral numerically:
    \f[ I_{det}(x_o, y_o, z_o)=\left[\int\limits_{-a/2}^{a/2}\int\limits_{-a/2}^{a/2} f(\sqrt{(x-x_o)^2+(y-y_o)^2}, z_o)\;\mathrm{d}x\;\mathrm{d}y\right]/\max\limits_{z}\left[\iint\limits_{-\infty}^\infty f(\sqrt{x^2+y^2}, z)\;\mathrm{d}x\;\mathrm{d}y\right] \f]
  Here a fluorescent point source is positioned at \f$ (x_o, y_o, z_o) \f$ in object space and the camera pixel is positioned symmetrically around
  \f$ (0,0,0) \f$ . If the function \f$ f(r,z) \f$ is not known analytically or very complex to calculate, it may be stored as an image
  (\a psf_rz_image_detection) with \a psf_rz_image_rwidth * \a psf_rz_image_zwidth pixels size and a resolution of \a psf_rz_image_rresolution
  micrometers in radial and \a psf_rz_image_zresolution micrometers in z-direction. Note that the image is symmetric with respect to z, so the
  real coordinates span the space \code -psf_rz_image_zresolution*(psf_rz_image_zwidth/2) ... psf_rz_image_zresolution*(psf_rz_image_zwidth/2)  \endcode
  Values between pixels are estimated using bilinear interpolation.
  The integral above is calculated analytically using a simple rectangualr grid on the pixel of size \a pixel_width (\f$ a \f$ ) micrometer. The grid-size is
  \a pixel_size_integrationdelta micrometers, but may be reduced so the pixel is subdivided into at least 10*10 grid points (see square_integrate() ).
  Then the approximation is:
    \f[ N_{grid}=\max\left(10, \frac{\mbox{pixel\_width}}{\mbox{pixel\_size\_integrationdelta}}\right) \f]
    \f[ \Delta a=\frac{\mbox{pixel\_width}}{N_{grid}} \f]
    \f[ \tilde{I}_{det}(x_o, y_o, z_o)=\left[\sum\limits_{x=0}^{N_{grid}-1}\sum\limits_{y=0}^{N_{grid}-1} f(\sqrt{(x-N_{grid}/2-x_o)^2+(y-N_{grid}/2-y_o)^2}, z_o)\cdot \Delta a^2\right]/\max\limits_{z}\left[\sum\limits_{r=0}^{r_{max}} f(r, z)\cdot \left(\pi(r+1)^2-\pi r^2\right)\cdot\mbox{psf\_rz\_image\_zresolution}^2 \right] \f]

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

        virtual void finalize_sim();

        GetSetMacro(double, expsf_r0);
        GetSetMacro(double, expsf_z0);
        GetSetMacro(double, detpsf_r0);
        GetSetMacro(double, detpsf_z0);
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
        GetSetMacro(unsigned int, detector_type);
        GetSetMacro(unsigned int, lindet_bits);
        GetSetMacro(double, lindet_gain);
        GetSetMacro(double, lindet_readnoise);
        GetSetMacro(double, lindet_var_factor);
    protected:

        /** \brief save the created data (time series and correlation function) into file
         *         ts001.dat, ts002.dat, corr001.dat, corr002.dat ... as comma-separated values
         */
        virtual void save();

        /** \brief run the actual FCS simulation */
        void run_fcs_simulation();

        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);
        /** \brief clear all internal data structures */
        void clear();

        /** \brief photon energy [Joule] (calculated)*/
        double Ephoton;

        /** \brief radius (1/e^2 width) of the gaussian excitation PSF in x-y-plane in [micrometers] */
        double expsf_r0;

        /** \brief size (1/e^2 width) of the gaussian excitation PSF in z-direction in [micrometers] */
        double expsf_z0;

        /** \brief radius (1/e^2 width) of the gaussian detection PSF in x-y-plane in [micrometers] */
        double detpsf_r0;

        /** \brief size (1/e^2 width) of the gaussian detection PSF in z-direction in [micrometers] */
        double detpsf_z0;

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

        /** \brief if this is \c true the object saves the timeseries  */
        bool save_timeseries;
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
        int32_t* timeseries;
        uint64_t timeseries_size;

        /** \brief this array holds the generated binned time series */
        int32_t* binned_timeseries;
        uint64_t binned_timeseries_size;

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
        MultiTauCorrelator<double, double>* correlator;

        /** \brief JanB's correlator */
        correlatorjb<double, double>* corrjanb;

        /** \brief end of current integration period */
        double endCurrentStep;

        /** \brief number of integrated photons in the current step */
        double nphot_sum;

        /** \brief counter for the current timestep */
        uint64_t current_timestep;

        /** \brief counts the photons in one display period */
        double display_temp;
        /** \brief used to calculate the binned time series (sum over bins) */
        int32_t bin_sum;
        /** \brief used to calculate the binned time series (how many bins have been summed up?) */
        uint32_t bin_counter;
        /** \brief index in the binned timeseries  */
        unsigned long long bin_i;
        /** \brief which correlator tu use: 0: JanK's multi-tau,  1:JanB's Multi-tau */
        unsigned int correlator_type;

        /** \brief this is the maximum number of photons that may be detected/counted during one measurement timestep. Before the correlation/storing,
         *         the number of detected photons is limited to this. */
        int32_t max_photons;

        /** \brief this is the minimum number of photons that may be detected/counted during one measurement timestep. Before the correlation/storing,
         *         the number of detected photons is limited to this. */
        int32_t min_photons;

        /** \brief minimum wavelength of the detection filter (switched off if <0) */
        double det_wavelength_min;
        /** \brief maximum wavelength of the detection filter (switched off if <0) */
        double det_wavelength_max;
        /** \brief background photon count rate in 1/second=Hz */
        double background_rate;
        /** \brief offset  photon count rate in [photons/(detection step)] */
        double offset_rate;
        /** \brief offset  photon counts standard deviation [photons] */
        double offset_std;
        /** \brief remove this offset photon counts [photons/(detection step)] from the acquired signal before correlation */
        double offset_correction;

        /** \brief plot PSF from 0 to psfplot_xmax in x direction in micrometers */
        double psfplot_xmax;
        /** \brief plot PSF from 0 to psfplot_ymax in y direction in micrometers */
        double psfplot_ymax;
        /** \brief plot PSF from 0 to psfplot_zmax in z direction in micrometers */
        double psfplot_zmax;

        /** \brief a list of other fcs objects of which the results shall be plotted in one Gnuplot file */
        std::vector<std::string> plot_with;
        std::map<int, std::vector<std::string> > plot_with_more;



        /*! \brief illumination intensity distribution

            possible values:
              - 0: gaussian \f[ I(x,y,z)=\exp\left(-2\cdot\frac{x^2+y^2}{\mbox{expsf\_r0}^2}-2\cdot\frac{z^2}{\mbox{expsf\_z0}^2}\right) \f]
              - 1: gaussian_SPIM (lateral: uniform, longitudinal: gaussian) \f[ I(x,y,z)=\exp\left(-2\cdot\frac{z^2}{\mbox{expsf\_z0}^2}\right) \f]
              - 2: slit_SPIM (lateral: uniform, longitudinal: \f[ I(x,y,z)=\left(\frac{\sin(\pi\cdot z/\mbox{expsf\_z0})}{\pi\cdot z/\mbox{expsf\_z0}}\right)^2 \f]
              - 3: gaussian_beam: \f[ I(x,y,z)=\left(\frac{\mbox{expsf\_r0}}{W(z, \mbox{expsf\_z0}, \mbox{expsf\_r0})}\right)^2\cdot\exp\left(-2\cdot\frac{x^2+y^2}{W(z, \mbox{expsf\_z0}, \mbox{expsf\_r0})^2}\right) \ \ \ \text{with}\ \ \  W(z, w_0, z_0)=w_0\cdot\sqrt{1+\frac{z^2}{z_0^2}} \f]
         */
         int ill_distribution;

         /*! \brief type of the detector

             possible values:
               - 0: photon counting detector
               - 1: digitized linear detector
          */
         unsigned int detector_type;

         /*! \brief for linear detectors: factor between detector noise variance and detected intensity

             example values for Andor iXon X3 860:
               <code>lindet_var_factor = EMGain/10</code>
          */
         double lindet_var_factor;

         /** \brief resolution of linear detector digitization (range is <code>0..2^lindet_bits-1</code>) */
         unsigned int lindet_bits;

         /** \brief gain of linear detector  */
         double lindet_gain;

         /** \brief readout noise standard deviation of linear detectors */
         double lindet_readnoise;

         /** \brief r-z-image to sample from for complex PSF of detection, this array represents a function \f$ f(r,z) \f$ */
         JKImage<double> psf_rz_image_detection;
         /** \brief this contains the maximal integral over each plane represented in psf_rz_image_detection, i.e. \f$ A_{max}=\max\limits_z\int\limits_0^{2\pi}\int\limits_0^\infty f(r,z)\cdot r\;\mathrm{d}r\;\mathrm{d}\varphi \f$ */
         double psf_rz_image_detection_integral_max;

         /** \brief r-z-image to sample from for complex PSF of illumination */
         JKImage<double> psf_rz_image_illumination;

        /** \brief estimates the contents of psf_rz_image_detection_integral
         *
         * This function calculates \f[ \max\limits_{z}\left[\iint\limits_{-\infty}^\infty f(\sqrt{x^2+y^2}, z)\;\mathrm{d}x\;\mathrm{d}y\right]=\max\limits_{z}\left[\sum\limits_{r=0}^{r_{max}} f(r, z)\cdot \left(\pi(r+1)^2-\pi r^2\right)\cdot\mbox{psf_rz_image_rresolution}^2 \right] \f]
         * The result is stored in psf_rz_image_detection_integral_max and the radial PSF \f$ f(r,z) \f$ is taken from psf_rz_image_detection.
         */
         void estimate_psf_integrals();

         /** \brief r-pixels in psf_rz_image_* */
         uint32_t psf_rz_image_rwidth;
         /** \brief z-pixels in psf_rz_image_* */
         uint32_t psf_rz_image_zwidth;

         /** \brief r-resolution (microns/pixel) in psf_rz_image_* */
         double psf_rz_image_rresolution;
         /** \brief z-resolution (microns/pixel) in psf_rz_image_* */
         double psf_rz_image_zresolution;

         /*! \brief integrate a radial function defined in image over a square patch (width/height a, resolved with into patches of size da*da or smaller)

             The pixel is centered around (dx,dy,dz). The radial function \f$ f(r,z) \f$ is encoded by \a image (a r-z-image, with resolution \a imaghe_dr and \ image_dz).
             The integration is performed by simply summing over the pixel's values. The radial function is linearly interpolated between two consecutive points.

                \f[ N_{grid}=\max\left(10, \frac{a}{\mbox{da}}\right) \f]
                \f[ \Delta a=\frac{a}{N_{grid}} \f]
                \f[ \tilde{I}_{det}(\mbox{dx}, \mbox{dy}, \mbox{dz})=\left[\sum\limits_{x=0}^{N_{grid}-1}\sum\limits_{y=0}^{N_{grid}-1} f(\sqrt{\Delta a^2\cdot(x-N_{grid}/2-\mbox{dx}/\Delta a+0.5)^2+\Delta a^2\cdot(y-N_{grid}/2-\mbox{dy}/\Delta a+0.5)^2}, \mbox{dz})\cdot \Delta a^2\right]/\max\limits_{z}\left[\sum\limits_{r=0}^{r_{max}} f(r, z)\cdot \left(\pi(r+1)^2-\pi r^2\right)\cdot\mbox{image\_dr}^2 \right] \f]

             The normalization constant \f$ \max\limits_{z}\left[\sum\limits_{r=0}^{r_{max}} f(r, z)\cdot \left(\pi(r+1)^2-\pi r^2\right)\cdot\mbox{image\_dr}^2 \right] \f$
             is not evaluated in every call to this function, but only once in estimate_psf_integrals() which is called towards the end of init().

             \note The integration square is divided into at least 10*10 grid points to guarantee a good integration accuracy!

         */
         double square_integrate(double dx, double dy, double dz, double a, double da, const JKImage<double>& image, double image_dr, double image_dz) const;


         std::string ill_distribution_to_str(int i) const;
         int str_to_ill_distribution(std::string i) const;

         /*! \brief detection propability distribution

            possible values:
              - 0: gaussian \f[ \mbox{PSF}(x,y,z)=\exp\left(-2\cdot\frac{x^2+y^2}{\mbox{detpsf\_r0}^2}-2\cdot\frac{z^2}{\mbox{detpsf\_z0}^2}\right) \f]
              - 1: square pixel (uses additional parameter pixel_size) \f[ \mbox{PSF}(x,y,z)=\exp\left(-2\cdot\frac{z^2}{\mbox{detpsf\_z0}^2}\right)\cdot  \frac{\mbox{erf}\left[\frac{\mbox{pixel\_size}-2x}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]+\mbox{erf}\left[\frac{\mbox{pixel\_size}+2x}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]}{2\cdot\mbox{erf}\left[\frac{\mbox{pixel\_size}}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]}\cdot\frac{\mbox{erf}\left[\frac{\mbox{pixel\_size}-2y}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]+\mbox{erf}\left[\frac{\mbox{pixel\_size}+2y}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]}{2\cdot\mbox{erf}\left[\frac{\mbox{pixel\_size}}{\sqrt{2}\cdot\mbox{detpsf\_r0}}\right]} \f]
              - 2: gaussian beam \f[ \mbox{PSF}(x,y,z)=\left(\frac{\mbox{detpsf\_r0}}{W(z, \mbox{detpsf\_z0}, \mbox{detpsf\_r0})}\right)^2\cdot\exp\left(-2\cdot\frac{x^2+y^2}{W(z, \mbox{detpsf\_z0}, \mbox{detpsf\_r0})^2}\right) \ \ \ \text{with}\ \ \  W(z, w_0, z_0)=w_0\cdot\sqrt{1+\frac{z^2}{z_0^2}} \f]
              - 3: gaussian beam pixel: the same as 2 ("gaussian beam") but with a square pixel (uses additional parameter pixel_size) \f[ \mbox{PSF}(x,y,z)=\frac{\iint\limits_{-\mbox{pixel_size}/2}^{+\mbox{pixel_size}/2}f(\sqrt{x^2+y^2}, z)\;\mathrm{d}x\;\mathrm{d}y}{\max\limits_z\iint\limits_{-\infty}^\infty f(\sqrt{x^2+y^2}, z)\;\mathrm{d}x\;\mathrm{d}y} \f] with \f[ f(r,z)=\left(\frac{\mbox{detpsf\_r0}}{W(z, \mbox{detpsf\_z0}, \mbox{detpsf\_r0})}\right)^2\cdot\exp\left(-2\cdot\frac{r^2}{W(z, \mbox{detpsf\_z0}, \mbox{detpsf\_r0})^2}\right) \f] Here the 2D-integral over the pixel is calculated as described in the main documentation of this class
         */
         int det_distribution;
         std::string det_distribution_to_str(int i) const;
         int str_to_det_distribution(std::string i) const;

         /** \brief size of possibly used detection pixels in [microns] */
         double pixel_size;
         /*! \brief size of an integration patch on the pixel for some detection PSF models [microns]

             If the program has to integrate over a pixel of width/height pixel_size, it actually calculates a sum over the
             value of the PSF taken on a grid with width pixel_size_integrationdelta on the pixel.
         */
         double pixel_size_integrationdelta;

         double detectionEfficiency(double dx, double dy, double dz) const;
         double illuminationEfficiency(double dx, double dy, double dz) const;
         int32_t getDetectedPhotons(double nphot_sum) const;

         /** \brief maximum avg. photon counts for detector test */
         double ndettest_max;
         /** \brief avg. photon counts step size for detector test */
         double ndettest_step;

         time_t start_time;

    private:
};

#endif // FCSMEASUREMENT_H
