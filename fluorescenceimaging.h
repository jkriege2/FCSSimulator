
#ifndef FLUORESCENCEIMAGING_H
#define FLUORESCENCEIMAGING_H

#include <cmath>
#include <iostream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include "browniandynamics.h"
#include "tools.h"
#include "jkiniparser2.h"
#include "fluorescencemeasurement.h"

/*! \brief create images of the fluorescence
    \ingroup diff4_measurement

   This class creates a series of images of the fluorescence in the simulational volume.
   To do so, it creates a rasterized plane (parallel to x-y-plane) where the fluorescence
   is recordet during the exposure time (summing up the photons) in every pixel. Between two
   subsequent images this class may propagate the system for a time imaging_delay.

   At every pixel the class counts those photons that are inside a detection volume with
   a gaussian shape (x-y-width \c psf_r0, z-width \c psf_z0 ). All pixels are measured in
   parallel. This class also takes into account the detection statistics by really counting
   a number of photons that is randomly choosen from a poissonian distribution with \f$ n_{phot} \f$
   as the mean value/parameter. The quantity \f$ n_{phot} \f$ is calculated as follows:
    \f[ n_{phot}(x_0,y_0,z_0)=\sum\limits_{i=1}^{N_{walker}}\frac{I_0\cdot \sigma_{abs}\cdot\Delta t_{sim}}{E_{photon}}\cdot q_{fluor}\cdot q_{detector}\cdot\exp\left(-2\cdot\frac{(x_i-x_0)^2+(y_i-y_0)^2}{\mbox{psf\_r0}^2}-2\cdot\frac{(z_i-z_0)^2}{\mbox{psf\_z0}^2}\right) \f]
   Here:
     - \f$ I_0 \f$ is the laser beam intensity (flat illumination!)
     - \f$ \sigma_{abs} \f$ is the absorbtion cross section
     - \f$ \Delta t_{sim} \f$ is the simulation timestep
     - \f$ E_{photon} \f$ is the energy of an excitation light photon
     - \f$ q_{fluor} \f$ is the quantum yeald of the fluorescence
     - \f$ q_{detector} \f$ is the detection efficiency
     - \f$ (x_i, y_i, z_i)^t \f$ is the position of the i-th walker
     - \f$ (x_0, y_0, z_0)^t \f$ is the position of the current pixel
     - \f$ \mbox{psf\_r0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the x- and y-direction
     - \f$ \mbox{psf\_z0} \f$ is the \f$ 1/e^2 \f$-width of the detection profile in the z-direction
   .

   The resulting series of images is saved as comma-separated values files called \c img001.date,
   \c img002.dat ... They may be read by Matlab into matrices for evaluation without further conversion.
 */
class FluorescenceImaging: public FluorescenceMeasurement
{
    private:
    protected:

        /** \brief photon energy [Joule] (calculated)*/
        double Ephoton;
        /** \brief radius of the gaussian PSF in x-y-plane in [micrometers] */
        double psf_r0;
        /** \brief size of the gaussian PSF in z-direction in [micrometers] */
        double psf_z0;
        /** \brief exposure time for imaging [seconds] */
        double exposure_time;
        /** \brief time between two conscutive images [seconds] */
        double imaging_delay;
        /** \brief wavelength of excitation light [nanometers] */
        double lambda_ex;
        /** \brief intensity I_0 of excitation laser [uW/m^2] */
        double I0;
        /** \brief absorbption cross section [m^2] */
        double sigma_abs;
        /** \brief quantum efficiency of fluorescence [0..1]*/
        double q_fluor;
        /** \brief quantum efficiency of detection [0..1]*/
        double q_det;
        /** \brief number of the current image */
        unsigned int current_image;
        /** \brief z-position of the imaging plane [microns] */
        double img_z0;
        /** \brief x start position of the imaging plane [microns] */
        double img_x0;
        /** \brief y start position of tha imaging plane [microns] */
        double img_y0;
        /** \brief image width in x direction [microns] */
        double img_dx;
        /** \brief image width in y direction [microns] */
        double img_dy;
        /** \brief pixel count in x direction */
        unsigned int img_xpixel;
        /** \brief pixel count in y direction */
        unsigned int img_ypixel;
        /** \brief simulation time when current state is completed */
        double timeStepComplete;
        /** \brief if \c false we are in delay between two images, if \c true we are currently acquiring an image */
        bool acquiringImage;

        /** \brief if set \c true the images will be saved in binary format (array of doubles/ints) */
        bool save_binary_images;

        /** \brief the current average image */
        gsl_matrix* avgImage;
        /** \brief the current count image */
        gsl_matrix* cntImage;

        /** \brief saves the current images to files */
        void save_current_images();

        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brief save the created matrices into files img001.dat, img002.dat ... as comma-separated values
         *
         * \param basename this string is prepended to all files that are generated by this routine. So if \c basename=".\results\r0"
         *                 the program will generate files \c .\results\r0config.ini, \c .\results\r0config.txt ...
         */
         //* \param inifilename filename (including path) to the INI file that configures the simulation. If this is present save() may store
         //*                    a copy of it together with the output.
        virtual void save();
    public:
        /** \brief class constructor */
        FluorescenceImaging(FluorophorManager* fluorophors, std::string objectname=std::string(""));
        /** \brief class destructor */
        virtual ~FluorescenceImaging();

        GET_SET_MACRO(double, psf_r0);
        GET_SET_MACRO(double, psf_z0);
        GET_SET_MACRO(double, exposure_time);
        GET_SET_MACRO(double, lambda_ex);
        GET_SET_MACRO(double, I0);
        //GET_SET_MACRO(double, sigma_abs);
        //GET_SET_MACRO(double, q_fluor);
        GET_SET_MACRO(double, q_det);
        GET_SET_MACRO(double, img_z0);
        GET_SET_MACRO(double, img_x0);
        GET_SET_MACRO(double, img_y0);
        GET_SET_MACRO(double, img_dx);
        GET_SET_MACRO(double, img_dy);
        GET_SET_MACRO(unsigned int, img_xpixel);
        GET_SET_MACRO(unsigned int, img_ypixel);
        GET_MACRO(unsigned int, current_image);


        /** \brief initialize the simulation */
        virtual void init() {
            init(exposure_time, imaging_delay);
        };

        /** \brief initialize the simulation to generate \a n_images images with an exposure time of \a exposure_time seconds
         *         each and a delay (with full brownian motion simulation) of \a imaging_delay seconds between two consecutive
         *         images.
         */
        virtual void init(double exposure_time, double imaging_delay);

        /** \brief run the full simulation after init() */
        virtual void propagate();

        /** \brief report the object state */
        virtual std::string report();

};

#endif // FLUORESCENCEIMAGING_H
