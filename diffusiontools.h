#ifndef DIFFUSIONTOOLS_H
#define DIFFUSIONTOOLS_H

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <exception>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <cstdlib>
#include <stdexcept>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

#include "../../../LIB/trunk/tools.h"
#include "../../../LIB/trunk/jkiniparser2.h"
#include "../../../LIB/trunk/datatable.h"



/** \brief error handling: exceptions for FluorophorDynamics class (and descendents)
 *  \ingroup diff4_dynamic
 */
class FluorophorException : public std::exception {
  private:
       std::string errormessage;
	public:
		/** \brief class constructors */
		FluorophorException() throw() { errormessage="unknown error in FluorophorException"; };

		/** \brief class destructor */
		virtual ~FluorophorException() throw() {};

		/** \brief constructor with supplied error message */
		FluorophorException(std::string msg) throw() { errormessage=msg; };

		/** \brief  returns the assigned errormessage */
		std::string getMessage() const { return errormessage; };

		/** \brief returns the error description as C string */
		virtual const char* what() const throw() { return getMessage().c_str(); };
};



/** \brief error handling: exceptions for FluorophorDynamics class (and descendents)
 *  \ingroup diff4_measurement
 */
class MeasurementException : public std::exception {
  private:
       std::string errormessage;
	public:
		/** \brief class constructors */
		MeasurementException() throw() { errormessage="unknown error in FluorophorMeasurement"; };

		/** \brief class destructor */
		virtual ~MeasurementException() throw() {};

		/** \brief constructor with supplied error message */
		MeasurementException(std::string msg) throw() { errormessage=msg; };

		/** \brief  returns the assigned errormessage */
		std::string getMessage() const { return errormessage; };

		/** \brief returns the error description as C string */
		virtual const char* what() const throw() { return getMessage().c_str(); };
};


/** \brief this class is used to manage fluorophor properties and spectra
 *  \ingroup diff4_tools
 *
 * There is also a system in this class that allows to lookup a spectral prefactor [0..1] in a known set of spectra.
 * The known spectra are loaded in the constructor and more may be added by calling add_spectrum(). The lookup
 * is implemented in get_spectral_efficiency(). It is sufficient to give some points of the spectrum, as
 * get_spectral_efficiency() interpolates between them. Spectra should be given in \c *.spec files that contain
 * two columns of data: The first is the wavelength in [nm] and the second is the absorption efficiency [0..1]
 * where 1 is maximum absorptiona and 0 is no absorption. The highest value in the spectrum should be 1.
 */

class FluorophorManager {
    public:
        /** \brief this struct contains all data neede to describe a spectrum */
        struct Spectrum {
            /** \brief the file the spectrum was loaded from */
            std::string filename_abs;
            /** \brief the wavelengths from the file [nanometers] */
            double* lambda;
            /** \brief the absorption efficiencies from the file [0..1] */
            double* eff_abs;
            /** \brief the fluorescence from the file [0..1] */
            double* eff_fl;
            /** \brief number of data points in lambda and eff */
            int val_count;
            /** \brief GSL lookup accelerator for absorbance */
            gsl_interp_accel *accel_abs;
            /** \brief this stores the spline calculated from the data in the file for absorbance */
            gsl_spline *spline_abs;
            /** \brief GSL lookup accelerator for fluorescence */
            gsl_interp_accel *accel_fl;
            /** \brief this stores the spline calculated from the data in the file for fluorescence */
            gsl_spline *spline_fl;
            /** \brief indicates whether this spectrum has been loaded and interpolated */
            bool loaded;
        };

        /** \brief specifies the properties of one fluorophor */
        struct FluorophorData {
            /** \brief the fluorescence efficiency of the fluorophor [0..1] */
            double fluorescence_efficiency;
            /** \brief fluorescence lifetime [sec] */
            double fluorescence_lifetime;
            /** \brief absorption cross section [meters^2] */
            double sigma_abs;
            /** \brief bleaching propability [0..1] */
            double bleaching_propability;
            /** \brief triplet lifetime [sec] */
            double triplet_lifetime;
            /** \brief triplet transition propability [0..1] */
            double triplet_propability;
            /** \brief spectrum ID of fluorophor */
            int spectrum;
        };

    protected:
        /** \brief specifies the type of interpolation for the spectra */
        const gsl_interp_type *spectral_interpolation_type;

        /** \brief if this is set \c true the program will output a set of values from each spectrum into a \c *.ispec file */
        bool test_spectra;

        /** \brief this vector stores all available spectra. It maps from the filename (without extension) to the Spectrum ID */
        std::map<std::string, int> spectra_map;

        /** \brief this vector contains all available spectra */
        std::vector<Spectrum> spectra;

        /** \brief load all \c *.spec file from the current directory */
        void init_spectra();

        /** \brief fluorophor database */
        std::map<std::string, FluorophorData> fluorophor_database;

        /** \brief parse the \c fluorophors.ini file and fill fluorophor_database
         *
         * The \c fluorophors.ini file contains a list of the fluorescence parameters of the different fluorophors.
         * Each entry is characterized by a name. Note that if you also have a spectrum for the dye the name of the
         * spectrum file has to be \c name.spec where \c name is the same name as in \c fluorophors.ini.
         *
         * Each entry has this form:
         * \verbatim
[Rho6G]
q_fluor=0.94
tau_fluor=3.9e-9
sigma_abs=3.6e-20 \endverbatim
         * Here \ Rho6G is the name and the spectrum is saved as \c rho6g.spec (case is ignored!). The property q_fluor
         * give the fluorescence efficiency (i.e. \f$ \frac{\mbox{fluorescence photons}}{\mbox{absorbed photons}} \f$) in
         * the range [0..1], \c tau_fluor gives the fluorescence lifetime in seconds and \c sigma_abs (\f$ \sigma_{abs} \f$)
         * the absorption cross  section in \f$ \mbox{meters}^2 \f$. Instead of \c sigma_abs you can also provide
         * c molar_extinction which is the molar extinction coefficient \f$ \epsilon_{mol} \f$ in 1/(cm*Molar)=liters/(mol*cm).
         * These two properties are connected by \f[ \sigma_{abs}=\frac{\epsilon_{abs}}{N_A}=\frac{\epsilon_{abs}[l/(\mbox{cm}\cdot\mbox{mol})]\cdot 10^{-1}}{6.022\cdot10^{23}}[\mbox{meters}^2] \f]
         */
        void init_fluorophor_database();

        /** \brief path of the fluorophor database */
        std::string database_path;

    private:

    public:
        /** \brief class constructor, initialises the database from the provided path. */
        FluorophorManager(std::string database_path=std::string(""));

        /** \brief class destructor */
        ~FluorophorManager();

        /** \brief load the specified spectrum, i.e. make sure that the specral data is loaded and the splines are calculated */
        void load_spectrum(int ID);

        /** \brief get the absorption efficiency [0..1] from the given spectrum at the given wavelength [nm]
         *
         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
         * from the given spectrum.
         */
        inline double get_spectral_absorbance(std::string spectrum, double wavelength) {
            return get_spectral_absorbance(spectra_map[spectrum], wavelength);
        }


        /** \brief get the absorption efficiency [0..1] from the given spectrum at the given wavelength [nm]
         *
         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
         * from the given spectrum.
         */
        inline double get_spectral_absorbance(int spectrum, double wavelength) {
            if (spectrum==-1) return 1.0;
            if (!spectra[spectrum].loaded) load_spectrum(spectrum);
            gsl_spline* spline=spectra[spectrum].spline_abs;
            gsl_interp_accel* acc=spectra[spectrum].accel_abs;
            return GSL_MIN(1.0, GSL_MAX(0.0, gsl_spline_eval(spline, wavelength, acc)));
        }

        /** \brief get the fluorescence value [0..1] from the given spectrum at the given wavelength [nm]
         *
         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
         * from the given spectrum.
         */
        inline double get_spectral_fluorescence(std::string spectrum, double wavelength) {
            return get_spectral_fluorescence(spectra_map[spectrum], wavelength);
        }

        /** \brief get the fluorescence value [0..1] from the given spectrum at the given wavelength [nm]
         *
         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
         * from the given spectrum.
         */
        inline double get_spectral_fluorescence(int spectrum, double wavelength) {
            if (spectrum==-1) return 1.0;
            if (!spectra[spectrum].loaded) load_spectrum(spectrum);
            gsl_spline* spline=spectra[spectrum].spline_fl;
            gsl_interp_accel* acc=spectra[spectrum].accel_fl;
            return GSL_MIN(1.0, GSL_MAX(0.0, gsl_spline_eval(spline, wavelength, acc)));
        }

        /** \brief get the fluorescence value [0..1] from the given spectrum at the given wavelength [nm]
         *
         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
         * from the given spectrum.
         */
        double get_spectral_fluorescence(int spectrum, double wavelength_start, double wavelength_end);

        /** \brief add a spectrum file to the internal database, but do NOT load the file
         *         loading is done internally by calling load_spectrum() */
        void add_spectrum(std::string filename);

        /** \brief returns the number of spectra */
        inline size_t getSpectraCount() { return spectra.size(); };

        /** \brief returns the number of fluorophor data */
        inline size_t getFluorophorCount() { return fluorophor_database.size(); };

        /** \brief return fluorophor dataset */
        inline FluorophorData& getFluorophorData(std::string f) {
            return fluorophor_database[f];
        }

        /** \brief returns \ true if fluorophor exists in database */
        inline bool fluorophorExists(std::string f) {
            //std::cout<<"trying to find '"<<f<<"'\ndatabase size "<<fluorophor_database.size()<<"\n";
            return (fluorophor_database.find(f)!=fluorophor_database.end());
        }

        /** \brief return thefilename of the specifiued spectrum */
        inline std::string getSpectrumFilename(int spectrum) {
            if (spectrum<0) return "none";
            return spectra[spectrum].filename_abs;
        }
};



#endif // DIFFUSIONTOOLS_H
