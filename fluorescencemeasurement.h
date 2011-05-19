#include "../../../LIB/trunk/jkINIParser2.h"
#include <string>
#include <vector>
#include "../../../LIB/trunk/tools.h"
#include "fluorophordynamics.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../../../LIB/trunk/ticktock.h"

#include "diffusiontools.h"

#ifndef FLUORESCENCEMEASUREMENT_H
#define FLUORESCENCEMEASUREMENT_H




/*! \brief abstract base class for all fluorescence measurement schemes
   \ingroup diff4_measurement

  This class provides the basic interface for measurement classes. It conatins a random number
  generator where the type is read from an ini file and a link to a dynamics class.

  You can overwrite the methods as you need it:
   - read_config() allows to read configuration parameters from an INI file
   - init() initializes the simulation
   - propagate() works through all the currently supplied fluorophor states
   - addWalkerState() adds a pointer to a FluorophorDynamics::walkerState array that should be
                      used for measurement
   - save() saves the result of the simulation into files. All filenames start with
            \a basename, so you can use it to store them in a specified directory and with
            a given prefix.
  .

  The simulation works as follows:
  A main program controls the  overall simulation, i.e. it manages a set of dynamics generation classes
  and calls their propagate() methods to generate a new set of fluorophor states. The measurement classes
  contain a list of dynamics classes that provide data to them. The main program subsequently calls the
  measurement classes's propagate() function which then read the current states from the dynamics classes
  and update their output.

 */
class FluorescenceMeasurement: public TickTock
{
    public:
        /** Default constructor */
        FluorescenceMeasurement(FluorophorManager* fluorophors, std::string objectname=std::string(""));

        /** Default destructor */
        virtual ~FluorescenceMeasurement();

        /** \brief read configuration from INI file
         *
         * \param parser the parser object to read the data from
         * \param group The configuration data is read from the given group in the  parser
         * \param supergroup If a supergroup is provided this method first readsthe configuration from the super group and afterwards from the group.
         *
         * The \a supergroup parameter allows to use a stacked configration scheme: A super group contains parameters that are used for all simulation
         * objects of the same type while the \a group may be used to further specialize these parameters for every single object.
         */
        virtual void read_config(jkINIParser2& parser, std::string group=std::string("measurement"), std::string supergroup=std::string(""));

        /** \brief initialize the simulation */
        virtual void init();

        /** \brief propagate the detection one step further */
        virtual void propagate();

        /** \brief clear the currently saved references to walker states */
        inline void clearDynamics() { dyn.clear(); };

        /** \brief add a new reference to a walker states */
        inline void addDynamics(FluorophorDynamics* d) { dyn.push_back(d); };

        /** \brief report the object state as human-readable string */
        virtual std::string report();

        /** \brief save the results of the measurement
         */
        virtual void save()=0;

        GetSetMacro(double, sim_timestep);
        GetMacro(double, sim_time);
        GetSetMacro(std::string, basename);
        GetSetMacro(std::string, object_name);

    protected:
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brief GSL helper object: random number generator type */
        const gsl_rng_type * rng_type;

        /** \brief GSL helper object: random number generator */
        gsl_rng * rng;

        /** \brief the simulation time step [seconds] */
        double sim_timestep;

        /** \brief the overall total simulation time [seconds] */
        double sim_time;

        /** \brief gives the duration of the measurement */
        double duration;

        /** \brief basename for the output files */
        std::string basename;

        /** \brief an object name */
        std::string object_name;

        /** \brief a pointer to a vector of FluorophorDynamics::walkerState arrays.
         *
         * The vector actually contains tuples where the first element contains the number of
         * state arrays and the second is a pointer to the first array item.
         */
        std::vector<FluorophorDynamics*> dyn;
        /** \brief pointer to the fluorophor database */
        FluorophorManager* fluorophors;
};

#endif // FLUORESCENCEMEASUREMENT_H
