/*
    Copyright (c) 2008-2015 Jan W. Krieger (<jan@jkrieger.de>), German Cancer Research Center + IWR, University Heidelberg

    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef FLUORESCENCEMEASUREMENT_H
#define FLUORESCENCEMEASUREMENT_H


#include "jkiniparser2.h"
#include <string>
#include <vector>
#include "tools.h"
#include "fluorophordynamics.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "ticktock.h"

#include "diffusiontools.h"


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
        virtual void save_results();

        /** \brief this function is called once after the simulation has finished.
         */
        virtual void finalize_sim();

        GET_SET_MACRO(double, sim_timestep);
        GET_MACRO(double, sim_time);
        GET_SET_MACRO(std::string, basename);
        GET_SET_MACRO(std::string, object_name);
        GET_SET_MACRO(std::string, description);
        GET_MACRO(std::string, group);
        GET_MACRO(std::string, supergroup);
        GET_SET_MACRO(int, object_number);


        virtual bool depends_on(const FluorescenceMeasurement* other) const;


        /*! \brief return DOT-code (for GraphViz) that represents the node without links (these are returned by dot_get_links() ).
        */
        virtual std::string dot_get_node(bool simple=false) ;
        /*! \brief return DOT label-code (for GraphViz) that represents the nodes properties
        */
        virtual std::string dot_get_properties() ;
        /*! \brief return DOT-code (for GraphViz) that represents the links of the nodes returned by dot_get_node()
        */
        virtual std::string dot_get_links() ;

    protected:

        /** \brief save the results of the measurement
         */
        virtual void save()=0;


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

        /** \brief a description that can e.g. be used in plots */
        std::string description;

        /** \brief number of the object */
        int object_number;
        /** \brief group name of this object */
        std::string group;
        /** \brief supergroup name of this object */
        std::string supergroup;
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
