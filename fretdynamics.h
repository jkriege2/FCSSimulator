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


#ifndef FRETDYNAMICS_H
#define FRETDYNAMICS_H

#include "fluorophordynamics.h"

/*! \brief this dynamics class takes the DYNAMIC trajectories (so it ignores the derived trajectories)
           of an input class and forwards them to the output, possibly adding and simulating additional
           fluorophores
    \ingroup diff4_dynamic


 */
class FRETDynamics : public FluorophorDynamics
{
    public:
        /** Default constructor */
        FRETDynamics(FluorophorManager* fluorophors, std::string object_name);
        /** Default destructor */
        virtual ~FRETDynamics();

        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();

        /** \brief output the next step from the trajectory
         *
         * \param boundary_check when \c boundary_check==false the walkers won't be reset when they reach
         *                       the border of the simulational box
         */
        virtual void propagate(bool boundary_check=true);

        /** \brief report the state of the object, as human-readable text */
        virtual std::string report();

    protected:
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);


        /** \brief name of the parent class to read from */
        std::string parent;
        /** \brief returns the parent class to read from */
        FluorophorDynamics* get_parent() const;

        /*! \brief returns \c true , if this object depends on the given \ other object.

            This function is used to determine the order in which the propagate() method of all dynamics objects in this
            simulation is called.
         */
        virtual bool depends_on(const FluorophorDynamics* other) const;

        virtual void handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores);

};

#endif // FRETDYNAMICS_H
