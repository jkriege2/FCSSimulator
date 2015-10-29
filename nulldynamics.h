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

#ifndef NULLDYNAMICS_H
#define NULLDYNAMICS_H


#include "fluorophordynamics.h"

/*! \brief this dynamics class does not output any trajectories
    \ingroup diff4_dynamic


 */
class NullDynamics : public FluorophorDynamics
{
    public:
        /** Default constructor */
        NullDynamics(FluorophorManager* fluorophors, std::string object_name);
        /** Default destructor */
        virtual ~NullDynamics();

        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();

        /** \brief output the next step from the trajectory
         *
         * \param boundary_check when \c boundary_check==false the walkers won't be reset when they reach
         *                       the border of the simulational box
         */
        virtual void propagate(bool boundary_check=true);




    protected:
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brife indicates whether there is a single (non-moving) walker in this dynamics object */
        bool has_walker;

        double x;
        double y;
        double z;


};


#endif // NULLDYNAMICS_H
