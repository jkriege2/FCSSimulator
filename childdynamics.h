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


#ifndef CHILDDYNAMICS_H
#define CHILDDYNAMICS_H

#include "fluorophordynamics.h"

/*! \brief this dynamics class takes the DYNAMIC trajectories (so it ignores the derived trajectories)
           of an input class and forwards them to the output, possibly adding and simulating additional
           fluorophores
    \ingroup diff4_dynamic


 */
class ChildDynamics : public FluorophorDynamics
{
    public:
        /** Default constructor */
        ChildDynamics(FluorophorManager* fluorophors, std::string object_name);
        /** Default destructor */
        virtual ~ChildDynamics();

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



        /*! \brief return the number of visible walkers in the simulational box

            \note Use this function to access the visible walkers, visible to a fluorescence detection. These may be less
                  than actually present. Use these methods together: get_visible_walker_count(), get_visible_walker_sigma_times_qfl(),
                  get_visible_walker_state(). And use them only in detection mode. In dynamics mode: use get_walker_count(),
                  get_walker_sigma_times_qfl(), get_walker_state()
         */
        virtual unsigned long get_visible_walker_count() ;
        virtual double get_visible_walker_sigma_times_qfl(unsigned long i);
        /*! \brief get pointer to array with all walker states

            \note Use this function to access the visible walkers, visible to a fluorescence detection. These may be less
                  than actually present. Use these methods together: get_visible_walker_count(), get_visible_walker_sigma_times_qfl(),
                  get_visible_walker_state(). And use them only in detection mode. In dynamics mode: use get_walker_count(),
                  get_walker_sigma_times_qfl(), get_walker_state()
        */
        virtual walkerState* get_visible_walker_state();

    protected:
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);


        /** \brief name of the parent class to read from */
        std::string parent;
        /** \brief returns the parent class to read from */
        FluorophorDynamics* get_parent() const;

        bool dont_copy_photophysics;

        /** \brief if set \c true, also the initial walker is visible */
        bool initial_walker_visible;

        /** \brief if set \c true [default], the exists state of the walker is copied from the parent */
        bool copy_existstate;
        /** \brief if set \c true [default], the QM state of the walker is reset at the boundary */
        bool reset_qmstate_on_boundary;

        /*! \brief returns \c true , if this object depends on the given \ other object.

            This function is used to determine the order in which the propagate() method of all dynamics objects in this
            simulation is called.
         */
        virtual bool depends_on(const FluorophorDynamics* other) const;

        virtual void handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores);
        virtual std::string dot_get_properties() ;
                /*! \brief return DOT-code (for GraphViz) that represents the links of the nodes returned by dot_get_node()
        */
        virtual std::string dot_get_links() ;
};

#endif // CHILDDYNAMICS_H
