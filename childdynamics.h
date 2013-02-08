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

#endif // CHILDDYNAMICS_H
