#ifndef GRIDRANDOMWALKDYNAMICS_H
#define GRIDRANDOMWALKDYNAMICS_H

#include <cmath>
#include <iostream>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>

#include "../../../LIB/trunk/tools.h"
#include "../../../LIB/trunk/jkiniparser2.h"
#include "fluorophordynamics.h"
#include "diffusiontools.h"


class GridRandomWalkDynamics : public FluorophorDynamics
{
    protected:
        /** \brief diffusion coeffizient in [µm^2/s] */
        double diff_coeff;

        /** \brief width of a grid cell in [micrometers] */
        double grid_constant;

        int64_t gridsize_X;
        int64_t gridsize_Y;
        int64_t gridsize_Z;

        uint8_t* obstacles;
        int64_t obstacle_size;

        double obstacle_fraction;
        int64_t numobstacles;

        double jump_propability;


        /** \brief indicates in which interval of sim timesteps to save the MSD (<=0 => no M;SD is saved) */
        int save_msd_every_n_timesteps;
        int save_msd_factor;
        /** \brief number of MSD values to save */
        int msd_size;

        /** \brief array used to store the msd, if save_msd_every_n_timesteps>0. Every array index \f$ i \f$ corresponds
         *         to a time \f$ i\cdot\mbox{save\_msd\_every\_n\_timesteps}\cdot\mbox{sim\_timestep} \$ . The size of this
         *         array is determined by msd_size */
        double* msd;
        double* msd2;
        double* msdg;
        double* msd2g;
        uint64_t* msd_count;


        /** \brief used to read configuration data from ini file ... overwrite this to read data in derived classes.
         *
         * This base class will care for entering the right group to read from! But keep in mind that you will HAVE to call the
         * corresponding function of the parent class in your implementation!
         */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brief calculate the new jump propability from the given diffusion coefficient, grid_size and timestep */
        void recalcJumpProp();

        /** \brief calculate the new grid size from the given diffusion coefficient and simulation timestep */
        void recalcGridSize();

        void init_obstacles();

    public:
        /** Default constructor */
        GridRandomWalkDynamics(FluorophorManager* fluorophors, std::string object_name=std::string(""));
        /** Default destructor */
        virtual ~GridRandomWalkDynamics();


        /** \brief set the diffusion coefficients for the simulation */
        virtual void set_diff_coeff(double value);

        /** \brief initialize the state of the i-th walker and put it to the given position. The walker step counter is reset to 0 */
        virtual void init_walker(unsigned long i, double x=0, double y=0, double z=0);


        /** \brief set the simulation timestep */
        virtual void set_sim_timestep(double value);

        void set_grid_constant(double gconst);


        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();

        /** \brief set the dimensions of the simulational sphere */
        virtual void set_sim_sphere(double rad);
        /** \brief set the dimensions of the simulational box */
        virtual void set_sim_box(double vx, double vy, double vz);

        /** \brief propagate some walkers [w_start .. w_end] in the simulation one timestep further
         *
         * \param w_start first walker to propagate
         * \param w_end last walker to propagate
         * \param boundary_check when \c boundary_check==false the walkers won't be reset when they reach
         *                       the border of the simulational box
         *
         * This method allows for simultaneous execution of different propagate steps, i.e. a multi-threaded/
         * parallelized simulation. Note that you only have to implement this method if you set multi_threadable
         * to \c true.
         *
         * This method will only propagate the walker \a w_start to \a w_end and leave the others untouched.
         * Be carefull when implementing this and setting multi_threadable to \c true, as the array walker_state
         * is not thread-safe, so maybe you will have to reimplement everything if you want e.g. interaction between
         * several walkers!
         */
        virtual void propagate(bool boundary_check=true);

        /** \brief perform a boundary check for the i-th walker and reset it to a random border position, if it left the sim box */
        virtual void perform_boundary_check(unsigned long i);

        /** \brief report the state of the object */
        virtual std::string report();
    protected:
    private:
};

#endif // GRIDRANDOMWALKDYNAMICS_H
