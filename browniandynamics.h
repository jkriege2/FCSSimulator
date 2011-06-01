#include <cmath>
#include <iostream>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>

#include "../../../LIB/trunk/tools.h"
#include "../../../LIB/trunk/jkINIParser2.h"
#include "fluorophordynamics.h"
#include "diffusiontools.h"

#ifndef BROWNIANDYNAMICS_H
#define BROWNIANDYNAMICS_H

/** \brief number of different diffusion coefficients */
#define DCOUNT 10

/*! \brief simulates independent brownian motion for all walkers
    \ingroup diff4_dynamic

   \section diff4_dynamic_transdiff Translational Diffusion
   Brownian translational dynamics is implemented in a very simple fassion. The diffusion speed is
   determined by a diffusion coefficient \f$ D \f$. From this we can calculate the mean jump width
   for one coordinate direction per simulation timestep \f$ \Delta t_{sim} \f$ as:
     \f[ \sigma_{jump}=\sqrt{2D\cdot\Delta t_{sim}} \f]
   Now we can choose three random numbers \f$ \Delta x, \Delta y, \Delta z \f$ from a gaussian distribution
   centered around 0 and with width \f$ \sigma_{jump} \f$. Those are the increments of the particle coordinates.


   \section diff4_dynamic_rotdiff Rotational Diffusion
   The implementation of the rotational diffusion is somewhat more difficult, but we try to stic to the formalism
   of translational diffusion. First we define a rotational diffusion coefficient \f$ D_r \f$, with:
     \f[ \sigma_{\alpha,jump}=\sqrt{2D_r\cdot\Delta t_{sim}} \f]
   where \f$ \sigma_{\alpha,jump} \f$ is the mean angle between the current direction \f$ \vec{\mu}_n \f$ and the
   next direction \f$ \vec{\mu}_{n+1} \f$ if the fluorophor dipole moment. The random walk process is described
   like a random walk on a sphere with radius 1:
     -# In every timestep we select a random direction for the next step, i.e. an angle \f$ \phi\in[-90^\circ..90^\circ] \f$
     -# We make a step of length \f$ \Delta\alpha \f$ randomly choosen from a \f$ N(0, \sigma_{\alpha,jump}) \f$
        distribution in the direction of step 1.
   .
   To do so, we first construct a perpendicular basis (2 vectors \f$ \vec{n}', \vec{n}'' \f$ ) of the plane perpendicular
   to the current \f$ \vec{\mu}_n \f$. Then we choose the angle \f$ \phi\in[-90^\circ..90^\circ] \f$ and get a vector
     \f[ \Delta\vec{\mu}=\cos(\phi)\cdot\vec{n}'+\sin(\phi)\cdot\vec{n}'' \f]
   which has a random direction perpendicular to \f$ \vec{\mu}_n \f$.
   Now the two vectors \f$ \vec{\mu}_{n} \f$ and \f$ \Delta\vec{\mu} \f$ define the plane in which the new \f$ \vec{\mu}_{n+1} \f$
   lies. Sow we can get \f$ \vec{\mu}_{n+1} \f$ as a rotation of \f$ \vec{\mu}_{n} \f$ araound an axis perpendicular to this plane,
   where we rotate for a randomly choosen vector:
     \f[ \vec{\mu}_{n+1}=\mathbf{R(\theta, \vec{\mu}_{n}\times\Delta\vec{\mu})}\;\vec{\mu}_{n} \f]
   From <a href="http://de.wikipedia.org/wiki/Drehmatrix">http://de.wikipedia.org/wiki/Drehmatrix</a> we get this formulation for
   a rotation matrix around an arbitrary axis (unit vector) \f$ \vec{v}=(v_1,v_2,v_3)^t \f$:
     \f[ \mathbf{R}(\alpha, \vec{v})=\begin{pmatrix} \cos \alpha +v_1^2 \left(1-\cos \alpha\right)   & v_1 v_2 \left(1-\cos \alpha\right) - v_3 \sin \alpha &  v_1 v_3 \left(1-\cos \alpha\right) + v_2 \sin \alpha \\   v_2 v_1 \left(1-\cos \alpha\right) + v_3 \sin \alpha  & \cos \alpha + v_2^2\left(1-\cos \alpha\right) &   v_2 v_3 \left(1-\cos \alpha\right) - v_1 \sin \alpha         \\ v_3 v_1 \left(1-\cos \alpha\right) - v_2 \sin \alpha &  v_3 v_2 \left(1-\cos \alpha\right) + v_1 \sin \alpha & \cos \alpha + v_3^2\left(1-\cos \alpha\right)\end{pmatrix}  \f]



   \section diff4_dynamic_bigparticles Multi-Fluorophore Particles
   It is also possible to simulate the motion of particles which consist of many fluorophores. To do so each particle is
   said to consist of exactly \f$ N_{fl} \f$ fluorophores which are positioned on its periphery. In the simples case such
   a particle is a sphere with radius \f$ r \f$ .
 */

class BrownianDynamics: public FluorophorDynamics
{
    protected:
        /** \brief rotational diffusion coefficient [rad^2/s] */
        double diff_rot;
        /** \brief diffusion coeffizient in [µm^2/s] */
        double diff_coeff[DCOUNT];
        /** \brief a second diffusion coeffizient in [µm^2/s] */
        //double diff_coeff1;

        /** \brief there are two areas of diffusion coefficients. One is in the x<diffarea_x0 half space (diff_coeff) and
         *         one in the x>=diffarea_x0 half space (diff_coeff1) */
        double diffarea_x0[DCOUNT];

        /** \brief width of jump length distribution (standard deviation of gaussian!) in [µm] for diff_coeff */
        double sigma_jump[DCOUNT];
        /** \brief width of jump length distribution (standard deviation of gaussian!) in [µm] for diff_coeff1 */
        //double sigma_jump1;
        /** \brief width of a single angle jump for rotational diffusion */
        double sigma_rotjump;
        /** \brief if set (\c true ) this class also uses rotational diffusion */
        bool use_rotational_diffusion;
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
        uint64_t* msd_count;

        /** \brief number of fluorophores per random walker/particle */
        unsigned int n_fluorophores;

        /** \brief used to read configuration data from ini file ... overwrite this to read data in derived classes.
         *
         * This base class will care for entering the right group to read from! But keep in mind that you will HAVE to call the
         * corresponding function of the parent class in your implementation!
         */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brief get the number of walkers for the current simulation settings */
        virtual unsigned long calc_walker_count();

    public:
        /** \brief class constructor with standard volume 30*30*30µm^3 and a concentration of 1nM */
        BrownianDynamics(FluorophorManager* fluorophors, std::string object_name=std::string(""));

        /** \brief class constructor, initialises with a box of the given dimensions */
        BrownianDynamics(FluorophorManager* fluorophors, double sim_x, double sim_y, double sim_z, double c_fluor, std::string object_name=std::string(""));

        /** \brief class constructor, initialises with a sphere of the given dimensions */
        BrownianDynamics(FluorophorManager* fluorophors, double sim_radius, double c_fluor, std::string object_name=std::string(""));

        /** \brief class destructor */
        virtual ~BrownianDynamics();

        /** \brief set the diffusion coefficients for the simulation */
        inline virtual void set_diff_coeff(double value, double value1) {
            diff_coeff[0]=value;
            diff_coeff[1]=value1;
            for (int i=0; i<DCOUNT; i++) {
                sigma_jump[i]=sqrt(2.0*diff_coeff[i]*sim_timestep);
            }
            sigma_rotjump=sqrt(2.0*diff_rot*sim_timestep);
        };

        /** \brief set the diffusion coefficients for the simulation */
        inline virtual void set_diff_coeffi(int n, double value) {
            diff_coeff[n]=value;
            for (int i=0; i<DCOUNT; i++) {
                sigma_jump[i]=sqrt(2.0*diff_coeff[i]*sim_timestep);
            }
            sigma_rotjump=sqrt(2.0*diff_rot*sim_timestep);
        };

        /** \brief set the rotational diffusion coefficients for the simulation */
        inline virtual void set_diffrot_coeff(double value) {
            diff_rot=value;
            sigma_rotjump=sqrt(2.0*diff_rot*sim_timestep);
        };

        /** \brief set the simulation timestep */
        inline virtual void set_sim_timestep(double value) {
            FluorophorDynamics::set_sim_timestep(value);
            for (int i=0; i<DCOUNT; i++) {
                sigma_jump[i]=sqrt(2.0*diff_coeff[i]*sim_timestep);
            }
            sigma_rotjump=sqrt(2.0*diff_rot*sim_timestep);
        };


        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();


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


        /** \brief report the state of the object */
        virtual std::string report();

        /** \brief tests the brownian dynamics by computing some walker steps and recording <Delta X^2>
         *
         * Generate a file \a filename which contains three data columns:
         *  -# time [microseconds]
         *  -# DeltaX^2
         *  -# 6*D*t
         *  -# DeltaAlpha^2
         *  -# 2*Drot*t
         *
         * If you plot the data of column 2 and 3 against column 1 the curves should be the same. This functions initialises
         * 100 walkers at 0,0,0 and propagates them for \a steps steps. So call init() afterwards!
         *
         * This also outputs the configuration information and a GnuPlot .plt file to plot the results.
         */
        virtual void test(unsigned int steps=1000, unsigned int walkers=100);

        /*GetMacro(double, diff_coeff);
        GetMacro(double, diff_coeff1);*/
        inline double get_diff_coeff() { return diff_coeff[0]; };
        inline double get_diff_coeff1() { return diff_coeff[1]; };
        GetMacro(double, diff_rot);
        GetSetMacro(bool,use_rotational_diffusion);

};

#endif // BROWNIANDYNAMICS_H
