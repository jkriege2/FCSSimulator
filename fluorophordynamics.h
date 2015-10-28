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




#ifndef FLUOROPHORDYNAMICS_H
#define FLUOROPHORDYNAMICS_H

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <exception>
#include <cstring>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <cstdlib>

#include "tools.h"
#include "jkiniparser2.h"
#include "diffusiontools.h"



#define N_FLUORESCENT_STATES 6


class RelativeAbsorbanceReader {
    public:
        virtual ~RelativeAbsorbanceReader() {};
        /** \brief returns the relative absorbance of the given walker in the give fluorophore dynamics. */
        double get_relative_absorbance_for(FluorophorDynamics* dyn, int i);
        virtual double get_relative_absorbance_for(FluorophorDynamics* dyn, int i, double x, double y, double z)=0;

};



/*! \brief this is the virtual base class for any class describing fluorophor dynamics
   \ingroup diff4_dynamic

  Basically this class provides an interface that allows to access a set of fluorescent particles which
  are described (see FluorophorDynamics::walkerState) by a position within the simulational volume and
  their internal (quantum-mechanical and photo-physical) properties (absorption cross section, fluorescence
  efficiency, QM state, ...).

  If you want to implement a dynamic you will have to write a class that publicly inherits FluorophorDynamics
  and mainly overwrites the propagate() and init() methods. If you don't wan to implement special features that
  go beyond the things implemented here, you won't have to overwrite other methods. You can signal other classes
  that your dnymaics implementation may be used in a multi-threaded way, by setting the protected data memeber
  multi_threadable to \c true in your constructur. In this case you will have to provide both varieties of the
  propagate() method, so another class may call \c propagate(w_start, w_end) with different parameters in parallel!

  If you need additional parameters from an ini file you will also have to reimplement read_config().

  \section FD_rng Random Number Generators
  This base class already instaciates a GSL random number generator (taus2 as default) in the data member rng
  which you may use for your implementation. By setting the rng property in the ini file you can select one of
  these generators (for details, see the GSL documentation):
  <center><tt>mt19937, ranlxs0, ranlxs1, ranlxs2, ranlxd1, ranlxd2, ranlux, ranlux389, cmrg, mrg, taus, taus2, gfsr4, rand, bsd, libc5, glibc2, rand48, ranf, ranmar, r250, tt800, minstd, knuthran2, knuthran2002, knuthran</tt></center>

  \section FD_init Initialisation of Simulation Environment
  By giving the volume \f$ V \f$ of the simulational box and the fluorophor concentration \f$ c_f \f$ the method
  change_walker_count() will initialize enough memory for \f$ \mbox{ceil}(v\cdot c_f) \f$ particles, set
  then at random positions in the volume and init their state with the given \c init_... values.


  \section FD_photophysics Photophysics Simulation
  This class implements basic photophysics, if activated with use_photophysics.
  The photophysics is implemented in propagate_photophysics() which you will have to call for every walker when implementing
  the dynamics.  If you want to implement your own photophysics simulation, overwrite this method in derved classes!

  The photophysics simulation implements only a stump that supportrs at most N_FLUORESCENT_STATES fluorescent states \f$ i \f$ ,
  with each a fluorescence efficiency \f$ \phi_i\geq0 \f$. The transition propabilities are defined in terms of a matrix
  for the propability to go from state \f$ i \f$ to state \f$ f \f$ : \f$ p_{if}=\mathbb{P}(i\rightarrow f) \f$. The current state
  is stored in qm_state. The function get_walker_qfluor() returns the fluorescence efficiency of the current walker in its current
  state.

  The photophysics may also depend on the fluorescence intensity emitted by a FluorescenceMeasurement object implementing the interface
  RelativeAbsorbanceReader. If this feature is switched on (photophysics_absorbance_dependent=true), the relative absorbance
  of the object abs_reader (if !=NULL) is read for every fluorophore and then is multiplied by photophysics_absorbance_factor.
  This factor \f$\alpha\f$ is used to scale (down) the jump probabilities \f$ p_{if} \forall i\neq f \f$ in the photophysics transition matrix!

  \section FD_depletion Reservoir Depletion Simulation
  This class also implements a simulation of reservoir depletion, where a random-walker is switched off (and left that way) with a certain propbability depletion_probability,
  every time it hits the simulation box boundary.


  \section FD_threads threading support

  By setting \a use_two_walkerstates to \c true the class uses two arrays of walker_state's which may be interchanged by calling
  swap_walker_states(). In this case the method propagate() writes to walker_state and all read functions use walker_stat_other
  so it is possible to fill walker_state while the old walker_state_other is beeing read. If \a use_two_walkerstates is set to \c false
  the two variables simply point to the same memory location.

  \section FD_nfluorophores Additional Fluorophores per Walker
  This feature allows each walker to have a set of  n_fluorophores fluorophores. The first fluorophore is at the center of the
  particle. Then n_fluorophores-1 additional walkers are updated every step: Their position is the central position + a (potentially
  random) displacement. The displacements are either fixed during the whole simulation or may be calculated by special rules in
  each setp (e.g. an oscillator). Each of the n_fluorophores fluorophores may have it's own photophysics.

  For the programmer this means that the fluorescence meausrement classes see more walkers than are actually simulated by the
  dynamics classes. The number of walkers that have to be simulated in terms of dynamics can be read using get_walker_count().
  The number of walkers that have to be taken into account for detection is returned by get_visible_walker_count().

  The organization in walker_state is like this:
  \verbatim
                                        get_visible_walker_count()
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      get_walker_count()
    **********************

    +----------------------------------------------------------------------------------------------------+
    | w1 | w2 | ... | wN || w12 | w22 | ... | wN2 | w13 | w23 | ... | wN3 |  ... | w1M | w2M | ... | wNM |
    +----------------------------------------------------------------------------------------------------+
      N: number of walkers
      M: number of fluorophores per walker
  \endverbatim


  There are two sets of function to access the walkers:
     - use get_visible_walker_count(), get_visible_walker_sigma_times_qfl(), get_visible_walker_state() to access the visible walkers in detection
     - use use get_walker_count(), get_walker_sigma_times_qfl(), get_walker_state() to acees the walkers when simulating dynamics!
  .
  Distinguishing this allows e.g. to have the above mentioned additional walker, but the original walkers invisble, if
  get_visible_walker_state() returns a pointer to w12 instead of w1 (see diagram above) and get_visible_walker_count()
  returns N*(M-1) instead of N*M. Note: By default this is NOT implemented and both sets of functions return the SAME VALUES!
  See ChildDynamics for an example of how this could work!

 */
class FluorophorDynamics
{
    public:
        /** \brief state of one random walker
         *
         * Each walker is described by an instance of this struct. In addition to the
         * position of the walker this struct contains members that describe the internal
         * state of the walker:
         *  - \c qm_state is an integer which may be used to store the (quantum-mechanic) state of the walker
         *             this can be used to also simulate photophysical phenomena, like triplet blinking ...
         *  - \c sigma_abs is the absorbtion cross section of the current walker
         *  - \c q_fluor is the quantum efficiency of the fluorescence, i.e. emitted photons/absorbed photons
         *  - \c p_x, \c p_y, \c p_z gives the direction vector of the dipole moment of the walker, note that
         *                           the program excepcts \f$ \|\vec{p}\|=\sqrt{p_x^2+p_y^2+p_z^2}=1 \f$
         */
        struct walkerState {
            /** \brief indicates whether this fluorophor exists (should be accounted for in simulation/detection) */
            bool exists;
            /** \brief current timestep */
            uint64_t time;
            /** \brief x-position of the walker [micrometers] */
            double x;
            /** \brief y-position of the walker [micrometers] */
            double y;
            /** \brief z-position of the walker [micrometers] */
            double z;
            /** \brief starting x-position of the walker [micrometers] */
            double x0;
            /** \brief starting y-position of the walker [micrometers] */
            double y0;
            /** \brief starting z-position of the walker [micrometers] */
            double z0;

            /** \brief x-position of the walker [gridpositions] */
            int32_t ix;
            /** \brief y-position of the walker [gridpositions] */
            int32_t iy;
            /** \brief z-position of the walker [gridpositions] */
            int32_t iz;
            /** \brief starting x-position of the walker [gridpositions] */
            int32_t ix0;
            /** \brief starting y-position of the walker [gridpositions] */
            int32_t iy0;
            /** \brief starting z-position of the walker [gridpositions] */
            int32_t iz0;

            /** \brief internal (quantum-mechanical) state of the walker: 0 (ground state), 1 (excited state) */
            int qm_state;
            /** \brief absorbption cross section [m^2] */
            double sigma_abs[N_FLUORESCENT_STATES];
            /** \brief quantum efficiency of fluorescence [0..1]*/
            double q_fluor[N_FLUORESCENT_STATES];
            /** \brief photophysics transition rates in [1/second] The transition probabilities will be calculated as
             *         rate*sim_timestep. So keep in mind, that all of these propabilitites have to sum up to 1 for one state */
            double photophysics_transition[N_FLUORESCENT_STATES*N_FLUORESCENT_STATES];
            /** \brief number of used QM states (has to be smaller than  N_FLUORESCENT_STATES */
            int used_qm_states;
            /** \brief x-component of dipole moment vector */
            double p_x;
            /** \brief y-component of dipole moment vector */
            double p_y;
            /** \brief z-component of dipole moment vector */
            double p_z;
            /** \brief if you have multiple types of fluorophores, you can specify the type here (not interpreted internally, for user only!!!) */
            int type;
            /** \brief this specifies the absorption spectrum to use for the fluorophor see get_spectral_efficiency() for details */
            int spectrum;

            /** \brief indicates, that this object has JUST been reset */
            bool was_just_reset;

            walkerState& operator=(walkerState& other);
            walkerState(walkerState& other);
        };

        /** \brief the possible shapes of the simulational volume */
        enum VolumeShape {
            Box,
            Ball
        };

        enum AdditionalWalkerPositions {
            SamePosition,
            InSphere
        };

    protected:
        /** \brief pointer to the fluorophor database */
        FluorophorManager* fluorophors;


        /** \brief an array which holds the states of all walkers for potentially multiple timesteps*/
        walkerState* walker_state;

        /** \brief an array which holds the states of all walkers for potentially multiple timesteps*/
        walkerState* walker_state_other;

        /** \brief x-displacement vector for additional walksers */
        double* walker_dx;
        /** \brief y-displacement vector for additional walksers */
        double* walker_dy;
        /** \brief z-displacement vector for additional walksers */
        double* walker_dz;

        /** \brief if set \c true this class supports two vector of walker_states that may be interchanged! otherwise they
         *         point to the same data */
        bool use_two_walkerstates;

        /** \brief number of simulated walkers */
        unsigned long walker_count;


        /** \brief fluorophor concentration in [nanomol/l] */
        double c_fluor;

        /** \brief x-width of simulation volume in [micrometers]
         *
         * When \c volume_shape==0 the simulational box is a rectangular box spanning the space from
         * \c [0..sim_x]*[0..sim_y]*[0..sim_z].
         */
        double sim_x;
        /** \brief y-width of simulation volume in [micrometers]
         *
         * When \c volume_shape==0 the simulational box is a rectangular box spanning the space from
         * \c [0..sim_x]*[0..sim_y]*[0..sim_z].
         */
        double sim_y;
        /** \brief z-width of simulation volume in [micrometers]
         *
         * When \c volume_shape==0 the simulational box is a rectangular box spanning the space from
         * \c [0..sim_x]*[0..sim_y]*[0..sim_z].
         */
        double sim_z;
        /** \brief radius of simulation volume in [micrometers]
         *
         * When \c volume_shape==1 the simulational box is a sphere around 0 with radius sim_radius
         */
        double sim_radius;
        /** \brief simulation timestep in [seconds] */
        double sim_timestep;
        /** \brief the overall absolute simulation time [seconds] */
        double sim_time;

        /** \brief gives the duration of the measurement */
        double duration;

        /** \brief if this is \c true the photophysics (bleaching, triplet, fluorescence lifetime) is also simulated */
        bool use_photophysics;

        /** \brief number of fluorophores attatched to each particle (each fluorophore may have it's own photophysics and own position (fixed, relative to position of central walker) */
        int n_fluorophores;


        /** \brief initial internal (quantum-mechanical) state of the walker: 0 (ground state), 1 (excited state), 2 (transition from 1 to 0), -1 (triplet state), -2 (bleached) */
        int init_qm_state;
        /** \brief initial value for number of used QM states */
        int init_used_qm_states;
        /** \brief initial type of the walker */
        int init_type;
        /** \brief initial absorption spectrum ID of the walker */
        int init_spectrum;
        /** \brief initial absorbption cross section [m^2] */
        double init_sigma_abs[N_FLUORESCENT_STATES];
        /** \brief initial quantum efficiency of fluorescence [0..1]*/
        double init_q_fluor[N_FLUORESCENT_STATES];
        /** \brief initial photophysics transition rates */
        double init_photophysics_transition[N_FLUORESCENT_STATES*N_FLUORESCENT_STATES];
        /** \brief initial x-component of dipole moment vector */
        double init_p_x;
        /** \brief initial y-component of dipole moment vector */
        double init_p_y;
        /** \brief initial z-component of dipole moment vector */
        double init_p_z;

        /** \brief initial fluorophor name */
        std::string init_fluorophor;

        /** \brief set \c true, while simulation is heating up */
        bool heating_up;


        /** \brief this string may be used to identify the object */
        std::string object_name;


        /** \brief shape of the simulational volume
         *
         * - \c 0 = box  \c [0..sim_x]*[0..sim_y]*[0..sim_z]
         * - \c 1 = sphere around 0 with radius sim_radius
         */
        VolumeShape volume_shape;

        /** \brief GSL helper object: random number generator type */
        const gsl_rng_type * rng_type;

        /** \brief GSL helper object: random number generator */
        gsl_rng * rng;

        virtual void handle_parent_walker_count_changed(unsigned long N_walker, unsigned long N_fluorophores);

        /** \brief all classes in this list are notified (by calling their handle_parent_walker_count_changed() method) if the walker count in this class changed. \note The hooks are called AFTER iit_walkers() !!! */
        std::vector<FluorophorDynamics*> notify_when_walkercount_changes;

        /** \brief if the photo-physics depends on the illumination intensity/relative absorbance, the readers in this list are used to calculate the relative absorbance of a fluorophore at a given position */
        std::string absorbance_reader;

        RelativeAbsorbanceReader* get_abs_reader() const;

        /** \brief indicates, whether the photo physics depends on the relative absorbance */
        bool photophysics_absorbance_dependent;
        /** \brief the relative absorbance, multiplied by this factor is the factor for the jump-probabilities  */
        double photophysics_absorbance_factor;


        /** \brief set the number of walkers and allocate the according amount of memory for walker_state */
        void change_walker_count(unsigned long N_walker, unsigned long N_fluorophores);
        /** \brief get the number of walkers for the current simulation settings */
        virtual unsigned long calc_walker_count();


        /** \brief if this is >0 the object will protocoll as many of the particle trajectories as given in protocol_trajectories.
         *
         * The propagate method will store the first protocol_trajectories trajectories in trajectories. The method save_trajectories()
         * may be used to write the stored trajectories into files.
         */
        unsigned int protocol_trajectories;

        /** \brief gives the maximum number of timesteps to store in a trajectory. if this is -1 all timesteps will be stored (default) */
        long long protocol_timestep_count;

        /** \brief file handles when saving trajectories */
        FILE** trajectoryFile;

        /** \brief the basename of the current simulation */
        std::string basename;

        /** \brief group name of this object */
        std::string group;
        /** \brief supergroup name of this object */
        std::string supergroup;
        /** \brief number of the object */
        int object_number;
        /** \brief trajectory store */
        std::vector<std::vector<walkerState> > trajectories;

        /** \brief stores the current walker states in the protocol for the first protocol_trajectories walkers. */
        void store_step_protocol();

        /** \brief stores the current walker states in the protocol for the first protocol_trajectories walkers that lie between \a w_start and \a w_end. */
        //void store_step_protocol(int w_start, int w_end);

        /** \brief used to read configuration data from ini file ... overwrite this to read data in derived classes.
         *
         * This base class will care for entering the right group to read from! But keep in mind that you will HAVE to call the
         * corresponding function of the parent class in your implementation!
         */
        virtual void read_config_internal(jkINIParser2& parser);

        /** \brief this is true as long, as this object may produce trajectories, set this to \c false if you reached the end of let's say
         *         a trajectory input file. */
        bool endoftrajectory;
        /** \brief depletion propability with this propability a particle does not re-enter the simulation box */
        double depletion_propability;
        /** \brief reset quantum state if particle leaves simulation box, this allows to perform bleaching with a reservoir of functional fluorophores or depletion
         *
         *  Usually bleaching can be implemented by a photophysics rate into a dark state that can n ot be left anymore. If this is set \c true, the
         *  quantum state (i.e. photophysics) of a particle is reset, when it leaves the simulation box and is therefore reinitialized at a new starting position.
         *  If this is set \c false (default), the fluorophores stays switched off if it leaves the volume in the off state and is never switched on again (unless
         *  the derived dynamics class does it explicitly, of course). This last case resembles the depletion of a reservoir of fluorophores.
         */
        bool reset_qmstate_at_simboxborder;

        /** \brief is this and use_photophysics are BOTH \c true, each additionalö walker will have it's own blinking dynamics */
        bool additional_walker_photophysics;
        /** \brief is this is \c true, the additional walkers are set non-existent if main-walker is non-existent */
        bool additional_walker_off_if_main_off;
        /** \brief how should the additional walkers be positioned */
        AdditionalWalkerPositions additional_walker_position_mode;
        /** \brief radius of sphere if additional_walker_position_mode==InSphere */
        double additional_walker_sphere_radius;

        /** \brief indicates whether walker statistics is stored */
        bool store_walker_statistics;

        /** \brief average statistics of this duration [s] for every entry */
        double walker_statistics_averageduration;
        double walker_statistics_nextsavetime;

        /** \brief number if simulation steps to heat up the simulation ... */
        int64_t heatup_steps;

        /** \brief struct for the walker statistics */
        struct walker_statistics_entry {
            walker_statistics_entry(double t=0.0) {
                clear(t);
            };
            inline void clear(double t=0.0) {
                time=t;
                count_all=0;
                count_existing=0;
                average_brightness=0;
                average_steps=0;
                posx=posy=posz=0;
                count_existing_reallyinside=0;
                for (int i=0; i<N_FLUORESCENT_STATES; i++) state_distribution[i]=0;
            }
            double time;
            uint64_t count_all;
            uint64_t count_existing;
            uint64_t count_existing_reallyinside;
            double average_brightness;
            double state_distribution[N_FLUORESCENT_STATES];
            uint64_t average_steps;
            double posx, posy, posz;
        };

        std::vector<walker_statistics_entry> walker_statistics;

        /** \brief save the results of the measurement
         */
        virtual void save();
    public:
        /** \brief class constructor with standard volume 30*30*30micron^3 and a concentration of 1nM */
        FluorophorDynamics(FluorophorManager* fluorophors, std::string object_name=std::string(""));


        /** \brief class constructor, initialises with a box of the given dimensions */
        FluorophorDynamics(FluorophorManager* fluorophors, double sim_x, double sim_y, double sim_z, double c_fluor, std::string object_name=std::string(""));

        /** \brief class constructor, initialises with a sphere of the given dimensions */
        FluorophorDynamics(FluorophorManager* fluorophors, double sim_radius, double c_fluor, std::string object_name=std::string(""));

        /** \brief class destructor */
        virtual ~FluorophorDynamics();

        /** \brief read configuration from INI file
         *
         * \param parser the parser object to read the data from
         * \param group The configuration data is read from the given group in the  parser
         * \param supergroup If a supergroup is provided this method first readsthe configuration from the super group and afterwards from the group.
         *
         * The \a supergroup parameter allows to use a stacked configration scheme: A super group contains parameters that are used for all simulation
         * objects of the same type while the \a group may be used to further specialize these parameters for every single object.
         */
        virtual void read_config(jkINIParser2& parser, std::string group=std::string("dynamics"), std::string supergroup=std::string(""));

        /** \brief set the simulation timestep */
        virtual void set_sim_timestep(double value);

        /** \brief set the fluorophor concentration */
        virtual void set_c_fluor(double value);

        /** \brief set the dimensions of the simulational box */
        virtual void set_sim_box(double vx, double vy, double vz);

        /** \brief set the dimensions of the simulational sphere */
        virtual void set_sim_sphere(double rad);

        /** \brief initialize the state of the i-th walker and put it to the given position. The walker step counter is reset to 0 */
        virtual void init_walker(unsigned long i, double x=0, double y=0, double z=0);

        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();
        /** \brief heats up the simulation */
        void run_heatup();

        /** \brief propagate all walkers in the simulation one timestep further
         *
         * \param boundary_check when \c boundary_check==false the walkers won't be reset when they reach
         *                       the border of the simulational box
         *
         */
        virtual void propagate(bool boundary_check=true);

        /** \brief implements the photophysics. Overwrite this if you want higher oder photophysics, other than fluorescence lifetime, triplet dynamics and photobleaching
         *
         * \param walker the walker to do the simulation for
         *
         * \note this function internally calls propagate_photophysics_scaled(), in the basic implementation with scale=1.0!
         *
         */
        virtual void propagate_photophysics(int walker);

        /** \brief implements the photophysics. Overwrite this if you want higher oder photophysics, other than fluorescence lifetime, triplet dynamics and photobleaching
         *
         * \param walker the walker to do the simulation for
         * \param scale a scaling factor for the jump probabilities, i.e. if scale<1, the jump probabilities will be smaller than initially and in the opposite case larger.
         *              make sure, that if scale>1, the summed probability does not get >1!!!
         *
         */
        virtual void propagate_photophysics_scaled(int walker, double scale=1.0);

        /** \brief propagates the additional walkers
         *
         */
        virtual void propagate_additional_walkers();


        /** \brief estimates teh runtime of the simulation, or return 0 if the simulation requires a given runtime */
        virtual double estimate_runtime();

        double get_walker_sigma_times_qfl(unsigned long i);


        /** \brief return the number of walkers in the simulational box */
        unsigned long get_walker_count() ;

        /** \brief get pointer to array with all walker states */
        inline walkerState* get_walker_state() { return walker_state_other; }


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


        /** \brief perform a boundary check for the i-th walker and reset it to a random border position, if it left the sim box */
        virtual void perform_boundary_check(unsigned long i);

        /** \brief returns true, if the given particle is inside the simulation box and not in its outer laiwer with width (rel_margin * diameter in each direction) */
        bool reallyInsideVolume(double x, double y, double z, double rel_margin=0.01) const;



        /** \brief returns \c true if the end of the possible trajectories is reached, i.e. as long as this object may
         *         suply data this is \c false */
        inline bool end_of_trajectory_reached() {
            return endoftrajectory;
        }

        /** \brief this function is called once after the simulation has finished.
         */
        virtual void finalize_sim();


        /** \brief report the state of the object, as human-readable text */
        virtual std::string report();

        /** \brief returns the simulational timestep */
        inline double get_sim_timestep() { return sim_timestep; };

        GET_MACRO(VolumeShape, volume_shape);
        GET_MACRO(double, sim_radius);
        GET_MACRO(double, c_fluor);
        GET_SET_MACRO(int, init_qm_state);
        inline double get_init_sigma_abs(int i) {
            if ((i>=0)&&(i<N_FLUORESCENT_STATES)) return init_sigma_abs[i];
            return 0;
        };
        inline double get_init_q_fluor(int i) {
            if ((i>=0)&&(i<N_FLUORESCENT_STATES)) return init_q_fluor[i];
            return 0;
        };
        GET_SET_MACRO(double, init_p_x);
        GET_SET_MACRO(double, init_p_y);
        GET_SET_MACRO(double, init_p_z);
        GET_SET_MACRO(double, depletion_propability);
        GET_SET_MACRO(int, init_type);
        GET_SET_MACRO(int, init_spectrum);
        //GET_SET_MACRO(bool, test_spectra);
        GET_SET_MACRO(bool, heating_up);
        GET_SET_MACRO(bool, use_photophysics);
        GET_SET_MACRO(std::string, basename);
        GET_SET_MACRO(std::string, object_name);
        GET_SET_MACRO(int, object_number);

        GET_MACRO(std::string, group);
        GET_MACRO(std::string, supergroup);
        GET_MACRO(bool, use_two_walkerstates);
        GET_MACRO(double, sim_time);
        GET_MACRO(std::string, absorbance_reader);

        GET_SET_MACRO(bool, additional_walker_photophysics)
        GET_SET_MACRO(bool, additional_walker_off_if_main_off)
        GET_SET_MACRO_I(AdditionalWalkerPositions, additional_walker_position_mode, init_walkers())
        GET_SET_MACRO_I(double, additional_walker_sphere_radius, init_walkers())


        inline void set_use_two_walkerstates(bool v) {
            use_two_walkerstates=v;
            change_walker_count(calc_walker_count(), n_fluorophores);
        };



        /** \brief this method saves the stored trajectories (see save_trajectories) in files \c traj001.dat ... */
        void save_trajectories();

        GET_SET_MACRO(unsigned int, protocol_trajectories);


//        /** \brief get the absorption efficiency [0..1] from the given walker at the given wavelength [nm] at timestep 0
//         *
//         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
//         * from the given spectrum.
//         */
//        inline double get_walker_spectral_absorbance(unsigned long walker, double wavelength) {
//            int spectrum=walker_state_other[walker].spectrum;
//            return fluorophors->get_spectral_absorbance(spectrum, wavelength);
//        }
//
//        /** \brief get the fluorescence value [0..1] from the given walker at the given wavelength [nm] at timestep 0
//         *
//         * If \a spectrum equals -1 this function simply returns 1, otherwise an interpolated value
//         * from the given spectrum.
//         */
//        inline double get_walker_spectral_fluorescence(unsigned long walker, double wavelength) {
//            int spectrum=walker_state_other[walker].spectrum;
//            return fluorophors->get_spectral_fluorescence(spectrum, wavelength);
//        }

        /** \brief swap the contents of the array that receives the new simulation data and the array used to read the fluorophor states.
         *         after a call the write array will contain the same data as the read array */
        inline void swap_walker_states() {
            walkerState* temp=walker_state_other;
            walker_state_other=walker_state;
            walker_state=temp;
            memcpy(walker_state, walker_state_other, walker_count*sizeof(walkerState));
        }

        /** \brief returns the estimated memory consumption in bytes */
        virtual inline double calc_mem_consumption() { return walker_count*sizeof(walker_state); };


        /** \brief tests the dynamics simulation by computing some walker steps and theroy results
         */
        virtual void test(unsigned int steps=1000, unsigned int walkers=100);

        /** \brief tests the photophysics simulation
         *
         *  This does NOT use the scaling by intensity feature!!!
         */
        virtual void test_photophysics(unsigned int steps=1000, unsigned int walkers=100);

        /** \brief goes through all walkers and loads the spectra used by these walkers. */
        virtual void load_all_used_spectra();

        /** \brief initialize all walkers to the default settings */
        void init_walkers();
        /** \brief initialize all additional walkers to the default settings */
        void init_additional_walkers();

        /** \brief save the results of the measurement
         */
        virtual void save_results();

        /*! \brief returns \c true , if this object depends on the given \ other object.

            This function is used to determine the order in which the propagate() method of all dynamics objects in this
            simulation is called.
         */
        virtual bool depends_on(const FluorophorDynamics* other) const;

        void ensure_dynamics_is_hooked(FluorophorDynamics* other);

        /*! \brief return DOT-code (for GraphViz) that represents the node without links (these are returned by dot_get_links() ).
        */
        virtual std::string dot_get_node(bool simple=false) ;
        /*! \brief return DOT-code (for GraphViz) that represents the links of the nodes returned by dot_get_node()
        */
        virtual std::string dot_get_links() ;
                /*! \brief return DOT label-code (for GraphViz) that represents the nodes properties
        */
        virtual std::string dot_get_properties() ;
};

#endif // FLUOROPHORDYNAMICS_H
