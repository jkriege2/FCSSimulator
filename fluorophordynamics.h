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

#include "../../../LIB/trunk/tools.h"
#include "../../../LIB/trunk/jkINIParser2.h"
#include "diffusiontools.h"


#ifndef FLUOROPHORDYNAMICS_H
#define FLUOROPHORDYNAMICS_H






/** \brief this is the virtual base class for any class describing fluorophor dynamics
 *  \ingroup diff4_dynamic
 *
 * Basically this class provides an interface that allows to access a set of fluorescent particles which
 * are described (see FluorophorDynamics::walkerState) by a position within the simulational volume and
 * their internal (quantum-mechanical and photo-physical) properties (absorption cross section, fluorescence
 * efficiency, QM state, ...).
 *
 * If you want to implement a dynamic you will have to write a class that publicly inherits FluorophorDynamics
 * and mainly overwrites the propagate() and init() methods. If you don't wan to implement special features that
 * go beyond the things implemented here, you won't have to overwrite other methods. You can signal other classes
 * that your dnymaics implementation may be used in a multi-threaded way, by setting the protected data memeber
 * multi_threadable to \c true in your constructur. In this case you will have to provide both varieties of the
 * propagate() method, so another class may call \c propagate(w_start, w_end) with different parameters in parallel!
 *
 * If you need additional parameters from an ini file you will also have to reimplement read_config().
 *
 * \section FD_rng Random Number Generators
 * This base class already instaciates a GSL random number generator (taus2 as default) in the data member rng
 * which you may use for your implementation. By setting the rng property in the ini file you can select one of
 * these generators (for details, see the GSL documentation):
 * <center><tt>mt19937, ranlxs0, ranlxs1, ranlxs2, ranlxd1, ranlxd2, ranlux, ranlux389, cmrg, mrg, taus, taus2, gfsr4, rand, bsd, libc5, glibc2, rand48, ranf, ranmar, r250, tt800, minstd, knuthran2, knuthran2002, knuthran</tt></center>
 *
 * \section FD_init Initialisation of Simulation Environment
 * By giving the volume \f$ V \f$ of the simulational box and the fluorophor concentration \f$ c_f \f$ the method
 * change_walker_count() will initialize enough memory for \f$ \mbox{ceil}(v\cdot c_f) \f$ particles, set
 * then at random positions in the volume and init their state with the given \c init_... values.
 *
 *
 * \section FD_photophysics Photophysics Simulation
 * This class implements basic photophysics, if activated with use_photophysics.
 * The photophysics is implemented in propagate_photophysics() which you will have to call for every walker when implementing
 * the dynamics.
 *
 * The photophysics simulation implements triplet dynamics and photo bleaching by means of the
 * qm_state property of a walker. These states are used:
 *   - \b 0: ground state (excitable to excited state and triplet state; bleachable
 *   - \b -1: in triplet state: may return to 0 after a triplet lifetime (exponential decay)
 *   - \b -2: bleached: may not return to any state
 * .
 * The impact of these states on the fluorescence is modeled by the function get_walker_state_factor() which returns a
 * number between 0 and 1. In the implementation herein it simply returns 1 if \c qm_state==0 and 0 otherwise. This
 * factor may be multiplied with the absorbtion propability. If it is 0 the current walker may NOT contribute to the
 * fluorescence signal.
 *
 * So basically this implementation assumes that the fluorescence dynamics, i.e. the transition from ground to excited
 * state and back is much faster than the simulation timestep. I.e. you won't be able to see the antibunching part of
 * the photon correlation.
 *
 *
 * \section FD_threads threading support
 *
 * By setting \a use_two_walkerstates to \c true the class uses two arrays of walker_state's which may be interchanged by calling
 * swap_walker_states(). In this case the method propagate() writes to walker_state and all read functions use walker_stat_other
 * so it is possible to fill walker_state while the old walker_state_other is beeing read. If \a use_two_walkerstates is set to \c false
 * the two variables simply point to the same memory location.
 *
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
            /** \brief internal (quantum-mechanical) state of the walker: 0 (ground state), 1 (excited state), 2 (transition from 1 to 0), -1 (triplet state), -2 (bleached) */
            int qm_state;
            /** \brief absorbption cross section [m^2] */
            double sigma_abs;
            /** \brief quantum efficiency of fluorescence [0..1]*/
            double q_fluor;
            /** \brief fluorescence lifetime [seconds] */
            double tau_fl;
            /** \brief bleaching propability [0..1] */
            double bleaching_propability;
            /** \brief triplet lifetime [sec] */
            double triplet_lifetime;
            /** \brief triplet transition propability [0..1] */
            double triplet_propability;
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
        };

        /** \brief the possible shapes of the simulational volume */
        enum VolumeShape {
            Box,
            Ball
        };

    protected:
        /** \brief pointer to the fluorophor database */
        FluorophorManager* fluorophors;


        /** \brief an array which holds the states of all walkers for potentially multiple timesteps*/
        walkerState* walker_state;

        /** \brief an array which holds the states of all walkers for potentially multiple timesteps*/
        walkerState* walker_state_other;

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


        /** \brief initial internal (quantum-mechanical) state of the walker: 0 (ground state), 1 (excited state), 2 (transition from 1 to 0), -1 (triplet state), -2 (bleached) */
        int init_qm_state;
        /** \brief initial type of the walker */
        int init_type;
        /** \brief initial absorption spectrum ID of the walker */
        int init_spectrum;
        /** \brief initial absorbption cross section [m^2] */
        double init_sigma_abs;
        /** \brief initial quantum efficiency of fluorescence [0..1]*/
        double init_q_fluor;
        /** \brief initial fluorescence lifetime [seconds]*/
        double init_tau_fl;
        /** \brief initial bleaching propability [0..1] */
        double init_bleaching_propability;
        /** \brief initial triplet lifetime [sec] */
        double init_triplet_lifetime;
        /** \brief initial triplet transition propability [0..1] */
        double init_triplet_propability;
        /** \brief initial x-component of dipole moment vector */
        double init_p_x;
        /** \brief initial y-component of dipole moment vector */
        double init_p_y;
        /** \brief initial z-component of dipole moment vector */
        double init_p_z;


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

        /** \brief set the number of walkrs and allocate the according amount of memory for walker_state */
        virtual void change_walker_count(unsigned long N_walker);

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

    public:
        /** \brief class constructor with standard volume 30*30*30µm^3 and a concentration of 1nM */
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
        inline void set_sim_timestep(double value) {
            sim_timestep=value;
        };

        /** \brief set the fluorophor concentration */
        inline void set_c_fluor(double value) {
            c_fluor=value;
            if (volume_shape==Box) {
                change_walker_count((unsigned long)round(c_fluor*1e-9*6.022e23*sim_x*1e-5*sim_y*1e-5*sim_z*1e-5));
            } else if (volume_shape==Ball) {
                change_walker_count((unsigned long)round(c_fluor*1e-9*6.022e23*4.0*M_PI/3.0*gsl_pow_3(sim_radius*1e-5)));
            }
        };

        /** \brief set the dimensions of the simulational box */
        inline void set_sim_box(double vx, double vy, double vz) {
            volume_shape=Box;
            sim_x=vx;
            sim_y=vy;
            sim_z=vz;
            change_walker_count((unsigned long)round(c_fluor*1e-9*6.022e23*sim_x*1e-5*sim_y*1e-5*sim_z*1e-5));
        };

        /** \brief set the dimensions of the simulational sphere */
        inline void set_sim_sphere(double rad) {
            volume_shape=Ball;
            sim_radius=rad;
            change_walker_count((unsigned long)round(c_fluor*1e-9*6.022e23*4.0*M_PI/3.0*gsl_pow_3(sim_radius*1e-5)));
        };

        /** \brief return the number of walkers in the simulational box */
        inline unsigned long get_walker_count() { return walker_count; };

        /** \brief initialize the simulation environment (random walker positions ... */
        virtual void init();

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
         */
        virtual void propagate_photophysics(int walker);

        /** \brief get state of i-th walker */
        inline walkerState get_walker_state(unsigned long i) { return walker_state_other[i]; };

        /** \brief get pointer to array with all walker states */
        inline walkerState* get_walker_state() { return walker_state_other; };

        /** \brief copy the array walker_state to the given location (this assumes that memory has already been reserved!)
         *
         * \return a pointer to the first memory location \b behind the copy
         */
        virtual walkerState* copy_walker_state(walkerState* start);

//        /** \brief set state of i-th walker*/
//        inline virtual void set_walker_state(unsigned long i, walkerState state) { walker_state_other[i]=state; };
//
//
//        /** \brief get position of i-th walker */
//        inline virtual void get_walker_position(unsigned long i, double* x, double* y, double* z) { *x=walker_state_other[i].x; *y=walker_state_other[i].y; *z=walker_state_other[i].z; };
//
//        /** \brief get dipole moment vector of i-th walker */
//        inline virtual void get_walker_dipole(unsigned long i, double* p_x, double* p_y, double* p_z) { *p_x=walker_state_other[i].p_x; *p_y=walker_state_other[i].p_y; *p_z=walker_state_other[i].p_z; };
//
//        /** \brief get the quantum mechanical state of the i-th walker at a give timestep */
//        inline virtual int get_walker_qm_state(unsigned long i) { return  walker_state_other[i].qm_state; };
//
//        /** \brief get the type of the i-th walker */
//        inline virtual int get_walker_type(unsigned long i) { return  walker_state_other[i].type; };
//
//        /** \brief get the spectrum of the i-th walker */
//        inline virtual int get_walker_spectrum(unsigned long i) { return  walker_state_other[i].spectrum; };
//
//        /** \brief get the absorption cross section of the i-th walker */
//        inline virtual double get_walker_sigma_abs(unsigned long i) { return  walker_state_other[i].sigma_abs; };
//
//        /** \brief get the quantum efficiency of the fluorescence of the i-th walker */
//        inline virtual double get_walker_q_fluor(unsigned long i) { return  walker_state_other[i].q_fluor; };
//
//        /** \brief get the fluorescence lifetime of the i-th walker */
//        inline virtual double get_walker_tau_fl(unsigned long i) { return  walker_state_other[i].tau_fl; };
//
//        /** \brief get the bleaching propability of the i-th walker [0..1] */
//        inline virtual double get_walker_bleaching_propability(unsigned long i) { return  walker_state_other[i].bleaching_propability; };
//
//        /** \brief get the triplet lifetim of the i-th walkere [sec] */
//        inline virtual double get_walker_triplet_lifetime(unsigned long i) { return  walker_state_other[i].triplet_lifetime; };
//
//        /** \brief get the triplet transition propability of the i-th walker [0..1] */
//        inline virtual double get_walker_triplet_propability(unsigned long i) { return  walker_state_other[i].triplet_propability; };
//
//        /** \brief returns whether the given walker really exists */
//        inline virtual double get_walker_exists(unsigned long i) { return  walker_state_other[i].exists; };
//
//        /** \brief get a factor [0..1] for the i-th walker which describes the dependence on the inner state of the walker */
//        inline virtual double get_walker_state_factor(unsigned long i) {
//            if (!use_photophysics) return 1.0;
//            return (walker_state_other[i].qm_state==0)?1.0:0.0;
//        };

        /** \brief returns \c true if the end of the possible trajectories is reached, i.e. as long as this object may
         *         suply data this is \c false */
        inline bool end_of_trajectory_reached() {
            return endoftrajectory;
        }


        /** \brief report the state of the object, as human-readable text */
        virtual std::string report();

        /** \brief returns the simulational timestep */
        inline double get_sim_timestep() { return sim_timestep; };

        GetMacro(VolumeShape, volume_shape);
        GetMacro(double, sim_radius);
        GetMacro(double, c_fluor);
        GetSetMacro(int, init_qm_state);
        GetSetMacro(double, init_sigma_abs);
        GetSetMacro(double, init_q_fluor);
        GetSetMacro(double, init_p_x);
        GetSetMacro(double, init_p_y);
        GetSetMacro(double, init_p_z);
        GetSetMacro(int, init_type);
        GetSetMacro(int, init_spectrum);
        //GetSetMacro(bool, test_spectra);
        GetSetMacro(bool, use_photophysics);
        GetSetMacro(std::string, basename);
        GetSetMacro(std::string, object_name);
        GetMacro(bool, use_two_walkerstates);
        GetMacro(double, sim_time);

        void set_use_two_walkerstates(bool v) {
            use_two_walkerstates=v;
            change_walker_count(walker_count);
        };



        /** \brief this method saves the stored trajectories (see save_trajectories) in files \c traj001.dat ... */
        void save_trajectories();

        GetSetMacro(unsigned int, protocol_trajectories);


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
        virtual void test(unsigned int steps=1000, unsigned int walkers=100) {};

        /** \brief goes through all walkers and loads the spectra used by these walkers. */
        virtual void load_all_used_spectra();
};

#endif // FLUOROPHORDYNAMICS_H
