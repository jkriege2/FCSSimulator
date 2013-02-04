#ifndef MSDMEASUREMENT_H
#define MSDMEASUREMENT_H

#include "fluorescencemeasurement.h"
#include "diffusiontools.h"
#include "../../../LIB/trunk/jkiniparser2.h"
#include <string>
#include <vector>
#include "../../../LIB/trunk/tools.h"
#include "fluorophordynamics.h"
#include "../../../LIB/trunk/multitau-msd.h"

/*! \brief MSD measurement class
    \ingroup diff4_measurement

   This class implements an evaluation that measures the MSD of the trajectories supplied to it.

 */
class MSDMeasurement : public FluorescenceMeasurement
{
    public:
        MSDMeasurement(FluorophorManager* fluorophors, std::string objectname=std::string(""));
        virtual ~MSDMeasurement();



        /** \brief initialize the simulation */
        virtual void init();

        /** \brief propagate the detection one step further */
        virtual void propagate();

        /** \brief report the object state */
        virtual std::string report();

    protected:

        virtual void save();
        /** \brief read configuration from INI file */
        virtual void read_config_internal(jkINIParser2& parser);
        /** \brief clear all internal data structures */
        void clear();


        /** \brief number of trajectories for which to calculate the MSD */
        int msd_for_trajectories;

        /** \brief decades in MSD */
        int msd_s;
        /** \brief time lags in one decade */
        int msd_p;
        /** \brief time factor between MSD decades */
        int msd_m;

        std::vector<MultiTauMSD<double>* > msds;
        struct trajectory_info {
            double sum_x;
            double sum2_x;
            double sum_y;
            double sum2_y;
            double sum_z;
            double sum2_z;
            uint64_t cnt;
            double xmin;
            double xmax;
            double ymin;
            double ymax;
            double zmin;
            double zmax;
        };
        std::vector<trajectory_info> trajectoryinfo;
    private:
};

#endif // MSDMEASUREMENT_H
