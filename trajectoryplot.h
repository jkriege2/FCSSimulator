#ifndef TRAJECTORYPLOT_H
#define TRAJECTORYPLOT_H

#include "fluorescencemeasurement.h"
#include "diffusiontools.h"
#include "../../../LIB/trunk/jkiniparser2.h"
#include <string>
#include <vector>
#include <map>
#include "../../../LIB/trunk/tools.h"
#include "fluorophordynamics.h"
#include "../../../LIB/trunk/multitau-msd.h"

/*! \brief measurement class that plots trajectories
    \ingroup diff4_measurement


 */
class TrajectoryPlotMeasurement : public FluorescenceMeasurement {
    public:
        TrajectoryPlotMeasurement(FluorophorManager* fluorophors, std::string objectname=std::string(""));
        virtual ~TrajectoryPlotMeasurement();



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


        /** \brief number of trajectories to plot */
        int trajectories_to_plot;

        /** \brief average over this many steps for each datapoint */
        int avg_steps;
        
        /** \brief read at most this number of steps from a walker */
        int max_steps;
        
        struct trajectory_info {
            float t;
            float x;
            float y;
            float z;
        };
        std::vector<std::vector<trajectory_info> > trajectories;
        std::map<int, trajectory_info> currentT;
        int avgCount;
    private:
};

#endif // TRAJECTORYPLOT_H
