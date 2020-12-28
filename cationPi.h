//
// Created by Yiming Tang on 2020/11/19.
//

#ifndef GMX_CATIONPI_CG_CATIONPI_H
#define GMX_CATIONPI_CG_CATIONPI_H

#include "gromacs/trajectoryanalysis.h"
#include <string>
#include <vector>


using namespace std;
using namespace gmx;

struct coordinate
{
    float x;
    float y;
    float z;
};

class cationPi: public TrajectoryAnalysisModule
{
public:
    cationPi();

    virtual void initOptions(IOptionsContainer          *options,
                             TrajectoryAnalysisSettings *settings);

    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation        &top);

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);

    virtual void finishAnalysis(int nframes);

    virtual void writeOutput();

private:

    static float calDistance(coordinate atom11, coordinate atom12, coordinate atom13,
                             coordinate atom2);


    gmx::AnalysisData             data_probability_;
    gmx::Selection                sel_ions_;
    gmx::Selection                sel_rings_;

    gmx::AnalysisNeighborhood     nb_;

    double cutoff_;

    unsigned long *probability_;

    std::string     fnPropability_;

    double          maxDistance_;
    double          stepDistance_;

    // The following parameters are only valid if the sel_ and ref_ selection groups are identical.
    // And, if the inter_molecule_ is true, the ring_number_in_molecule must be specified.
    // Vice versa, if the ring_number_in_molecule_ is specified, the inter_molecule_ must be true.

    bool            inter_molecule_;
    bool            intra_molecule_;
    int             ring_number_in_molecule_;
    int             ion_number_in_molecule_;

    t_topology     *top_;
    t_atoms         atoms_;

    bool **ring_has_contact_;

};


#endif //GMX_CATIONPI_CG_CATIONPI_H
