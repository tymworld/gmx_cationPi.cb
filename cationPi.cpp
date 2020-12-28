//
// Created by Yiming Tang on 2020/11/19.
//

#include "cationPi.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <gromacs/fileio/matio.h>
#include <gromacs/fileio/gmxfio.h>

using namespace std;
using namespace gmx;


float cationPi::calDistance(const coordinate atom11, const coordinate atom12, const coordinate atom13,
                            const coordinate atom2)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x) / 3;
    y_1 = (atom11.y + atom12.y + atom13.y) / 3;
    z_1 = (atom11.z + atom12.z + atom13.z) / 3;

    x_2 = atom2.x;
    y_2 = atom2.y;
    z_2 = atom2.z;

    return (float)sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}


cationPi::cationPi(): cutoff_(1.5), maxDistance_(1.2), stepDistance_(0.01)
{
    registerAnalysisDataset(&data_probability_, "probability");
}

void cationPi::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
    {
            "Distance probability (PDF) between ions and center of centroids of rings."
    };

    settings->setHelpText(desc);


    options->addOption(FileNameOption("raw-probability").filetype(eftGenericData).outputFile()
                               .store(&fnPropability_)
                               .description("Raw Probability Data for further analysis"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                               .description("Ring and ions whose maximum atom-wise distance are beyond this cutoff is not considered contacted"));

    options->addOption(DoubleOption("maxD").store(&maxDistance_).defaultValue(1.2)
                               .description("Maxmimum Centroid Distance to output."));

    options->addOption(DoubleOption("stepD").store(&stepDistance_).defaultValue(0.01)
                               .description("Step of Centroid Distance to output."));


    options->addOption(SelectionOption("ions").store(&sel_ions_).required()
                               .description("Groups containing ions to calculate angle and distance from. One atoms per ion."));

    options->addOption(SelectionOption("rings").store(&sel_rings_).required()
                               .description("Groups containing rings to calculate angle and distance from. Three atoms per ion."));


    options->addOption(BooleanOption("inter_molecule").store(&inter_molecule_).defaultValue(false)
                               .description("Whether only contact between different molecule is calculated. MUST SPECIFY num_ring!!!"));

    options->addOption(BooleanOption("intra_molecule").store(&intra_molecule_).defaultValue(false)
                               .description("Whether only contact between the same molecule is calculated. MUST SPECIFY num_ring!!!"));

    options->addOption(IntegerOption("num_ring").store(&ring_number_in_molecule_).defaultValue(1)
                               .description("Number of rings in each molecule. Specify this without enabling inter_molecule is useless."));

    options->addOption(IntegerOption("num_ion").store(&ion_number_in_molecule_).defaultValue(1)
                               .description("Number of ions in each molecule. Specify this without enabling inter_molecule is useless."));

}

void cationPi::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{
    probability_ = new unsigned long[(int)(maxDistance_ / stepDistance_)];

    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        probability_[i] = 0;
    }


    if(sel_rings_.atomCount() % 3)
    {
        cerr << "ERROR: Rinh selection contains atoms not a multiple of three." << endl;
        exit(1);
    }

    ring_has_contact_ = new bool*[sel_ions_.atomCount()];
    for(int i = 0; i < sel_ions_.atomCount(); i++)
    {
        ring_has_contact_[i] = new bool[sel_rings_.atomCount() / 3];
    }

    data_probability_.setDataSetCount(1);
    data_probability_.setColumnCount(0,1);

    if((not inter_molecule_ and not intra_molecule_) and (ring_number_in_molecule_ > 1 or ion_number_in_molecule_ > 1))
    {
        cerr << "ERROR: Specify num_ring without -inter_molecule or -intra_molecule" << endl;
        exit(1);
    }


    if(inter_molecule_)
    {
        cout << "*****************************************" << endl;
        cout << "**        VERY IMPORTANT NOTICE        **" << endl;
        cout << "** You have specified -inter-molecule_ **" << endl;
        cout << "** You should specify -num_ring & ion  **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ring_number_in_molecule_ << " rings!! **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ion_number_in_molecule_ << "  ions!! **" << endl;
        cout << "*****************************************" << endl;
    }

    if(intra_molecule_)
    {
        cout << "*****************************************" << endl;
        cout << "**        VERY IMPORTANT NOTICE        **" << endl;
        cout << "** You have specified -intra-molecule_ **" << endl;
        cout << "** You should specify -num_ring & ion  **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ring_number_in_molecule_ << " rings!! **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ion_number_in_molecule_ << "  ions!! **" << endl;
        cout << "*****************************************" << endl;
    }


    top_   = top.topology();
    atoms_ = top.topology()->atoms;

    nb_.setCutoff(cutoff_);

}

void cationPi::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dhProbability = pdata->dataHandle(data_probability_);
    dhProbability.startFrame(frnr, fr.time);


    // We first clear all true values in ring_has_contact_ parameter.

    for(int i = 0; i < sel_ions_.atomCount(); i++)
    {
        for(int j = 0; j < sel_rings_.atomCount() / 3; j++)
        {
            ring_has_contact_[i][j] = false;
        }
    }

    // Begin analysis

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel_ions_);

    gmx::ArrayRef<float const[3]>::iterator iter_coordinate_sel;
    for(iter_coordinate_sel = sel_rings_.coordinates().begin(); iter_coordinate_sel != sel_rings_.coordinates().end(); iter_coordinate_sel++)
    {
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(*iter_coordinate_sel);
        AnalysisNeighborhoodPair pair;

        while(pairSearch.findNextPair(&pair))
        {
            int atom_index_in_ion = pair.refIndex();
            int atom_index_in_ring = pair.testIndex();
            int ion_index  = atom_index_in_ion;
            int ring_index = floor(atom_index_in_ring / 3.0);
            ring_has_contact_[ion_index][ring_index] = true;
        }
    }


    // We now iterate through contacted rings to calculate their centroid distances and angles.

    for(int ion_index = 0; ion_index < sel_ions_.atomCount(); ion_index++)
    {
        for(int ring_index = 0; ring_index < sel_rings_.atomCount() / 3; ring_index++)
        {
            if(inter_molecule_
               and floor(ring_index / float(ring_number_in_molecule_))
                   == floor(ion_index / float(ion_number_in_molecule_)))
            {
                continue;
            }

            if(intra_molecule_
               and floor(ring_index / float(ring_number_in_molecule_))
                   != floor(ion_index / float(ion_number_in_molecule_)))
            {
                continue;
            }

            coordinate molecule1[3];
            coordinate molecule2;

            for(int i = 0; i < 3; i++)
            {
                molecule1[i].x = sel_rings_.coordinates()[ring_index * 3 + i][0];
                molecule1[i].y = sel_rings_.coordinates()[ring_index * 3 + i][1];
                molecule1[i].z = sel_rings_.coordinates()[ring_index * 3 + i][2];
            }

            molecule2.x = sel_ions_.coordinates()[ion_index][0];
            molecule2.y = sel_ions_.coordinates()[ion_index][1];
            molecule2.z = sel_ions_.coordinates()[ion_index][2];

            float tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2],
                                             molecule2);


            if (tempDistance <= maxDistance_ )
            {
                probability_[(int) floor(tempDistance / stepDistance_)] += 1;
            }
        }
    }

    dhProbability.finishFrame();
}

void cationPi::finishAnalysis(int /*nframes*/) {}

void cationPi::writeOutput()
{

    // Construct matrix probability
    real *matProbability = new real[(int)(maxDistance_ / stepDistance_)];

    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matProbability[i] = probability_[i] / (real)data_probability_.frameCount();
    }


    if (!fnPropability_.empty())
    {
        if (fnPropability_.compare(".dat"))    {}
        else                                        {fnPropability_ += ".dat";}

        FILE *fpProbability;
        fpProbability = fopen(fnPropability_.c_str(), "w" );

        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            fprintf(fpProbability, "%f\t%f ", i * stepDistance_, matProbability[i]);
            fprintf(fpProbability, "\n");
        }
        fclose(fpProbability);

    }

}




