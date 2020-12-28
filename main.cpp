#include <iostream>
#include "cationPi.h"

int main(int argc, char * argv[]) {
    // insert code here...
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<cationPi>(argc, argv);

}

