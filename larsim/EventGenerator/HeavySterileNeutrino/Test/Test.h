#ifndef STERILE_TEST_H_
#define STERILE_TEST_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

namespace hsngen {

	void StupidFunction();
	void TestFunction(CLHEP::HepRandomEngine& engine);

}

#endif