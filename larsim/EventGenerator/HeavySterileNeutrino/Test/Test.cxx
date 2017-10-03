#include "Test.h"

namespace hsngen {

	void TestFunction(CLHEP::HepRandomEngine& engine)
	{
		CLHEP::RandFlat flat(engine);
		printf("Printing number from function engine: %.3f\n", flat());
		return;
	}

	void StupidFunction()
	{
		printf("This is a very stupid function.\n");
	}

}