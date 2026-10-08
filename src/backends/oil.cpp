// Rename only the entry point; preserve the upstream scientific implementation.
#define main siliconelab_oil_main
#include "../../vendor/Silicone_Oil/Generator/oil_generator.cpp"
#undef main
