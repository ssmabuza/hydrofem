
#include "Hydrofem_Mesh_gtests.hpp"
#include "Hydrofem_Linear_gtests.hpp"
#include "Hydrofem_Disc_FE_gtests.hpp"


int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc,argv);
  return RUN_ALL_TESTS();
}