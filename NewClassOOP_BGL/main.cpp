#include "boltzmann_alpha0.h"

void make_dir()
{
  char subfile[300];
  sprintf(subfile, "mkdir %s", Subdirectory);
  system(subfile);
}

int main(int argc, char const *argv[])
{
  make_dir();
  boltzmann_alpha0 BGL_Eqs;
  BGL_Eqs.evolution();
  return 0;
}

