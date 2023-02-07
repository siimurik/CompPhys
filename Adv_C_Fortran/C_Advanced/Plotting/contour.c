#include <stdio.h>
#include "pbPlots.h"
#include "supportLib.h"

int main(int argc, char** argv) {
  // Create a pbplot object and configure it
  pbplot plot;
  pbplot_init(&plot);
  pbplot_set_title(&plot, "My Contour Plot");
  pbplot_set_xlabel(&plot, "X Axis");
  pbplot_set_ylabel(&plot, "Y Axis");

  // Create a 2D data array for the contour plot
  double data[][3] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };

  // Add the contour plot to the pbplot object
  pbplot_add_contour(&plot, (double*)data, 3, 3, 0, 3, 0, 3);

  // Show the plot
  pbplot_show(&plot);

  // Free the memory used by the pbplot object
  pbplot_free(&plot);

  return 0;
}
