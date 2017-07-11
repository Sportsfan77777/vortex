#ifndef BOUNDARY_CONDITIONS_H_   /* Include guard */
#define BOUNDARY_CONDITIONS_H_

/// Fixed ///
void FixedBoundary(const Data *d, RBox *box, int side, Grid *grid);

/// Evanescent ///
void EvanescentBoundary(const Data *d, RBox *box, int side, Grid *grid);


#endif // BOUNDARY_CONDITIONS_H_