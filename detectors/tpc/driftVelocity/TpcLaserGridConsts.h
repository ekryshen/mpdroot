#ifndef TPCLASERGRIDVELOCITYCONSTS_H
#define TPCLASERGRIDVELOCITYCONSTS_H

class TpcLaserGridConsts {
protected:
   const double        z_pos_offset = 30.;
   const double        eps_z        = 8.;
   const int           num_layers   = 4;
   std::vector<double> loc_y_offset = {
      /*10., 20., 35.,*/ 85.}; // last element must exist and should be a bit greater than max local R of sector

   std::vector<double> z_layers; //(num_layers);
};

#endif
