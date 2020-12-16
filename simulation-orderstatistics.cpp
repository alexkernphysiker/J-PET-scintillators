// this file is distributed under
// MIT license
#include <iostream>
#include <list>
#include <math_h/tabledata.h>
#include <math_h/randomfunc.h>
#include <gnuplot_wrap.h>
#include <RectScin/scintillator.h>
#include <RectScin/sensitive.h>
#include <RectScin/photon2signal.h>
#include <RectScin/signal_processing.h>
#include <RectScin/signal_statistics.h>
#include "model.h"

using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
using namespace RectangularScintillator;

//Time resolution
const double tts=0.068;//ns
const auto DOI = value<>(0.0,0.063);//ns

//Scintillator
const double opt_dens=1.58;
const auto sizeX=make_pair(0.0,6.0);//milimeters
const auto sizeY=make_pair(0.0,30.0);//milimeters
const auto sizeZ=make_pair(-150.0,150.0);//milimeters
double absorption(const double&lambda){
  // I discussed this coefficient with Paweł Moskal two years ago
  // Please ask him if it is needed. Maybe you should just update table data
  return polyester_absorp(lambda)*1.8;
}

//Where the gamma is registered
const RandomUniform<> x_ph(sizeX.first,sizeX.second);
const RandomUniform<> y_ph(sizeY.first,sizeY.second);
const double z_ph = 0.0;


int main(){

  //Scintillator
  auto scin1=MakeScintillator(
            {sizeZ,sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5),
            make_shared<DistribTable>(LinearInterpolation(BC420_lambda.clone())),absorption
  );
  scin1->Configure(Scintillator::Options(1,50));//1 thread, max 50 reflections

  //The block gathering statistics of the time difference between two signals
  vector<shared_ptr<SignalStatictics>> time_differences;
      {
        //Creating photosensor covering all butt of the scintillator
        // the parameters: covered surface, the efficiency of the optical glue, quantum efficiency, time resolution
        auto photosensor=[](){return Photosensor({sizeX,sizeY},1.0,tube_QE.func(),tts);};
        //glue photosensors to the scintillator butts
        auto left_sensor=photosensor(), right_sensor=photosensor();
        scin1->Surface(0,RectDimensions::Left) >> left_sensor;
        scin1->Surface(0,RectDimensions::Right) >> right_sensor;

        for(size_t index=0; index<50; index++){//cycle over different order statistics
            // virtual wires to conduct the "signal" corresponding to needed poton's registration time
            auto left_wire=make_shared<Signal>(), right_wire=make_shared<Signal>();

            // TimeSignal({make_pair(index,1)}) : the time of signal is obtained for the order statistics of index
            left_sensor >> ( TimeSignal({make_pair(index,1)}) >> left_wire );
            right_sensor >> ( TimeSignal({make_pair(index,1)}) >> right_wire );

            //Invert time of signal from the right photomultiplier
            auto inv_right=SignalInvert();
            right_wire >> inv_right;

            //Add left signal with the inverted right signal and connect the sum to the block gathering statistics
            auto time_difference = make_shared<SignalStatictics>();
            (make_shared<SignalSumm>()<<left_wire<<inv_right)>>time_difference;

            //pushing this statistics to vector
            time_differences.push_back(time_difference);
        }

      }

   //repeat simulating registering gammas
   for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
        scin1->RegisterGamma( {z_ph,x_ph(),y_ph()}, N_photons);
   }

   //empty curve to plot
   SortedPoints<> res_by_order_statistics;
   for(size_t i=0,n=time_differences.size(); i<n; i++){
       res_by_order_statistics << make_point( (double)i, (time_differences[i]->data() + DOI).uncertainty() );
   }

  //Plotting the curve using gnuplot
  Plot("time resolution").Line(res_by_order_statistics,"1. Phm.tube + absorption + DOI")<<"set key on";
  return 0;
}
