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
double absorption(const double&lambda){
  // I discussed this coefficient with Paweł Moskal in Sept-Oct 2018
  // Please ask him if it is needed. Maybe you should just update table data
  return polyester_absorp(lambda)*1.8;
}
const std::shared_ptr<Scintillator> create_scintillator(const double&L){
  auto res=MakeScintillator(
    {make_pair(-L/2,L/2),sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5),
    make_shared<DistribTable>(BC420_lambda.clone()),absorption
  );
  res->Configure(Scintillator::Options(1,50));//1 threads, max 50 reflections
  return res;
};

//Where the gamma is registered
const RandomUniform<> x_ph(sizeX.first,sizeX.second);
const RandomUniform<> y_ph(sizeY.first,sizeY.second);


int main(){
  //list of scintillator lengths
  list<double> lengths={30, 50, 75, 100, 150, 200, 300, 400, 500, 1000, 1400, 2000};

  //empty curve to plot
  SortedPoints<> curve1;

  //Varying scintillator length
  for(auto L:lengths){
    {
      cout<<L<<"mm"<<endl;
      //Scintillator
      auto scin1=create_scintillator(L);

      //The block gathering statistics of the time difference between two signals
      auto time_difference1=make_shared<SignalStatictics>();
      {
        //Creating photosensor covering all buttofthe scintillator
        // the parameters: covered surface, the efficiency of the optical glue, quantum efficiency, time resolution
        auto photosensor=[](){return Photosensor({sizeX,sizeY},1.0,tube_QE.func(),tts);};

        //creating channels to transfer signal from left and right photomultipliers
        auto left=make_shared<Signal>(),right=make_shared<Signal>();

        //glue photosensors to the scintillator butts and connect them to the signal channels
        // TimeSignal({make_pair(0,1)}) : the time of signal is obtained for photon order statistics 0 with weight 1
        // To calculate signal time using times of several photons: TimeSignal({make_pair(0,0.5),make_pair(1,0.5)}) will use average of first two photons times
        scin1->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>left));
        scin1->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>right));

        //Invert time of signal from the right photomultiplier
        auto inv_right=SignalInvert();
        right>>inv_right;

        //Add left signal with the inverted one and connect the sum to the block gathering statistics
        (make_shared<SignalSumm>()<<left<<inv_right)>>time_difference1;
      }

      //repeat simulating registering gammas
      for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
        scin1->RegisterGamma({0.0,x_ph(),y_ph()},N_photons);
      }

      // add the statistics for this length tothe curve
      curve1 << make_point( L, (time_difference1->data() + DOI).uncertainty() );
    }
  }

  //Plotting the curve using gnuplot
  Plot("time resolution").Line(curve1,"1. Phm.tube + absorption + DOI")<<"set key on";
  return 0;
}
