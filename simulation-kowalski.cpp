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

const double tts_tube=0.068;//ns
const auto DOI = value<>(0.0,0.063);//ns
const double opt_dens=1.58;
const auto sizeX=make_pair(0.0,6.0);//milimeters
const auto sizeY=make_pair(0.0,30.0);//milimeters
const RandomUniform<> x_ph(sizeX.first,sizeX.second);
const RandomUniform<> y_ph(sizeY.first,sizeY.second);
double absorption(const double&lambda){
  return polyester_absorp(lambda)*1.8;
}
const std::shared_ptr<Scintillator> withabsorption(const double&l){
  auto res=MakeScintillator(
    {make_pair(-l/2,l/2),sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5),
    make_shared<DistribTable>(BC420_lambda),absorption
  );
  res->Configure(Scintillator::Options(4,50));//4 threads, max 50 reflections
  return res;
};
int main(int , char **){
  list<double> lengths={30, 50, 75, 100, 150, 200, 300, 400, 500, 1000, 1400, 2000};
  SortedPoints<> curve1;
  for(auto L:lengths){
    {//photomultipliers
      cout<<L<<"mm"<<endl;
      auto scin1=withabsorption(L);
      auto time_difference1=make_shared<SignalStatictics>();
      {
        auto photosensor=[](){return Photosensor({sizeX,sizeY},1.0,tube_QE,tts_tube);};
        auto left=make_shared<Signal>(),right=make_shared<Signal>();
        //TimeSignal({make_pair(0,1)}) : photon order statistics 0, weight 1
        // To calculate signal time using times of several photons:
        // TimeSignal({make_pair(0,0.5),make_pair(1,0.5)})
        // will use average of first two photons times
        // (Indeed, order statistics is not measured)
        scin1->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>left));
        scin1->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>right));
        auto inv_right=SignalInvert();
        right>>inv_right;
        (make_shared<SignalSumm>()<<left<<inv_right)>>time_difference1;
      }
      for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
        scin1->RegisterGamma({0.0,x_ph(),y_ph()},N_photons);
      }
      curve1<<make_point(L,(time_difference1->data() + DOI).uncertainty());
    }
  }
  Plot("kowalski")
  .Line(curve1,"1. Phm.tube + absorption + DOI","kowalski")
  <<"set key on";
  return 0;
}
