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
const double tts=0.128;//ns
// I discussed DOI with Paweł Moskal two years ago
const auto DOI = value<>(0.0,0.063);//ns

//Scintillator
const double opt_dens=1.58;
const auto sizeX=make_pair(0.0,12.0);//milimeters
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
            make_shared<DistribTable>(BC420_lambda.clone()),absorption
  );
  scin1->Configure(Scintillator::Options(1,50));//1 thread, max 50 reflections

  vector<shared_ptr<SignalStatictics>> time_differences;
  {
        auto left_sorter = make_shared<SignalSort>(), right_sorter = make_shared<SignalSort>();
        for(size_t i=0; i<5; i++)for(size_t j=0; j<2; j++){
            auto photosensor=[i,j](){return Photosensor(
                            {make_pair( 6.0*j, 6.0*(j+1) ),make_pair( 6.0*i, 6.0*(i+1) )},
                            1.0,Si_Photo_QE.func(),tts);};
            auto left_sensor=photosensor(), right_sensor=photosensor();

            scin1->Surface(0,RectDimensions::Left) >> left_sensor;
            scin1->Surface(0,RectDimensions::Right) >> right_sensor;

            auto left_wire=make_shared<Signal>(), right_wire=make_shared<Signal>();

            left_sensor >> ( TimeSignal({make_pair(0,1)}) >> left_wire );
            right_sensor >> ( TimeSignal({make_pair(0,1)}) >> right_wire );

            left_sorter << left_wire;
            right_sorter << right_wire;
        }
        for(size_t index=0; index<(2*5); index++){
            auto left_wire=make_shared<Signal>(), right_wire=make_shared<Signal>();

            left_sorter >> left_wire;
            right_sorter >> (SignalInvert() >> right_wire);

            auto time_difference = make_shared<SignalStatictics>();
            ( make_shared<SignalSumm>()<<left_wire<<right_wire ) >> time_difference;
            time_differences.push_back(time_difference);
        }
   }

   for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
        scin1->RegisterGamma( {z_ph,x_ph(),y_ph()}, N_photons);
   }

   SortedPoints<> res_by_order_statistics, efficiency_by_order_statistics;
   for(size_t i=0,n=time_differences.size(); i<n; i++){
       if (time_differences[i]->data().Sample().count() > 1) // means that we can obtain resolution
            res_by_order_statistics << make_point( (double)(i+1), (time_differences[i]->data() + DOI).uncertainty() );
   }

   Plot("matrix_resolution").Line(res_by_order_statistics)<<"set key on"
                                                       << "set xlabel 'photosensor order statistics'"
                                                       << "set ylabel 'time resolution'" ;
  return 0;
}
