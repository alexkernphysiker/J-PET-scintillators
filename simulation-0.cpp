// this file is distributed under 
// MIT license
#include <iostream>
#include <list>
#include <math_h/tabledata.h>
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
const double tts=0.128;//ns
const auto DOI = value<>(0.0,0.063);//ns
const double opt_dens=1.58;
const auto sizeX=make_pair(0.0,7.0);//milimeters
const auto sizeY=make_pair(0.0,19.0);//milimeters
const double x_ph=3.5;//milimeters
const double y_ph=9.5;//milimeters
double absorption(const double&lambda){
  return polyester_absorp(lambda)*1.8;
}
const std::shared_ptr<Scintillator> withabsorption(const double&l){
    return MakeScintillator(
      {make_pair(-l/2,l/2),sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5),
      make_shared<DistribTable>(BC420_lambda),absorption
    );
};
const vector<vector<pair<double,double>>> si_phm_matrix={
  {make_pair(0.0,3.0),make_pair(0.0,3.0)},
  {make_pair(0.0,3.0),make_pair(4.0,7.0)},
  {make_pair(0.0,3.0),make_pair(8.0,11.0)},
  {make_pair(0.0,3.0),make_pair(12.0,15.0)},
  {make_pair(0.0,3.0),make_pair(16.0,19.0)},
  {make_pair(4.0,7.0),make_pair(0.0,3.0)},
  {make_pair(4.0,7.0),make_pair(4.0,7.0)},
  {make_pair(4.0,7.0),make_pair(8.0,11.0)},
  {make_pair(4.0,7.0),make_pair(12.0,15.0)},
  {make_pair(4.0,7.0),make_pair(16.0,19.0)},
};
int main(int , char **){
    list<double> lengths={30, 50, 75, 100, 150, 200, 300, 400, 500, 1000, 1400, 2000};
    SortedPoints<> curve1,curve2,curve3;
    for(auto L:lengths){
        {//photomultipliers
          cout<<L<<"mm"<<endl;
          cout<<"Creating system 1"<<endl;
          auto scin1=withabsorption(L);
          auto time_difference1=make_shared<SignalStatictics>();
          {
                auto photosensor=[](){return Photosensor({sizeX,sizeY},1.0,tube_QE,tts_tube);};
                auto left=make_shared<Signal>(),right=make_shared<Signal>();
                scin1->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>left));
                scin1->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>right));
                auto inv_right=SignalInvert();
                right>>inv_right;
                (make_shared<SignalSumm>()<<left<<inv_right)>>time_difference1;
          }
          cout<<"Creating system 2"<<endl;
          auto scin2=withabsorption(L);
          auto time_difference2=make_shared<SignalStatictics>();
          {
            auto left=make_shared<SignalSortAndSelect>(2),right=make_shared<SignalSortAndSelect>(2);
            for(const auto& size_m:si_phm_matrix){
                auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE,tts);};
                auto l=make_shared<Signal>(),r=make_shared<Signal>();
                scin2->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
                scin2->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
                left<<l;right<<r;
            }
            auto inv_right=SignalInvert();
            right>>inv_right;
            (make_shared<SignalSumm>()<<left<<inv_right)>>time_difference2;
          }
          cout<<"Creating system 3-a"<<endl;
          auto scin3=withabsorption(L);
          vector<shared_ptr<SignalStatictics>> time_differences3={};
          {
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<si_phm_matrix.size();i++){
              time_differences3.push_back(make_shared<SignalStatictics>());
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              left>>l;right>>r;
              auto inv_r=SignalInvert();
              r>>inv_r;
              (make_shared<SignalSumm>()<<l<<inv_r)>>time_differences3[i];
            }
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE,tts);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin3->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin3->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }
          cout<<"Run 1"<<endl;
          for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
                scin1->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin2->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin3->RegisterGamma({0.0,x_ph,y_ph},N_photons);
          }
          cout<<"Getting points 1,2"<<endl;
          curve1<<make_point(L,(time_difference1->data()+DOI).uncertainty());
          curve2<<make_point(L,(time_difference2->data()+DOI).uncertainty());
          cout<<"Creating system 3-b"<<endl;
          auto scin3_final=withabsorption(L);
          auto time_difference3=make_shared<SignalStatictics>();
          {
            double norm=0;
            auto Sl=make_shared<SignalSumm>(),Sr=make_shared<SignalSumm>();
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<time_differences3.size();i++){
              time_differences3.push_back(make_shared<SignalStatictics>());
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              if(time_differences3[i]->data().Sample().count()>2){
                double w=pow(time_differences3[i]->data().uncertainty(),2);
                norm+=w;
                left >>(SignalMultiply(time_differences3[i]->data().Sample().count())>>l);
                right>>(SignalMultiply(time_differences3[i]->data().Sample().count())>>r);
                Sl<<l;Sr<<r;
              }
            }
            auto inv_r=SignalInvert();
            Sr>>inv_r;
            ((make_shared<SignalSumm>()<<Sl<<inv_r)>>SignalMultiply(1.0/norm))>>time_difference3;
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE,tts);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin3_final->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin3_final->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }
          cout<<"Run 2"<<endl;
          for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
                scin3_final->RegisterGamma({0.0,x_ph,y_ph},N_photons);
          }
          cout<<"Getting point 3"<<endl;
          curve3<<make_point(L,time_difference3->data().uncertainty());
        }
    }
    Plot("0")
    .Line(curve1,"tube + absorption + DOI","0-1")
    .Line(curve2,"2x5 matrix (3rd) + absorption + DOI","0-2")
    .Line(curve3,"2x5matrix(1-3 wieghted) + absorption","0-3")
    <<"set key on";
    return 0;
}

