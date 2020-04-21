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
const auto DOI = value<>(0.0,0.063);//ns
const double opt_dens=1.58;
const auto sizeX=make_pair(0.0,7.0);//milimeters
const auto sizeY=make_pair(0.0,19.0);//milimeters
const double x_ph=3.5;//milimeters
const double y_ph=9.5;//milimeters
const std::shared_ptr<Scintillator> absorptionless(const double&l){
    auto res=MakeScintillator_absorptionless(
      {make_pair(-l/2,l/2),sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5)
      );
    res->Configure(Scintillator::Options(4,50));//2 threads, max 50 reflections
    return res;
};
double absorption(const double&lambda){
  return polyester_absorp(lambda)*0.552;
}
const std::shared_ptr<Scintillator> withabsorption(const double&l){
    auto res=MakeScintillator(
      {make_pair(-l/2,l/2),sizeX,sizeY},opt_dens,TimeDistribution2(0.005,0.2,1.5),
      make_shared<DistribTable>(LinearInterpolation(BC420_lambda.clone())),absorption
    );
    res->Configure(Scintillator::Options(4,50));//4 threads, max 50 reflections
    return res;
};
const vector<vector<pair<double,double>>> si_phm_matrix={
  {make_pair(0.0,3.5),make_pair(0.0,3.5)},
  {make_pair(0.0,3.5),make_pair(3.5,7.5)},
  {make_pair(0.0,3.5),make_pair(7.5,11.5)},
  {make_pair(0.0,3.5),make_pair(11.5,15.5)},
  {make_pair(0.0,3.5),make_pair(15.5,19.0)},
  {make_pair(3.5,7.0),make_pair(0.0,3.5)},
  {make_pair(3.5,7.0),make_pair(3.5,7.5)},
  {make_pair(3.5,7.0),make_pair(7.5,11.5)},
  {make_pair(3.5,7.0),make_pair(11.5,15.5)},
  {make_pair(3.5,7.0),make_pair(15.5,19.0)},
};
int main(int , char **){
    list<double> lengths={30, 50, 75, 100, 150, 200, 300, 400, 500, 1000, 1400, 2000};
    SortedPoints<> curve1,curve2,curve3,curve4,curve5;
    for(auto L:lengths){
        {//photomultipliers
          cout<<L<<"mm"<<endl;
          auto scin1=absorptionless(L);
          vector<shared_ptr<SignalStatictics>> time_differences1={};
          {
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<si_phm_matrix.size();i++){
              time_differences1.push_back(make_shared<SignalStatictics>());
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              left>>l;right>>r;
              auto inv_r=SignalInvert();
              r>>inv_r;
              (make_shared<SignalSumm>()<<l<<inv_r)>>time_differences1[i];
            }
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,[](const double&){return 1.0;},0);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin1->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin1->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

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
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.040);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin3->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin3->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          auto scin4=withabsorption(L);
          vector<shared_ptr<SignalStatictics>> time_differences4={};
          {
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<si_phm_matrix.size();i++){
              time_differences4.push_back(make_shared<SignalStatictics>());
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              left>>l;right>>r;
              auto inv_r=SignalInvert();
              r>>inv_r;
              (make_shared<SignalSumm>()<<l<<inv_r)>>time_differences4[i];
            }
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.080);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin4->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin4->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          auto scin5=withabsorption(L);
          vector<shared_ptr<SignalStatictics>> time_differences5={};
          {
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<si_phm_matrix.size();i++){
              time_differences5.push_back(make_shared<SignalStatictics>());
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              left>>l;right>>r;
              auto inv_r=SignalInvert();
              r>>inv_r;
              (make_shared<SignalSumm>()<<l<<inv_r)>>time_differences5[i];
            }
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.128);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin5->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin5->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          cout<<"Run 1"<<endl;
          for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
                scin1->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin3->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin4->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin5->RegisterGamma({0.0,x_ph,y_ph},N_photons);
          }


          auto scin1_final=withabsorption(L);
          auto time_difference1=make_shared<SignalStatictics>();
          {
            double norm=0;
            auto Sl=make_shared<SignalSumm>(),Sr=make_shared<SignalSumm>();
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<time_differences1.size();i++){
              if(time_differences1[i]->data().Sample().count()>2){
                double w=1.0/pow(time_differences1[i]->data().uncertainty(),2);
                double w0=1.0/pow(time_differences1[0]->data().uncertainty(),2);
                if(w>=w0){
                  norm+=w;
                  auto l=make_shared<Signal>(),r=make_shared<Signal>();
                  left >>(SignalMultiply(w)>>l);
                  right>>(SignalMultiply(w)>>r);
                  Sl<<l;Sr<<r;
                }
              }
            }
            auto inv_r=SignalInvert();
            Sr>>inv_r;
            (make_shared<SignalSumm>()<<Sl<<inv_r)>>(SignalMultiply(1.0/norm)>>time_difference1);
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,[](const double&){return 1.0;},0);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin1_final->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin1_final->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          auto scin3_final=withabsorption(L);
          auto time_difference3=make_shared<SignalStatictics>();
          {
            double norm=0;
            auto Sl=make_shared<SignalSumm>(),Sr=make_shared<SignalSumm>();
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<time_differences3.size();i++){
              if(time_differences3[i]->data().Sample().count()>2){
                double w=1.0/pow(time_differences3[i]->data().uncertainty(),2);
                double w0=1.0/pow(time_differences3[0]->data().uncertainty(),2);
                if(w>=w0){
                  norm+=w;
                  auto l=make_shared<Signal>(),r=make_shared<Signal>();
                  left >>(SignalMultiply(w)>>l);
                  right>>(SignalMultiply(w)>>r);
                  Sl<<l;Sr<<r;
                }
              }
            }
            auto inv_r=SignalInvert();
            Sr>>inv_r;
            (make_shared<SignalSumm>()<<Sl<<inv_r)>>(SignalMultiply(1.0/norm)>>time_difference3);
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.040);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin3_final->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin3_final->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          auto scin4_final=withabsorption(L);
          auto time_difference4=make_shared<SignalStatictics>();
          {
            double norm=0;
            auto Sl=make_shared<SignalSumm>(),Sr=make_shared<SignalSumm>();
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<time_differences4.size();i++){
              if(time_differences4[i]->data().Sample().count()>2){
                double w=1.0/pow(time_differences4[i]->data().uncertainty(),2);
                double w0=1.0/pow(time_differences4[0]->data().uncertainty(),2);
                if(w>=w0){
                  norm+=w;
                  auto l=make_shared<Signal>(),r=make_shared<Signal>();
                  left >>(SignalMultiply(w)>>l);
                  right>>(SignalMultiply(w)>>r);
                  Sl<<l;Sr<<r;
                }
              }
            }
            auto inv_r=SignalInvert();
            Sr>>inv_r;
            (make_shared<SignalSumm>()<<Sl<<inv_r)>>(SignalMultiply(1.0/norm)>>time_difference4);
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.080);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin4_final->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin4_final->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          auto scin5_final=withabsorption(L);
          auto time_difference5=make_shared<SignalStatictics>();
          {
            double norm=0;
            auto Sl=make_shared<SignalSumm>(),Sr=make_shared<SignalSumm>();
            auto left=make_shared<SignalSort>(),right=make_shared<SignalSort>();
            for(int i=0;i<time_differences5.size();i++){
              if(time_differences5[i]->data().Sample().count()>2){
                double w=1.0/pow(time_differences5[i]->data().uncertainty(),2);
                double w0=1.0/pow(time_differences5[0]->data().uncertainty(),2);
                if(w>=w0){
                  norm+=w;
                  auto l=make_shared<Signal>(),r=make_shared<Signal>();
                  left >>(SignalMultiply(w)>>l);
                  right>>(SignalMultiply(w)>>r);
                  Sl<<l;Sr<<r;
                }
              }
            }
            auto inv_r=SignalInvert();
            Sr>>inv_r;
            (make_shared<SignalSumm>()<<Sl<<inv_r)>>(SignalMultiply(1.0/norm)>>time_difference5);
            for(const auto& size_m:si_phm_matrix){
              auto photosensor=[&size_m](){return Photosensor(size_m,1.0,Si_Photo_QE.func(),0.128);};
              auto l=make_shared<Signal>(),r=make_shared<Signal>();
              scin5_final->Surface(0,RectDimensions::Left)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>l));
              scin5_final->Surface(0,RectDimensions::Right)>>(photosensor()>>(TimeSignal({make_pair(0,1)})>>r));
              left<<l;right<<r;
            }
          }

          cout<<"Run 2"<<endl;
          for(unsigned int cnt=0;cnt<virtual_experiments_count;cnt++){
                scin1_final->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin3_final->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin4_final->RegisterGamma({0.0,x_ph,y_ph},N_photons);
                scin5_final->RegisterGamma({0.0,x_ph,y_ph},N_photons);
          }
          cout<<"Getting point 3"<<endl;
          curve1<<make_point(L,time_difference1->data().uncertainty());
          curve2<<make_point(L,time_difference3->data().uncertainty());
          curve3<<make_point(L,(time_difference3->data()+DOI).uncertainty());
          curve4<<make_point(L,(time_difference4->data()+DOI).uncertainty());
          curve5<<make_point(L,(time_difference5->data()+DOI).uncertainty());
        }
    }
    Plot("5")
    .Line(curve1,"1. 2x5 solid, 100% eff, tts=0","5-1")
    .Line(curve2,"2. 2x5 solid, tts=0.040 + absorption","5-2")
    .Line(curve3,"3. 2x5 solid, tts=0.040 + absorption + DOI","5-3")
    .Line(curve4,"4. 2x5 solid, tts=0.080 + absorption + DOI","5-4")
    .Line(curve5,"5. 2x5 solid, tts=0.128 + absorption + DOI","5-5")
    <<"set key on";
    return 0;
}
