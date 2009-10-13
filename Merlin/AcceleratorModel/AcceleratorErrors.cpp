//---------------
//  AcceleratorErrors
//  DK 1.12.08

#include "AcceleratorErrors.h"
using namespace std;
namespace {
   class Errors {
   public:
        Errors( double  vvx, double  vvy, double  vvz,
	        double meanx, double meany, double meanz, 
                const string& p, bool clear, bool trans, ostream* l)
              : vx(vvx),vy(vvy),vz(vvz), 
                mx(meanx),my(meany),mz(meanz),
                c(clear),t(trans),pat("*."+p),log(l)  {};

        void operator()(LatticeFrame* frame) const {
            if(frame && pat((*frame).GetQualifiedName())) {

                if(c) frame->ClearLocalFrameTransform();

                double ex= vx==0 ? mx : RandomNG::normal(mx,vx);
                double ey= vy==0 ? my : RandomNG::normal(my,vy);
                double ez= vz==0 ? mz : RandomNG::normal(mz,vz);
                if(t) {
                   frame->Translate(ex,ey,ez);
                   if(log) (*log)<<(*frame).GetQualifiedName()<<" translate: "
                                 <<ex<<" "<<ey<<" "<<ez<<endl;
                } else {
                   if(ex) frame->RotateX(ex);
                   if(ey) frame->RotateY(ey);
                   if(ez) frame->RotateZ(ez);
                   if(log) (*log)<<(*frame).GetQualifiedName()<<" rotate: "
                                 <<ex<<" "<<ey<<" "<<ez<<endl;
                }
            }
        };
   private:
        double vx,vy,vz;
        double mx,my,mz;
        bool c,t;
        StringPattern pat;
        ostream* log;
   };
};
void AcceleratorErrors::ApplyShifts(AcceleratorModel::Beamline& b, const string& p){
             for_each(b.begin(),b.end(),Errors(vx,vy,vz,mx,my,mz,p,clear,true,log));
};
void AcceleratorErrors::ApplyRotations(AcceleratorModel::Beamline& b, const string& p){
             for_each(b.begin(),b.end(),Errors(vx,vy,vz,mx,my,mz,p,clear,false,log));
};
