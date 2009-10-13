/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:52 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#include "AcceleratorModel/Supports/SupportStructure.h"

SupportStructure::SupportStructure (const string& id, Type type)
        : SequenceFrame(id)
{
    sup1 = new AcceleratorSupport();
    sup2 = (type==girder) ? new AcceleratorSupport() : 0;
}

SupportStructure::SupportStructure (const SupportStructure& rhs)
        : SequenceFrame(rhs)
{
    sup1 = new AcceleratorSupport();
    if(rhs.sup2)
        sup2 = new AcceleratorSupport();
}

SupportStructure::~SupportStructure ()
{
    delete sup1;
    if(sup2)
        delete sup2;
}

int SupportStructure::ExportSupports (AcceleratorSupportList& supports)
{
    supports.push_back(sup1);
    if(sup2)
        supports.push_back(sup2);
    return sup2 ? 2:1;
}

Transform3D SupportStructure::GetLocalFrameTransform () const
{
    UpdateSupportTransform();
    return LatticeFrame::GetLocalFrameTransform()*Ts;
}

void SupportStructure::ConsolidateConstruction ()
{
    SequenceFrame::ConsolidateConstruction();

    Transform3D gT = GetPhysicalTransform(GLOBAL_FRAME);
    double s0 = GetPosition(GLOBAL_FRAME);

    if(!sup2) {
        // support is at origin, so this is now trivial
        sup1->SetPosition(s0,gT.X().x,gT.X().z);
        Rg = gT.R();
    }
    else {
        // Here life gets a little more complicated, as we need to transform
        // the origin information to the entrance and exit of the frame,
        // as that's where the supports are located.
        AcceleratorGeometry::Extent ext=GetLocalGeometryExtent();

        Transform3D t1 = GetGeometryTransform(AcceleratorGeometry::entrance)*gT;
        sup1->SetPosition(s0+ext.first,t1.X().x,t1.X().z);
        Rg = t1.R();

        t1 = GetGeometryTransform(AcceleratorGeometry::exit)*gT;
        sup2->SetPosition(s0+ext.second,t1.X().x,t1.X().z);
    }
}

void SupportStructure::UpdateSupportTransform () const
{

    // Check if cached state needs updating, otherwise
    // return.

    if(!sup1->modified || (sup2 && !sup2->modified))
        return;

    // The local frame transformation t0 must now be modified
    // by any 'ground motion' like transformation from the
    // AcceleratorSupports.

    if(!sup2) {
        // This is the trivial case, with the translation applied
        // to the centre of the frame
        Vector3D X = Rg(sup1->GetOffset());
        Ts = Transform3D::translation(X.x,X.y,X.z);
    }
    else {
        // This is much more complicated, expecially as
        // large offset of the two supports will effectively
        // stretch the geometry. We will assume that the offsets
        // are small, and that the stretching is negligble.

        Vector3D X1 = Rg(sup1->GetOffset());
        Vector3D X2 = Rg(sup2->GetOffset());

        // check to see if we have a simple translation
        if(X1==X2)
            Ts = Transform3D::translation(X1.x,X1.y,X1.z);
        else {

            // Approximate rotations about x- and y-axis.
            Vector3D dX=X2-X1;
            double s0 = sup1->DistanceTo(*sup2);
            double phix = -dX.y/s0;
            double phiy =  dX.x/s0;

            // Calculate the approximate transformation caused by the two offsets
            // about the entrance plane reference frame.

            Rotation3D R = Rotation3D::rotationX(phix)*Rotation3D::rotationY(phiy);
            Transform3D T1 = Transform3D(Point3D(X1.x,X1.y,X1.z),R);

            // Convert Ts to a transformation about the local frame origin.
            Transform3D Tin = GetGeometryTransform(AcceleratorGeometry::entrance);
            Ts = Tin.inv()*T1*Tin;
        }
    }

    sup1->modified = false;
    if(sup2)
        sup2->modified = false;
}

// Class GirderMount


const string& GirderMount::GetType () const
{
    _TYPESTR(GirderMount);
}

ModelElement* GirderMount::Copy () const
{
    return new GirderMount(*this);
}

// Class SimpleMount


const string& SimpleMount::GetType () const
{
    _TYPESTR(SimpleMount);
}

ModelElement* SimpleMount::Copy () const
{
    return new SimpleMount(*this);
}

