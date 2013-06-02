// PO_Blur.C

// Copyright (c) 2009 The Foundry Visionmongers Ltd.  All Rights Reserved.
// Permission is granted to reuse portions or all of this code for the
// purpose of implementing Nuke plugins, or to demonstrate or document
// the methods needed to implemente Nuke plugins.

// This is the name that Nuke will use to store this operator in the
// scripts. So that nuke can locate the plugin, this must also be the
// name of the compiled plugin (with .machine added to the end):
static const char* const CLASS = "PO_Blur";

// This is the text displayed by the [?] button in the control panel:
static const char* const HELP =
  "This is a demonstration of a Nuke plugin that moves pixels by "
  "use of Tile. In this case blocks of pixels are averaged together "
  "to produce the result. Notice that this implementation is much "
  "slower than necessary as it does not reuse calculations done for "
  "adjacent lines."
  "\n\n"
  "The source code has a lot of comments "
  "inserted into it to demonstrate how to write a plugin.";

////////////////////////////////////////////////////////////////

// Most plugins will need these include files. The include files are in
#include "DDImage/Iop.h"
#include "DDImage/Row.h"
#include "DDImage/Tile.h"
#include "DDImage/Knobs.h"
#include "DDImage/Thread.h"
using namespace DD::Image;

//PO inlcudes
#include <TruncPoly/TruncPolySystem.hh>
#include <OpticalElements/OpticalMaterial.hh>
#include <OpticalElements/Spherical5.hh>
#include <OpticalElements/Cylindrical5.hh>
#include <OpticalElements/Propagation5.hh>
#include <OpticalElements/TwoPlane5.hh>
#include <OpticalElements/FindFocus.hh>

#define cimg_display 0
#include <include/CImg.h>
#include <include/spectrum.h>
using namespace cimg_library;

#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;

// Disable some warnings on Microsoft VC++ compilers.
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4305)
#endif

// Define the C++ class that is the new operator. This may be a subclass
// of Iop, or of some subclass of Iop such as Blur:

class PO_Blur : public Iop
{
  // These are the locations the user interface will store into:
  double _knob_sample_mul;
  bool _knob_fixed_samples;
  int _knob_num_lambdas;
  double _knob_focus_point;
  double _knob_obj_distance;
  double _knob_aperture;
  bool _knob_useDepthChannel;
  bool _knob_useApertureMask;

  bool _firstTime; 
  Lock _lock;
  CImg <float> * _img_out;
  CImg <float> * _img_in;
  float _magnification;
  float _sensor_width;
  float _pixel_size;

public: 
  // You must implement these functions:
  PO_Blur(Node*);
  void _validate(bool);
  void _open();
  void _request(int, int, int, int, ChannelMask, int);
  void engine(int y, int x, int r, ChannelMask, Row &);

  virtual void knobs(Knob_Callback);
  virtual int knob_changed(Knob*);

  const char* Class() const { return CLASS; }
  const char* node_help() const { return HELP; } 

  static const Description d;
  // You may need to implement the destructor if you allocate memory:
  //~PO_Blur();

private:
  System44f get_system(float, float,  int = 3);
  bool traceRow(System44f*);
};

// The constructor must initialize the user controls to their default
// values:
PO_Blur::PO_Blur(Node* node) : Iop(node)
{
  _knob_sample_mul = 100.0f; 
  _knob_fixed_samples = false;
  _knob_num_lambdas = 1;
  _knob_focus_point = 5000;
  _knob_obj_distance = 5000;  
  _knob_aperture = 5.6f;
  _knob_useDepthChannel = false;
  _knob_useApertureMask = false;

  _firstTime = true;
  _img_in = NULL;
  _img_out = NULL;
}


// The knobs function describes each user control so the user interface
// can be built. The possible controls are listed in DDImage/Knobs.h:
void PO_Blur::knobs(Knob_Callback f)
{
  // WH_knob(f, &width, "size");
  Int_knob(f, &_knob_num_lambdas, "lambdas", "Lambdas");
  Float_knob(f, &_knob_sample_mul, "samples", "Samples"); 
  SetRange(f, 10, 5000);
  Tooltip(f, "Number of samples to shoot for an input pixel of value 1.");
  Newline(f);
  Bool_knob(f,&_knob_fixed_samples, "fixed_samplesequence", "Fixed Samplesequence"); 
  Tooltip(f, "If enabled uses a fixed samplesequence for all pixel.\nWill result in less noise but may introduce visible sampling patterns.");
  Divider(f);

  Float_knob(f, &_knob_focus_point, "focuspoint", "Focuspoint");
  SetRange(f, 1, 5000);
  Float_knob(f, &_knob_obj_distance, "objectdistance", "Objectdistance");
  Bool_knob(f,&_knob_useDepthChannel, "use_depth_channel", "Use Depth Channel"); 
  Newline(f);
  Float_knob(f, &_knob_aperture, "aperture", "Aperture");
  SetRange(f, 1, 22);
  Bool_knob(f,&_knob_useApertureMask, "use_aperture_mask", "Use Aperture Mask"); 
  Newline(f);
}

// Called whenever a knob is changed
int PO_Blur::knob_changed(Knob* k)
{
  //called when the properties tab is drawn for the fist time
  if(k == &Knob::showPanel) {
     knob("objectdistance")->enable(1 - _knob_useDepthChannel);
     return 1;
  }

  if(k->is("use_depth_channel")) {
     knob("objectdistance")->enable(1 - _knob_useDepthChannel);
     return 1;
  }
  // call parent class
  return Iop::knob_changed(k);
}

// This is a function that creates an instance of the operator, and is
// needed for the Iop::Description to work:
static Iop* PO_Blur_c(Node* node) { return new PO_Blur(node); }

// The Iop::Description is how Nuke knows what the name of the operator is,
// how to create one, and the menu item to show the user. The menu item
// is ignored in recent versions of Nuke, but you might want to set it anyway.
const Iop::Description PO_Blur::d(CLASS, "Filter/PO_Blur", PO_Blur_c);


void PO_Blur::_open()
{
  _firstTime = true;
}

// When the operator runs, "open" is called first. This function
// copies the "info" data from the input operator(s), modifies it
// according to what this operator does, and saves this information so
// that operators connected to this one can see it. The boolean flag
// is so that Nuke can call a "fast" and "slow" version. The fast
// version should avoid opening any files and just make the best
// guess. The slow (argument==true) version will always be called
// before request() or engine() are called:
void PO_Blur::_validate(bool for_real)
{
  copy_info(); // copy bbox channels etc from input0, which will validate it.
}

// After open is done, "request" is called. This is passed a "viewport"
// which describes which channels and a rectangular "bounding box" that will
// be requested. The operator must translate this to the area that will
// be requested from it's inputs and call request() on them. If you don't
// implement this then a default version requests the entire area of
// pixels and channels available from the input. If possible you should
// override this so that the resulting caches are smaller, this also helps
// Nuke display what portions of the tree are being used.
void PO_Blur::_request(int x, int y, int r, int t, ChannelMask channels, int count)
{
  // request all input input as we are going to search the whole input area
  ChannelSet readChannels = input0().info().channels();
  input(0)->request( readChannels, count );
}

// This is the operator that does all the work. For each line in the area
// passed to request(), this will be called. It must calculate the image data
// for a region at vertical position y, and between horizontal positions
// x and r, and write it to the passed row structure. Usually this works
// by asking the input for data, and modifying it:
void PO_Blur::engine ( int y, int l, int r, ChannelMask channels, Row& row )
{
  //cout << "Line " << y << endl;
  {
    Guard guard(_lock);
    if ( _firstTime ) {
      cout << "First time " << endl;
      
      ////////////////////////////////////////
      // the nuke stuff 
      Format format = input0().format();
      // these useful format variables are used later
      const int fx = format.x();
      const int fy = format.y();
      const int fr = format.r();
      const int ft = format.t();
      ChannelSet readChannels = input0().info().channels();
      Interest interest( input0(), fx, fy, fr, ft, readChannels, true );
      interest.unlock();

      ////////////////////////////////////////
      // the PO stuff
      const int height = ft - fy ;
      const int width = fr - fx ;

      // sone predefined constants for now
      const int degree = 3;
      const int filter_size = 1;
      _sensor_width = 36;
      const float lambda_from = 440;
      const float lambda_to = 660;      

      //the output is stored here
      //first delete the old image
      delete _img_out;      
      _img_out = new CImg<float>(width, height, 1, 3, 0);
       
      // Focus on 550nm
      System44f lens = get_system(550, _knob_obj_distance, degree);
      System44f tempsystem = get_system(550, _knob_focus_point, degree);

      // Determine back focal length from degree-1 terms (matrix optics)
      float d3 = find_focus_X(tempsystem);
      cout << "Focus: " << d3 << endl;
      // Compute magnification and output equation system
      _magnification = get_magnification_X(lens >> propagate_5(d3));
      
      // Support of an input image pixel in world plane
      _pixel_size = _sensor_width/(float)width/_magnification;
      cout << "Magnification: " << _magnification << endl;

      // Add the depth to the filmgate propagation
      Transform4f prop = propagate_5(d3, degree);
      lens = lens >> prop;    

      //int samplesFired = 0;
      
      //the framebuffer we save the input in
      _img_in = new CImg<float>(width, height, 1, 3, 0);

      const float sensor_scaling = width / _sensor_width;


      //generate a sample sequence
      // float* x_ap = new float[_knob_sample_mul * 3];
      // float* y_ap = new float[_knob_sample_mul * 3];
      // for (int sample = 0; sample < _knob_sample_mul * 3; ++sample) {
      //   // Rejection-sample points from lens _knob_aperture:
      //   // float x_ap, y_ap;
      //   do {
      //     x_ap[sample] = (rand() / (float)RAND_MAX - 0.5) * 2 * _knob_aperture;
      //     y_ap[sample] = (rand() / (float)RAND_MAX - 0.5) * 2 * _knob_aperture;
      //   } while (x_ap[sample] * x_ap[sample] + y_ap[sample] * y_ap[sample] > _knob_aperture * _knob_aperture);
      // }

      // fetch each row
      for ( int ry = fy; ry < ft; ry++) {
        const int j = ry-fy;

        ////////////////////////////////////////
        // the nuke stuff 
        progressFraction( ry, ft - fy );
        Row row( fx, fr );
        row.get( input0(), ry, fx, fr, readChannels );
        if ( aborted() )
          return;

        ChannelSet done;
        foreach( z, channels )  {
          // skip if we did it as a side-effect of another channel or this is not an rgb channel
          if ( (done & z) || (z>3) || (z<1)    )      
            continue;             
          // Find the rgb channels that belongs to the set this channel is in.
          // Add them all to "done" so we don't run them a second time:
          Channel rchan = brother(z, 0);
          done += rchan;
          Channel gchan = brother(z, 1);
          done += gchan;
          Channel bchan = brother(z, 2);
          done += bchan;

          const float* rIn = row[rchan];
          const float* gIn = row[gchan];
          const float* bIn = row[bchan];
       
          ////////////////////////////////////////
          // the PO stuff
          const float y_sensor = ((j - height/2)/(float)width) * _sensor_width;
          const float y_world = y_sensor / _magnification;
          
#pragma omp parallel for
        // loop through x
          for (int x = fx ; x < fr ; x++)
          {
            const int i = x - fx;

            const float x_sensor = (i / (float)width - 0.5) * _sensor_width;
            const float x_world = x_sensor / _magnification;          
            
            const float rgbin[3] = { rIn[x], gIn[x], bIn[x] };

            // Quasi-importance sampling: 
            // pick number of samples according to pixel intensity, gamma corrected, different weights for rgb amd max 3
            const float impValue =  min(3.0f, pow ( ( (rgbin[0] * 2 + rgbin[1] * 3 + rgbin[2])/ 6), 0.45f) );
            const int num_samples = max(1, (int)(impValue * _knob_sample_mul));
            const float sample_weight = 1.0f / num_samples;


            // With that, we can now start sampling the _knob_aperture:
            for (int sample = 0; sample < num_samples; ++sample) {
              // Rejection-sample points from lens _knob_aperture:
              float x_ap, y_ap;
              do {
                x_ap = (rand() / (float)RAND_MAX - 0.5) * 2 * _knob_aperture;
                y_ap = (rand() / (float)RAND_MAX - 0.5) * 2 * _knob_aperture;
              } while (x_ap * x_ap + y_ap * y_ap > _knob_aperture * _knob_aperture);
            
              float in[5], out[4];

              // Fill in variables and evaluate systems:
              in[0] = x_world;// + _pixel_size * (rand()/(float)RAND_MAX - 0.5);
              in[1] = y_world;
              in[2] = x_ap;
              in[3] = y_ap;

              lens.evaluate(in,out); 
              //samplesFired++;
              // Scale to pixel size:
              out[0] = out[0] * sensor_scaling + width/2;
              out[1] = out[1] * sensor_scaling + height/2;

              // out[2] contains one minus square of Lambertian cosine
              float lambert = sqrt(1 - out[2]);
              if (lambert != lambert) lambert = 0; // NaN check

              _img_out->set_linear_atXY(lambert * sample_weight * rgbin[0], out[0],out[1],0,0, true);
              _img_out->set_linear_atXY(lambert * sample_weight * rgbin[1], out[0],out[1],0,1, true);
              _img_out->set_linear_atXY(lambert * sample_weight * rgbin[2], out[0],out[1],0,2, true);

            }// sample loop            
          }// column x loop          
        }//channel loop       
      }// rowloop


      _firstTime = false;
      //cout << " done " << samplesFired << endl;

    } // end first time
  } // end lock

  // write the output
  Row in( l,r);
  in.get( input0(), y, l, r, channels );
  if ( aborted() )
    return;

  foreach( z, channels ) {
    float *CUR = row.writable(z) + l;
    const float* inptr = in[z] + l;
    const float *END = row[z] + r;
    int x = l;
    while ( CUR < END ) {
      int chanNo = z - 1;
      if (chanNo < 3) { 
        *CUR++ = _img_out->atXY(x++,y,0,chanNo);
      }
      else
      {
        *CUR++ = *inptr++ ;
      }
    }
  } // end loop channel

  //cout << "Done Writing " << endl;
}


////////////////////////////
//The PO Code

System44f PO_Blur::get_system(float lambda, float d0, int degree) {
 // Let's simulate Edmund Optics achromat #NT32-921:
  /* Clear _knob_aperture CA (mm)   39.00
     Eff. Focal Length EFL (mm)   120.00
     Back Focal Length BFL (mm)   111.00
     Center Thickness CT 1 (mm)   9.60
     Center Thickness CT 2 (mm)   4.20
     Radius R1 (mm)   65.22
     Radius R2 (mm)   -62.03
     Radius R3 (mm)   -1240.67
     Substrate  N-SSK8/N-SF10 
  */

  OpticalMaterial glass1("N-SSK8", true);
  OpticalMaterial glass2("N-SF10", true);
  
  // Also try: const float d0 = 5000; // Scene is 5m away
  // const float d0 = 5000000; // Scene is 5km away
  const float R1 = 65.22;
  const float d1 = 9.60;
  const float R2 = -62.03;
  const float d2 = 4.20;
  const float R3 = -1240.67;

  return two_plane_5(d0, degree)
    >> refract_spherical_5(R1,1.f,glass1.get_index(lambda), degree)
    >> propagate_5(d1, degree)
    >> refract_spherical_5(R2,glass1.get_index(lambda),glass2.get_index(lambda), degree)
    >> propagate_5(d2, degree)
    >> refract_spherical_5(R3,glass2.get_index(lambda), 1.f, degree);
}

