//
// $Id$
//
/*
    Copyright (C) 2005- $Date$, Hammersmith Imanet Ltd
    Internal GE use only
*/
/*!
  \file
  \ingroup motion_utilities
  \brief Utility to report RMSE w.r.t. reference position and within the frames

  \author Kris Thielemans
  $Date$
  $Revision$
  
  \par Usage
\verbatim
  report_movement \\
     [--frame_num_to_process number]\\
     [par_file]
\endverbatim
  See class documentation for stir::ReportMovement for more info, including the format
  of the par_file. 

  Command line switches override any values in the par_file.

*/

#include "local/stir/motion/TimeFrameMotion.h"
#include "stir/Succeeded.h"
#include "stir/is_null_ptr.h"
#include "stir/CartesianCoordinate3D.h"

START_NAMESPACE_STIR

/*! \ingroup motion
  \brief A class for reporting the movement within the frame w.r.t. to the reference position

  \par Example par file
  \see TimeFrameMotion for other parameters
  \verbatim
  ReportMovement Parameters:=

  ; parameters from TimeFrameMotion

  reference point 1:={50,0,0} 
  reference point 2:={100,0,30} 
  reference point 3:={100,0,-30} 

  END :=
\endverbatim
*/  
class ReportMovement : public TimeFrameMotion
{
private:
  typedef TimeFrameMotion base_type;
public:
  ReportMovement(const char * const par_filename);

  virtual Succeeded process_data();

protected:
  std::vector<CartesianCoordinate3D<float> > reference_points;
  
  //! parsing functions
  virtual void set_defaults();
  virtual void initialise_keymap();
  virtual bool post_processing();

};

void 
ReportMovement::set_defaults()
{
  base_type::set_defaults();
  reference_points.resize(3, make_coordinate(0.F,0.F,0.F));
}

void 
ReportMovement::initialise_keymap()
{
  parser.add_start_key("ReportMovement Parameters");

  base_type::initialise_keymap();
  parser.add_key("reference point 1", &reference_points[0]);
  parser.add_key("reference point 2", &reference_points[1]);
  parser.add_key("reference point 3", &reference_points[2]);

  parser.add_stop_key("END");
}

ReportMovement::
ReportMovement(const char * const par_filename)
{
  set_defaults();
  if (par_filename!=0)
    {
      if (parse(par_filename)==false)
	exit(EXIT_FAILURE);
    }
  else
    ask_parameters();

}

bool
ReportMovement::
post_processing()
{
  if (base_type::post_processing() == true)
    return true;

  return false;
}


Succeeded 
ReportMovement::
process_data()
{

  const unsigned int min_frame_num =
    this->get_frame_num_to_process()==-1
    ? 1 : this->get_frame_num_to_process();
  const unsigned int max_frame_num =
    this->get_frame_num_to_process()==-1 
    ? this->get_time_frame_defs().get_num_frames() 
    : this->get_frame_num_to_process();

  for (unsigned int current_frame_num = min_frame_num;
       current_frame_num<=max_frame_num;
       ++current_frame_num)
    {
      const double start_time = 
	this->get_frame_start_time(current_frame_num);
      const double end_time = 
	this->get_frame_end_time(current_frame_num);

      cerr << "\nDoing frame " << current_frame_num
	   << ": from " << start_time << " to " << end_time << endl;

      //set_frame_num_to_process(current_frame_num);

      const RigidObject3DTransformation frame_transformation_to_reference =
	compose(this->get_rigid_object_transformation_to_reference(),
		this->get_motion().
		compute_average_motion_in_scanner_coords_rel_time(start_time, end_time));

      // now go through tracker data for this frame
      {
	const std::vector<double> sample_times =
	  this->get_motion().
	  get_rel_time_of_samples(start_time, end_time);

	if (sample_times.size() == 0)
	  error("No tracker samples between %g and %g (relative to scan start)",
		start_time, end_time);

	for (std::vector<double>::const_iterator iter=sample_times.begin();
	     iter != sample_times.end();
	     ++iter)
	  {
	    const RigidObject3DTransformation current_transformation_to_reference =
	      compose(this->get_rigid_object_transformation_to_reference(),
		      this->get_motion().
		      get_motion_in_scanner_coords_rel_time(*iter));

	    const RigidObject3DTransformation current_transformation_to_frame_ref =
	      compose(frame_transformation_to_reference,
		      current_transformation_to_reference.inverse());

	    const double RMSE_within_frame =
	      RigidObject3DTransformation::
	      RMSE(current_transformation_to_frame_ref,
		   this->reference_points.begin(), this->reference_points.end(),
		   this->reference_points.begin());

	    const double RMSE_from_ref =
	      RigidObject3DTransformation::
	      RMSE(current_transformation_to_reference,
		   this->reference_points.begin(), this->reference_points.end(),
		   this->reference_points.begin());

	    std::cout << *iter << " " << RMSE_within_frame << " " << RMSE_from_ref << '\n';
	  } // end of loop over tracker samples
      }
    } // end of loop over frames

  return Succeeded::yes;
}


END_NAMESPACE_STIR



USING_NAMESPACE_STIR

int main(int argc, char * argv[])
{
  bool move_to_reference=true;
  bool set_move_to_reference=false;
  bool set_frame_num_to_process=false;
  int frame_num_to_process=-1;
  while (argc>=2 && argv[1][1]=='-')
    {
      if (strcmp(argv[1], "--frame_num_to_process")==0)
	{
	  set_frame_num_to_process=true;
	  frame_num_to_process=atoi(argv[2]);
	  argc-=2; argv+=2;
	}
      else
	{
	  warning("Wrong option\n");
	  exit(EXIT_FAILURE);
	}
    }

  if (argc!=1 && argc!=2) {
    cerr << "Usage: " << argv[0] << " \\\n"
	 << "\t[--frame_num_to_process number]\\\n"
	 << "\t[par_file]\n";
    exit(EXIT_FAILURE);
  }
  ReportMovement application(argc==2 ? argv[1] : 0);
  if (set_move_to_reference)
    application.move_to_reference(move_to_reference);
  if (set_frame_num_to_process)
    application.set_frame_num_to_process(frame_num_to_process);
  Succeeded success =
    application.process_data();

  return success == Succeeded::yes ? EXIT_SUCCESS : EXIT_FAILURE;
}