/*
    Copyright (C) 2011, Hammersmith Imanet Ltd
    Copyright (C) 2013, University College London
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    See STIR/LICENSE.txt for details
*/
/*!

  \file
  \ingroup recon_test
  
  \brief Test program for stir::PoissonLogLikelihoodWithLinearModelForMeanAndProjData

  \par Usage

  <pre>
  test_PoissonLogLikelihoodWithLinearModelForMeanAndProjData [proj_data_filename [ density_filename ] ]
  </pre>
  where the 2 arguments are optional. See the class documentation for more info.

  \author Kris Thielemans
*/

#include "stir/RunTests.h"
#include "stir/num_threads.h"
#include <iostream>

#include "stir/IO/OutputFileFormat.h"
#include "stir/recon_buildblock/distributable_main.h"
#include "stir/recon_buildblock/PresmoothingForwardProjectorByBin.h"
#include "stir/recon_buildblock/PostsmoothingBackProjectorByBin.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingSeparateProjectors.h"
#include "local/stir/motion/Transform3DObjectImageProcessor.h"
#include "local/stir/motion/NonRigidObjectTransformationUsingBSplines.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingProjMatrixByBin.h"
#include "stir/ProjData.h"
#include "stir/DiscretisedDensity.h"
#include "stir/HighResWallClockTimer.h"
#include "stir/ProjDataInMemory.h"

START_NAMESPACE_STIR


/*!
  \ingroup test
  \brief Test class for PoissonLogLikelihoodWithLinearModelForMeanAndProjData

*/
class ProjectorsWithMotionTests : public RunTests
{
public:
  //! Constructor that can take some input data to run the test with
  ProjectorsWithMotionTests();

  void run_tests();
protected:
  shared_ptr<ProjectorByBinPair> _projector_pair_sptr;
  std::string _disp_field_x;
  std::string _disp_field_y;
  std::string _disp_field_z;
  int _bspline_order;

  //! run the test
  /*! Note that this function is not specific to PoissonLogLikelihoodWithLinearModelForMeanAndProjData */
  void run_test_motion_fwrd_and_back();
};

ProjectorsWithMotionTests::
ProjectorsWithMotionTests()
{
    _projector_pair_sptr.reset(new ProjectorByBinPairUsingProjMatrixByBin);
    _projector_pair_sptr->parse("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/test_projectors.par");
    _bspline_order = 3;
    _disp_field_x = "/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/4_PET_regis/FMV_g4d1.nii";
    _disp_field_y = "/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/4_PET_regis/FMV_g4d2.nii";
    _disp_field_z = "/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/4_PET_regis/FMV_g4d3.nii";
}

void
ProjectorsWithMotionTests::
run_test_motion_fwrd_and_back() 
{
    // Create non rigid transformations
    shared_ptr<NonRigidObjectTransformationUsingBSplines<3,float> > fwrd_non_rigid;

    fwrd_non_rigid.reset(
                new NonRigidObjectTransformationUsingBSplines<3,float>(
                    _disp_field_x,
                    _disp_field_y,
                    _disp_field_z,
                    _bspline_order));
    // Create image processors
    shared_ptr<Transform3DObjectImageProcessor<float> > fwrd_transform, back_transform;
    fwrd_transform.reset( new Transform3DObjectImageProcessor<float>(fwrd_non_rigid) );
    back_transform.reset( new Transform3DObjectImageProcessor<float>(*fwrd_transform) );
    back_transform->set_do_transpose(!fwrd_transform->get_do_transpose());
/*
    shared_ptr<DiscretisedDensity<3,float> > input;
    input.reset( DiscretisedDensity<3,float>::read_from_file("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/3_NAC/NAC_g24.nii"));
    shared_ptr<DiscretisedDensity<3,float> > output_fwrd, output_back;
    output_fwrd.reset(input->clone());
    output_back.reset(input->clone());

    std::cout << "\nDoing the forward projection...\n";
    output_fwrd->fill(0.);
    fwrd_transform->apply(*output_fwrd,*input);

    std::cout << "OK!\nDoing the backward projection...\n";
    output_back->fill(0.);
    back_transform->apply(*output_back,*output_fwrd);
    std::cout << "OK!\n";

    shared_ptr<OutputFileFormat<DiscretisedDensity<3,float> > > output_file_format_sptr;
    output_file_format_sptr = OutputFileFormat<DiscretisedDensity<3,float> >::default_sptr();
    output_file_format_sptr->write_to_file("/Users/rich/Desktop/output_fwrd", *output_fwrd);
    output_file_format_sptr->write_to_file("/Users/rich/Desktop/output_back", *output_back);*/
/*
    // Create projectors
    shared_ptr<PresmoothingForwardProjectorByBin> fwrd_projector;
    shared_ptr<PostsmoothingBackProjectorByBin>   back_projector;
    fwrd_projector.reset(
                new PresmoothingForwardProjectorByBin(
                    this->_projector_pair_sptr->get_forward_projector_sptr(),
                    fwrd_transform));
    back_projector.reset(
                new PostsmoothingBackProjectorByBin(
                    this->_projector_pair_sptr->get_back_projector_sptr(),
                    back_transform));
    // Finally create the projector pair and store it in the vector
    _projector_pair_sptr.reset(
                new ProjectorByBinPairUsingSeparateProjectors(
                    fwrd_projector,
                    back_projector));
*/
    shared_ptr<ProjData> projdata_in = ProjData::read_from_file("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/1_sinograms/sinogram_g24.hs");
    shared_ptr<DiscretisedDensity<3,float> > im;
    im.reset( DiscretisedDensity<3,float>::read_from_file("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/3_NAC/NAC_g24.nii"));

    ProjDataInMemory projdata(projdata_in->get_exam_info_sptr(), projdata_in->get_proj_data_info_sptr());

    _projector_pair_sptr->set_up(projdata.get_proj_data_info_sptr(), im);
    _projector_pair_sptr->get_forward_projector_sptr()->set_input(im);
    std::cout << "\nthe max of the input image is: " << im->sum() << "\n";
    _projector_pair_sptr->get_back_projector_sptr()->start_accumulating_in_new_image();
    std::cout << "\nthe max of the input image is (should be unchanged): " << im->sum() << "\n";

    // Start the timer
    HighResWallClockTimer t;
    t.reset();
    t.start();

    std::cout << "\ndoing the forward projection...\n";
    _projector_pair_sptr->get_forward_projector_sptr()->forward_project(projdata);
    std::cout << "OK!\n";

    std::cout << "\ndoing the backwards projection...\n";
    _projector_pair_sptr->get_back_projector_sptr()->back_project(projdata);
    _projector_pair_sptr->get_back_projector_sptr()->get_output(*im);
    std::cout << "OK!\n";

    t.stop();
    std::cout << "Total Wall clock time: " << t.value() << " seconds" << std::endl;

    shared_ptr<OutputFileFormat<DiscretisedDensity<3,float> > > output_file_format_sptr;
    output_file_format_sptr = OutputFileFormat<DiscretisedDensity<3,float> >::default_sptr();
    output_file_format_sptr->write_to_file("/Users/rich/Desktop/test", *im);
}

void
ProjectorsWithMotionTests::
run_tests()
{
  std::cerr << "Tests for ProjectorsWithMotionTests\n";

  this->run_test_motion_fwrd_and_back();
} 

END_NAMESPACE_STIR


USING_NAMESPACE_STIR

#ifdef STIR_MPI
int stir::distributable_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
  set_default_num_threads();

  ProjectorsWithMotionTests tests;
  tests.run_tests();

  return tests.main_return_value();
}
