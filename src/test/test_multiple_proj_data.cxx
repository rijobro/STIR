//
//
/*!

  \file
  \ingroup test

  \brief Test program for stir::ProjDataInfo hierarchy

  \author Sanida Mustafovic
  \author Kris Thielemans
  \author PARAPET project

*/
/*
    Copyright (C) 2000 PARAPET partners
    Copyright (C) 2000- 2011, Hammersmith Imanet Ltd
    Copyright (C) 2018, University College London
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

#include "stir/num_threads.h"
#include "stir/RunTests.h"
#include "stir/MultipleProjData.h"
#include "stir/DynamicProjData.h"
#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearKineticModelAndDynamicProjectionData.h"
#include "stir/recon_buildblock/PoissonLogLikelihoodWithLinearModelForMeanAndGatedProjDataWithMotion.h"

#ifndef STIR_NO_NAMESPACES
using std::cerr;
using std::setw;
using std::endl;
using std::min;
using std::max;
using std::size_t;
#endif

START_NAMESPACE_STIR

/*!
  \ingroup test
  \brief Test class for MultipleProjData
*/

class MultipleProjDataTests: public RunTests
{
public:  
  void run_tests();
};

void
MultipleProjDataTests::run_tests()

{ 
  std::cout << "-------- Testing MultipleProjData --------\n";
  {/*
    // Test on the empty constructor
    std::cout << "\n\nTesting empty constructor.\n";
    MultipleProjData ob1;
    std::cout << "OK!\n";
*/
    // Test with parser
    std::cout << "\n\nTesting MultipleProjData with parser.\n";
    shared_ptr<MultipleProjData> obj2;
    obj2 = MultipleProjData::read_from_file("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/1_sinograms/sinogram.txt");
    std::cout << "OK!\n";

    for (int i=0;i<10;i++) {
        std::cout << "\n\nTesting DynamicProjData with parser (number " << i << ").\n";
        shared_ptr<DynamicProjData> obj3;
        obj3 = DynamicProjData::read_from_file("/Users/rich/Documents/OneDrive-UCL/Data/MCIR/Output-MCIR/1_sinograms/sinogram.txt");
        std::cout << "num frames: " << std::flush << obj3->get_num_gates() << "\n";
        std::cout << "start time: " << std::flush << obj3->get_start_time_in_secs_since_1970() << "\n";
        std::cout << "OK!\n";
    }
  }
}

END_NAMESPACE_STIR


USING_NAMESPACE_STIR

int main()
{
  set_default_num_threads();

  {
    MultipleProjDataTests tests;
    tests.run_tests();
    if (!tests.is_everything_ok())
      return tests.main_return_value();
  }
}
