//
// $Id$
//
/*!
  \file 
 
  \brief defines the warning() function

  \author Kris Thielemans
  \author PARAPET project

  $Date$

  $Revision$

*/
/*
    Copyright (C) 2000 PARAPET partners
    Copyright (C) 2000- $Date$, IRSL
    See STIR/LICENSE.txt for details
*/
#include "stir/common.h"

#include <cstdarg>
#include <stdlib.h>

START_NAMESPACE_STIR

void warning(const char *const s, ...)
{
  va_list ap;
  va_start(ap, s);
  vfprintf(stderr, s, ap);
}

END_NAMESPACE_STIR
