//
// $Id$
//
/*!

  \file
  \ingroup buildblock
  \brief Declaration of class 

  \author Kris Thielemans
  \author Sanida Mustafovic

  $Date$
  $Revision$
*/
/*
    Copyright (C) 2000- $Date$, IRSL
    See STIR/LICENSE.txt for details
*/

#ifndef __stir_RegisteredParsingObject_H__
#define __stir_RegisteredParsingObject_H__


#include "stir/ParsingObject.h"
#include <string>

#ifndef STIR_NO_NAMESPACE
using std::string;
#endif


START_NAMESPACE_STIR


template <typename Base>
class AddParser : public Base, public ParsingObject
{};



/*!
  \brief Parent class for all leaves in a RegisteredObject hierarchy that
  do parsing of parameter files.

  \see RegisteredObject
  
  RegisteredParsingObject::read_from_stream is implemented in terms of
  ParsingObject::parse.

  Requirements on the class Base:
  - It needs to be derived from RegisteredObject<Base>

  Requirements on the class Derived:
  - It needs to have a static member static const char * const registered_name
  - It needs to have a default constructor
  - It needs to be derived from RegisteredParsingObject<Derived,Base,Parent>

  Requirements on the class Parent:
  - It needs to be derived from ParsingObject
  - It needs to be derived from Base
*/
template <typename Derived, typename Base, typename Parent = AddParser<Base> >
class RegisteredParsingObject : public Parent
{
public:
  //! Construct a new object (of type Derived) by parsing the istream
  /*! When the istream * is 0, questions are asked interactively. 
  
      Currently, the return value is a Base*. Preferably, it should be a 
      Derived*, but it seems the registration machinery would need extra 
      (unsafe) reinterpret_casts to get that to work.
      (TODO find a remedy).
  */
  inline static Base* read_from_stream(istream*); 

  //! Returns  Derived::registered_name
  inline string get_registered_name() const;
  //! Returns a string with all parameters and their values, in a form suitable for parsing again
  inline string parameter_info();

protected:
#ifdef _MSC_VER
public:
#endif
  //! A helper class to allow automatic registration.
  struct RegisterIt
  {
    RegisterIt()
    {
      //std::cerr << "Adding " << Derived::registered_name <<" to registry"<<std::endl;
      registry().add_to_registry(Derived::registered_name, read_from_stream);  
    }
    ~RegisterIt()
    {
#if 0
      // does not work yet, as registry might be destructed before this
      // RegisterIt object. A solution to this problem is coming up.
      cerr << "In RegisterIt destructor for " << Derived::registered_name<<endl;
      cerr <<"Current keys: ";
      registry().list_keys(cerr);
      registry().remove_from_registry(Derived::registered_name);
#endif
    }
  };
  // RegisterIt needs to be a friend to have access to registry()
  friend struct RegisterIt;
  
};


END_NAMESPACE_STIR

#include "stir/RegisteredParsingObject.inl"

#endif

