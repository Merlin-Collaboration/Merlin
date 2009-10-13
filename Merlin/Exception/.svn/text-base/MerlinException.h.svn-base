/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:54 $
// $Revision: 1.3 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef MerlinException_h
#define MerlinException_h 1

#include <string>

//	Root class for all Merlin exceptions.

class MerlinException 
{
  public:

      explicit MerlinException (const std::string& s);
      MerlinException ();
      const char* Msg () const;

  protected:

      void SetMsg (const std::string& s);
      void AppendMsg (const std::string& s);

  private:
      std::string msg;
};

inline MerlinException::MerlinException (const std::string& s)
  : msg(s)
{}

inline MerlinException::MerlinException ()
{}

inline const char* MerlinException::Msg () const
{
	return msg.c_str();
}

inline void MerlinException::SetMsg (const std::string& s)
{
	msg=s;
}

inline void MerlinException::AppendMsg (const std::string& s)
{
	msg+=s;
}

#endif
