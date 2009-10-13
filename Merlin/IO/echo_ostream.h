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
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef _h_echo_ostream 
#define _h_echo_ostream 1

#include <iostream>

// template class basic_echo_ostream
//
// Use to create an ostream which re-directs
// output to two additional ostream objects.
//
// Note that these implementation is probably
// very inefficient, and could be made better
// by somebody who understands the std lib better.
//

template<class T, class Tr=char_traits<T> >
class basic_echo_ostream : public std::basic_ostream<T,Tr> {
public:

class echo_streambuf : public std::basic_streambuf<T,Tr> {
    public:

        echo_streambuf(std::basic_ostream<T,Tr>& os1, std::basic_ostream<T,Tr>& os2)
                : _os1(os1),_os2(os2) {}

    protected:

        virtual int_type overflow(int_type c) {
            _os1<<T(c);
            _os2<<T(c);
            return char_traits<T>::not_eof(c);
        }

    private:

        std::basic_ostream<T,Tr>& _os1;
        std::basic_ostream<T,Tr>& _os2;
    };

    // special buffer to redirect output

    basic_echo_ostream(std::basic_ostream<T,Tr>& first, std::basic_ostream<T,Tr>& second)
            : std::basic_ostream<T,Tr>(&buff),buff(first,second)
    {}

private:

    echo_streambuf buff;
};

typedef basic_echo_ostream<char> echo_ostream;

#endif
