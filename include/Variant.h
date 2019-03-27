#pragma once
#include<sstream>
#include<vector>
#include<stdexcept>
namespace XFITTER_PARS{
/*
 Variant is a class used for transferring parameters to reations
 It can encapsulate one of several types:
 *const double*
 *std::string
 *int
 *vector<const double*>

 We use pointers for double-typed parameters because it is assumed
 that any of them can change during the minimization

 Variant is used mainly as return type of TermData::getParam
*/
class Variant{
public:
  enum Type{None=0,DoublePtr=1,String=2,Int=3,Array=4};
  Variant();
  Variant(const Variant&);
  Variant(const double*);
  Variant(const std::string&);
  Variant(const char*);
  Variant(int);
  Variant(const std::vector<const double*>&);
  ~Variant();
  Variant&operator=(const Variant&);
  operator const double*()const;
  operator std::string()const;
  operator int()const;
  operator const std::vector<const double*>&()const;
  Type type()const;
  bool isNone()const;
  bool isDoublePtr()const;
  bool isString()const;
  bool isInt()const;
  bool isArray()const;
  friend std::ostream&operator<<(std::ostream&,const Variant&);
  class bad_cast;
private:
  Type _type;
  union{
    const double*_ptr;
    std::string _string;
    int _int;
    std::vector<const double*>_array;
  };
};
std::ostream&operator<<(std::ostream&,const Variant&);
std::ostream&operator<<(std::ostream&,Variant::Type);
const char*to_string(Variant::Type);
class Variant::bad_cast:public std::runtime_error{
public:
  bad_cast(const Variant&,Variant::Type);
};
}
