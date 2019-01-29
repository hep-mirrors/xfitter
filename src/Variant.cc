#ifndef XFITTER_VARIANT
#define XFITTER_VARIANT
#include"Variant.h"
namespace XFITTER_PARS{
Variant::Variant():                                  _type{None}{}
Variant::Variant(const double*p):                    _type{DoublePtr},_ptr{p}{}
Variant::Variant(const std::string&s):               _type{String}   ,_string{s}{}
Variant::Variant(const char*s):                      _type{String}   ,_string{s}{}
Variant::Variant(int i):                             _type{Int}      ,_int{i}{}
Variant::Variant(const std::vector<const double*>&v):_type{Array}    ,_array{v}{}
Variant::Variant(const Variant&o):_type{o._type}{
	switch(_type){
		case None:break;
		case DoublePtr:_ptr=o._ptr;break;
		//One can't simply assign to initialize a string in union
		//because assignment requires left value to be valid
		//That's why we have this weird-looking construct-in-place
		case String:new(&_string)std::string(o._string);break;
		case Int:_int=o._int;break;
		case Array:new(&_array)std::vector<const double*>(o._array);break;
	}
}
Variant::~Variant(){
switch(_type){
	case String:_string.std::string::~string();break;
	case Array:_array.std::vector<const double*>::~vector<const double*>();break;
	default:break;
}
}
Variant&Variant::operator=(const Variant&o){
  
	if(_type==String){
		if(o._type==String){
			_string=o._string;
			return *this;
    //When type changes from string to something, we need to destroy the old string
		}else _string.std::string::~string();
	}else if(_type==Array){
		if(o._type==Array){
			_array=o._array;
			return *this;
    //When type changes from array to something, we need to destroy the old array
		}else _array.std::vector<const double*>::~vector<const double*>();
	}
	_type=o._type;
	switch(_type){
		case DoublePtr:_ptr=o._ptr;break;
		case String:new(&_string)std::string(o._string);break;
		case Int:_int=o._int;break;
		case Array:new(&_array)std::vector<const double*>(o._array);break;
		default:break;
	}
	return *this;
}
std::string bad_cast_message(const Variant&v,Variant::Type castTo){
	std::ostringstream ss;
	ss<<"Failed to cast Variant "<<v<<" from "<<to_string(v.type())<<" to "<<to_string(castTo);
	return ss.str();
}
Variant::bad_cast::bad_cast(const Variant&v,Variant::Type to):std::runtime_error(bad_cast_message(v,to)){}
Variant::operator const double*()const{
	if(_type!=DoublePtr)throw bad_cast(*this,DoublePtr);
	return _ptr;
}
Variant::operator std::string()const{
	if(_type==String)return _string;
	if(_type==Int)return std::to_string(_int);
	throw bad_cast(*this,String);
};
Variant::operator int()const{
	switch(_type){
		case Int:return _int;
		case String:
			try{return std::stoi(_string);
			}catch(std::invalid_argument&ex){
				break;
			}catch(std::out_of_range&ex){
				break;
			}
		default:break;
	}
	throw bad_cast(*this,Int);
}
Variant::operator const std::vector<const double*>&()const{
	if(_type!=Array)throw bad_cast(*this,Array);
	return _array;
}
const char*to_string(Variant::Type t){
	static const char*a[]={"None","DoublePtr","String","Int","Array"};
	return a[int(t)];
}
std::ostream&operator<<(std::ostream&os,const Variant&v){
  if(v._type==Variant::String)return os<<v._string;
  if(v._type==Variant::Int)return os<<v._int;
  if(v._type==Variant::DoublePtr)return os<<v._ptr;
  return os<<to_string(v._type);
}
std::ostream&operator<<(std::ostream&os,Variant::Type t){
  return os<<to_string(t);
}
Variant::Type Variant::type()const{return _type;}
bool Variant::isNone()     const{return _type==None;}
bool Variant::isDoublePtr()const{return _type==DoublePtr;}
bool Variant::isString()   const{return _type==String;}
bool Variant::isInt()      const{return _type==Int;}
bool Variant::isArray()    const{return _type==Array;}
}
#endif
