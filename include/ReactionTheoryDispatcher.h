#pragma once

#include <string>
#include <list>

class ReactionTheory;

using namespace std::string;
using namespace std::list;

/**
  @class ReactionTheoryDispatcher

  @brief Class manages reaction theory objects

  It creates new objects for reaction theory corresponding to it's type. This is a singleton class.

  @author A.Sapronov <sapronov@ifh.de>

  @version 0.1
  @date 2016/01/21
  */

class ReactionTheoryDispatcher 
{
 public:
  static ReactionTheoryDispatcher &getInstance(){
    static ReactionTheoryDispatcher rtdInstance;
    return rtdInstance;
  }

 private:
  ReactionTheoryDispatcher() {};
  ~ReactionTheoryDispatcher(){};

  ReactionTheoryDispatcher(const ReactionTheoryDispatcher &);
  void operator =(const ReactionTheoryDispatcher &);

 public:
  ReactionTheory &getReactionTheory(const string &reaction_type);
  
 private:
  /// list of created theories
  list<ReactionTheory*> _rt_list;
};
