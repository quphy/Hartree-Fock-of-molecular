#ifndef mixing_h
#define mixing_h

#include "scf.h"

namespace hf::scf{
template<class State>
struct simple_mixing{
double alpha;
State state;
 //当调用比如simping_mixing updater(previous_state)的时候,捕获 old_state
 std::pair<scf::Update<State>, simple_mixing> operator()(const State & state) {
    simple_mixing<State> new_struct{alpha, state};
    //将previous state存成old_state，方便后续的闭包调用
    const auto & old_state = state;
    //renewed.first(new_state)的时候的时候，混合
    const scf::Update<State> updater = [old_state, alpha = this->alpha]
        (const State & new_state) -> State {
      return alpha* new_state+(1.0-alpha) * old_state;
    };
    return  {updater, new_struct};
  }

};
}






#endif