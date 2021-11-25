#include "tau/tau.h"
#include "types.h"
#include "utils.h"


TEST(Cec2017State, Init) {

  cec_state_t state = cec_mk_state(CEC_2017, 10, NONE);

  REQUIRE_EQ(state.version_, CEC_2017);
  REQUIRE_EQ(state.suite_, NONE);
  REQUIRE_EQ(state.data_.rotates_->problem_nums_, 30);
  REQUIRE_EQ(state.data_.shifts_->problem_nums_, 30);
  REQUIRE_EQ(state.data_.shuffles_->problem_nums_, 30);
  REQUIRE(state.data_.biases_ == NULL);

}



