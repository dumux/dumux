# check support of constexpr
CHECK_CXX_SOURCE_COMPILES("
  int main(void) {
    constexpr double g = 9.81;
    return 0; }"
  HAVE_CONSTEXPR)
