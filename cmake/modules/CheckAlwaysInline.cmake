# check support for __attribute__((always_inline))
# this test is deprecated and will be removed after DuMuX 2.4
CHECK_CXX_SOURCE_COMPILES("
   void __attribute__((always_inline)) foo(void) {}
   int main(void)
   { foo(); return 0; }"
   HAVE_ATTRIBUTE_ALWAYS_INLINE)
