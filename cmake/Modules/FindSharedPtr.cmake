# check if make_shared works
macro(CHECK_MAKE_SHARED)
  include(CheckIncludeFileCXX)
  set(CMAKE_REQUIRED_FLAGS "--std=c++0x -DHAVE_MEMORY=${HAVE_MEMORY} -DHAVE_TR1_MEMORY=${HAVE_TR1_MEMORY}")
  
  if(SHARED_PTR_NAMESPACE EQUAL "boost")
    CHECK_INCLUDE_FILE_CXX("boost/make_shared.hpp" HAVE_BOOST_MAKE_SHARED_HPP)
  endif(SHARED_PTR_NAMESPACE EQUAL "boost")
  
  CHECK_CXX_SOURCE_COMPILES("
    #if defined(HAVE_MEMORY)
    # include <memory>
    #endif
    #if defined(HAVE_TR1_MEMORY)
    # include <tr1/memory>
    #endif
    #if defined(HAVE_BOOST_SHARED_PTR_HPP) && defined(HAVE_BOOST_MAKE_SHARED_HPP)
    # include <boost/shared_ptr.hpp>
    # include <boost/make_shared.hpp>
    #endif
    
    int main(void)
    {
       ${SHARED_PTR_NAMESPACE}::make_shared<int>(3);
       return 0;
    }
    " HAVE_MAKE_SHARED)
endmacro(CHECK_MAKE_SHARED)

# check location of shared_ptr header file and the necessary namespace
set(CMAKE_REQUIRED_FLAGS "--std=c++0x")
include(CheckCXXSourceCompiles)

# search namespace
foreach(SHARED_PTR_NAMESPACE_ "std" "tr1" "std::tr1" "boost")
  CHECK_CXX_SOURCE_COMPILES("
    #include <memory>
    #include <string>

    using ${SHARED_PTR_NAMESPACE_}::shared_ptr;
    using namespace std;

    int main(void)
    {
      shared_ptr<string> test_ptr(new string(\"test string\"));
      return 0;
    }"
    SHARED_PTR_NAMESPACE_FOUND)

  if(SHARED_PTR_NAMESPACE_FOUND)
    #search header name
    foreach(SHARED_PTR_HEADER_ "<memory>" "<tr1/memory>" "<boost/shared_ptr.hpp>")
      CHECK_CXX_SOURCE_COMPILES("
        # include ${SHARED_PTR_HEADER_}
        #include <string>

        using ${SHARED_PTR_NAMESPACE_}::shared_ptr;
        using namespace std;

        int main(void)
        {
          shared_ptr<string> test_ptr(new string(\"test string\"));
          return 0;
        }"
        SHARED_PTR_HEADER_FOUND)
    
      if(SHARED_PTR_HEADER_FOUND)
        # save result
        set(SHARED_PTR_NAMESPACE ${SHARED_PTR_NAMESPACE_})
        set(SHARED_PTR_HEADER ${SHARED_PTR_HEADER_})
        CHECK_MAKE_SHARED()
        return()
      endif(SHARED_PTR_HEADER_FOUND)
    endif(SHARED_PTR_NAMESPACE_FOUND)
  endforeach(SHARED_PTR_HEADER_ "memory" "tr1/memory" "boost/shared_ptr.hpp")
endforeach(SHARED_PTR_NAMESPACE_ "std" "tr1" "std::tr1" "boost")
