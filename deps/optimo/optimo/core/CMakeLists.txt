# Tests:
IF (${BUILD_TESTING})
  GTEST(objects_tests)
  GTEST(numerical_gradient_tests)
  GTEST(numerical_hessian_tests)
ENDIF (${BUILD_TESTING})
