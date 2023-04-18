# Create a library called "libTpcAlignment" .

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
	CorrectionMatrix.cxx
	Debug.cxx
	EigenCalculator.cxx
	EigenCoefficients.cxx
	FileHelper.cxx
	Linalg.cxx
	Point.cxx
	Rays.cxx
	Runner.cxx
	SystemOfLinearEquations.cxx
	TpcAlignmentTask.cxx
	TpcCorrectionTask.cxx
	Track.cxx
)

Set(HEADERS)

Set(LINKDEF tpcAlignmentLinkDef.h)
Set(LIBRARY_NAME TpcAlignment)
Set(DEPENDENCIES Core Base )

GENERATE_LIBRARY()


