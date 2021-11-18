 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# Create a library called "libCluster" which includes the source files given in
# the array

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/clustering
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
) 

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
ClusterElement.cxx
Cluster.cxx
ClusterCollection.cxx
ClManager.cxx
)

Set(HEADERS)
Set(LINKDEF ClusterLinkDef.h)
Set(LIBRARY_NAME Cluster)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
