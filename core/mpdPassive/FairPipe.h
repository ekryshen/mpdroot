/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             * 
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *  
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef PIPE_H
#define PIPE_H

#include "FairModule.h"                 // for FairModule

#include "Rtypes.h"                     // for FairPipe::Class, ClassDef, etc
#include <string>                       // for string

class FairPipe : public FairModule
{
  public:
    FairPipe(const char* name, const char* Title="PND Pipe");
    FairPipe();

    virtual ~FairPipe();
    void ConstructGeometry();
    void ConstructASCIIGeometry();

    ClassDef(FairPipe,1) //PNDPIPE

};

#endif //PIPE_H

