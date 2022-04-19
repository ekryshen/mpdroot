#ifndef MPDGENERATORFACTORY_H
#define MPDGENERATORFACTORY_H

#include <map>
#include <memory>
#include <functional>

#include <TString.h>
#include <TSystem.h>

#include "MpdGeneratorType.h"

#include "FairGenerator.h"
#include "FairRunSim.h"

class MpdGenericGenCreator;
class MpdUrqmdGenCreator;
class MpdVHLLEGenCreator;
class MpdFLUIDGenCreator;
// class MpdPARTGenCreator;
// class MpdIONGenCreator;
// class MpdBOXGenCreator;
class MpdHSDGenCreator;
class MpdLAQGSMGenCreator;
class MpdFactoryMadeGenerator;

class MpdGeneratorsFactory {
public:
   MpdGeneratorsFactory();
   virtual ~MpdGeneratorsFactory();

   std::shared_ptr<MpdFactoryMadeGenerator> create(MpdGenerator &Gen);

protected:
   typedef std::map<MpdGeneratorType, MpdGenericGenCreator *> CreatorsMap;
   CreatorsMap                                                creators;

   template <class C>
   void add(const MpdGeneratorType &id);
};

class MpdFactoryMadeGenerator : public std::enable_shared_from_this<MpdFactoryMadeGenerator> {
public:
   MpdFactoryMadeGenerator(std::shared_ptr<FairGenerator> fg, std::function<Int_t(FairRunSim *)> pa, TString dataFile)
      : f_gen(fg), f_pa(pa), f_dataFile(dataFile)
   {
   }
   std::shared_ptr<MpdFactoryMadeGenerator> getptr() { return shared_from_this(); }

   virtual const std::shared_ptr<FairGenerator> &GetGeneratorPtr() const { return f_gen; }
   virtual Int_t                                 PostActions(FairRunSim *fRun) { return f_pa(fRun); }
   // virtual void SetPostActions(std::function<Int_t(FairRunSim*)> pa){f_pa = pa;}

protected:
   std::shared_ptr<FairGenerator>     f_gen;
   TString                            f_dataFile = "";
   std::function<Int_t(FairRunSim *)> f_pa;
};

class MpdGenericGenCreator {
public:
   virtual FairGenerator                     *create(TString &inFile, Int_t &nStartEvent, Int_t &nEvents) = 0;
   virtual std::function<Int_t(FairRunSim *)> postActions(TString &)
   {
      return [](FairRunSim *) { return 0; };
   };
};

class MpdUrqmdGenCreator : public MpdGenericGenCreator {
public:
   virtual FairGenerator *create(TString &inFile, Int_t &nStartEvent, Int_t &nEvents) override;
};

class MpdVHLLEGenCreator : public MpdGenericGenCreator {
public:
   virtual FairGenerator *create(TString &inFile, Int_t &, Int_t &) override;
};

class MpdFLUIDGenCreator : public MpdGenericGenCreator {
public:
   virtual FairGenerator *create(TString &inFile, Int_t &nStartEvent, Int_t &) override;
};

/*class MpdPARTGenCreator : public MpdGenericGenCreator
{
public:
   virtual FairGenerator * create (TString &, Int_t &, Int_t &) override;
};*/

/*class MpdIONGenCreator : public MpdGenericGenCreator
{
public:
   virtual FairGenerator * create (TString &, Int_t &, Int_t &) override;
};*/

/*class MpdBOXGenCreator : public MpdGenericGenCreator
{
public:
   virtual FairGenerator * create (TString &, Int_t &, Int_t &) override;
};*/

class MpdHSDGenCreator : public MpdGenericGenCreator {
public:
   virtual FairGenerator *create(TString &inFile, Int_t &nStartEvent, Int_t &nEvents) override;
};

class MpdLAQGSMGenCreator : public MpdGenericGenCreator {
public:
   virtual FairGenerator                         *create(TString &inFile, Int_t &nStartEvent, Int_t &nEvents) override;
   virtual std::function<Int_t(FairRunSim *fRun)> postActions(TString &inFile) override;
};

#endif
