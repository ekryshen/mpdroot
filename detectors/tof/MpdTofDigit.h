//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPDTOFDIGIT_H
#define	__MPDTOFDIGIT_H 1

#include "TObject.h"
#include <iostream>
//------------------------------------------------------------------------------------------------------------------------
class MpdTofDigit : public TObject
{
    Float_t fAmplitude = -1.f;
    Float_t fTime = -1.f;
    Short_t fSector = -1;
    Short_t fPlane = -1;
    Short_t fStrip = -1;
    Short_t fSide = -1;

public:
    MpdTofDigit() = default;
    MpdTofDigit(Short_t sector, Short_t plane, Short_t strip, Short_t side, Float_t t, Float_t a);
    MpdTofDigit(const MpdTofDigit*, Float_t t, Float_t a);
    virtual ~MpdTofDigit() = default;

    Short_t GetSector()    const { return fSector;    }
    Short_t GetPlane()     const { return fPlane;     }
    Short_t GetStrip()     const { return fStrip;     }
    Short_t GetSide()      const { return fSide;      }
    Float_t GetAmplitude() const { return fAmplitude; }
    Float_t GetTime()      const { return fTime;      }

    void SetSector   (Short_t v) { fSector = v; }
    void SetPlane    (Short_t v) { fPlane = v; }
    void SetStrip    (Short_t v) { fStrip = v; }
    void SetSide     (Short_t v) { fSide = v; }
    void SetAmplitude(Float_t v) { fAmplitude = v; }
    void SetTime     (Float_t v) { fTime = v; }

    void print(const char* comment = nullptr, std::ostream& os = std::cout)const;

    bool operator< ( const MpdTofDigit& rhs)const // sort by strip guid(fSector + fPlane + fStrip) 
    {
	if(this->fSector != rhs.fSector) return (this->fSector < rhs.fSector);	

	if(this->fPlane != rhs.fPlane) return (this->fPlane < rhs.fPlane);	
	
	return (this->fStrip < rhs.fStrip);					
    }

    ClassDef(MpdTofDigit, 2);
};
//------------------------------------------------------------------------------------------------------------------------
#endif
