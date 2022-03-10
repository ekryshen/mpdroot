#ifndef MPD_BASE_HEADER_INCLUDED
#define MPD_BASE_HEADER_INCLUDED

namespace SI_units { // basic SI units taken from CLHEP project to avoid linking the library
static constexpr double millimeter = 1.;
static constexpr double centimeter = 10. * millimeter;
static constexpr double meter      = 1000. * millimeter;
static constexpr double micrometer = 1.e-6 * meter;
}; // namespace SI_units

#endif // #ifndef MPD_BASE_HEADER_INCLUDED
