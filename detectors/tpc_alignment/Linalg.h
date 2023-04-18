#ifndef LINALG_HH
#define LINALG_HH
#include <vector>

namespace TpcAlignment{

/// <summary>
/// Algebraic programs
/// </summary>
class Linalg
{
public:
Linalg() = delete;

/// <summary>
/// Solving equation Ax=B
/// Matrices are arranged in columns
/// </summary>
/// <param name="vA">left side </param>
/// <param name="B">right side</param>
/// <param name="X">calculated X</param>
static void lsolve(std::vector<double> & vA, 
					std::vector<double> & B, 
					std::vector<double> & X);

};

}

#endif