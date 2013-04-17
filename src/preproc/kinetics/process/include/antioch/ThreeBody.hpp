#ifndef _THREE_BODY_REACTION_
#define _THREE_BODY_REACTION_

#include "antioch/Process.hpp"
#include <map>

namespace Antioch{
/*!\file ThreeBody.hpp
 * \brief contains the corresponding class
 *
 * \class ThreeBody
 * \ingroup kinproc
 * \brief Three-body reaction
 *
 * Three body reactions are typically decompositions.
 *
 *  For a three-body reaction:
 *  \f[
 *      \ce{A + B + M -> AB + M}
 *  \f]
 *  we have the rate:
 *  \f[
 *     \frac{\mathrm{d}[\ce{A}]}{\mathrm{d}t} = k_{\ce{AB}} [\ce{M}] [\ce{A}] [\ce{B}]
 *  \f]
 *  with \f$[\ce{M}]\f$ defined as:
 *  \f[
 *      [\ce{M}] = \sum_i \alpha_i [\ce{S}_i]
 *  \f]
 *  with \f$[\ce{S}_i]\f$ the concentration of species \f$\ce{S}_i\f$
 *  By default, \f$\alpha_i = 1\f$.
 */

class ThreeBody:public Process
{
   public:
        ThreeBody(){}
        ~ThreeBody(){}

        void addOneCoefficient(const std::string &molecule, double coef) {alpha[molecule] = coef;}
        double getCoefficient(const std::string &molecule)         const;

   private:
       std::map<std::string,double> alpha;
};
}
#endif
