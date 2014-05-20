//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
#ifndef ANTIOCH_CHEMKIN_PARSER_H
#define ANTIOCH_CHEMKIN_PARSER_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/string_utils.h"
#include "antioch/parsing_enum.h"

//ChemKin

//C++
#include <sstream>
#include <string>
#include <vector>
#include <map>

namespace Antioch{

  template <typename NumericType = double>
  class ChemKinParser{
        public:
          ChemKinParser(const std::string &filename);
          ~ChemKinParser();

//// first local pointers
         /*! Read header of file, go to interesting part*/
         bool initialize();

         /*! go to next reaction*/
         bool reaction();

         /*! go to next rate constant*/
         bool rate_constant(const std::string & kinetics_model);

         /*! return true if there's a Troe block*/
         bool Troe();

         /*! return reaction id, 0 if not provided*/
         const std::string reaction_id() const;

         /*! return reaction equation */
         const std::string reaction_equation() const;

         /*! return reaction chemical process*/
         const std::string reaction_chemical_process() const;

         /*! return reaction kinetics model*/
         const std::string reaction_kinetics_model(const std::vector<std::string> &kinetics_models) const;

         /*! return reversible state*/
         bool reaction_reversible() const;

         /*! return pairs of reactants and stoichiometric coefficients*/
         bool reactants_pairs(std::vector<std::pair<std::string,int> > & reactants_pair) const;

         /*! return pairs of products and stoichiometric coefficients*/
         bool products_pairs(std::vector<std::pair<std::string,int> > & products_pair) const;

         /*! return true if "name" attribute is found with value "k0"*/
         bool is_k0(unsigned int nrc, const std::string & kin_model) const;

         /*! return index of k0 (0 or 1)*/
         unsigned int where_is_k0(const std::string & kin_model) const;

         /*! return true if pre exponentiel coefficient*/
         bool rate_constant_preexponential_parameter(    NumericType & A,    std::string & A_unit,    std::string & def_unit) const;

         /*! return true if beta coefficient*/
         bool rate_constant_power_parameter(             NumericType & b,    std::string & b_unit,    std::string & def_unit) const;

         /*! return true if activation energie*/
         bool rate_constant_activation_energy_parameter(NumericType & Ea,   std::string & Ea_unit,   std::string & def_unit) const;

         /*! return true if D coefficient*/
         bool rate_constant_Berthelot_coefficient_parameter(NumericType & D,    std::string & D_unit,    std::string & def_unit) const;

         /*! return true if Tref*/
         bool rate_constant_Tref_parameter(              NumericType & Tref, std::string & Tref_unit, std::string & def_unit) const;

         /*! return true if lambda*/
         bool rate_constant_lambda_parameter(       std::vector<NumericType> & lambda, std::string & lambda_unit, std::string & def_unit) const;

         /*! return true if sigma*/
         bool rate_constant_cross_section_parameter(std::vector<NumericType> & sigma,  std::string & sigma_unit,  std::string & def_unit) const;

         /*! return true if a Kooij is called Arrhenuis*/
         bool verify_Kooij_in_place_of_Arrhenius() const;

         /*! return true if efficiencies are found*/
         bool efficiencies(std::vector<std::pair<std::string,double> > & par_values) const;

         /*! return true is alpha*/
         bool Troe_alpha_parameter(NumericType & alpha, std::string & alpha_unit, std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T1_parameter(   NumericType & T1,    std::string & T1_unit,    std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T2_parameter(   NumericType & T2,    std::string & T2_unit,    std::string & def_unit) const;

         /*! return true is alpha*/
         bool Troe_T3_parameter(   NumericType & T3,    std::string & T3_unit,    std::string & def_unit) const;

        private:

          /*! Never use default constructor*/
          ChemKinParser();
          std::ofstream                    _doc;
          bool                             _reversible;
          unsigned int                     _nrates;
          NumericType                      _A;
          NumericType                      _b;
          NumericType                      _Ea;
          NumericType                      _Tref;
          NumericType                      _Troe_alpha;
          NumericType                      _Troe_T1;
          NumericType                      _Troe_T2;
          NumericType                      _Troe_T3;
          std::map<ParsingKey,std::string> _map;
          std::map<ParsingKey,std::string> _default_unit;

  };


  template <typename NumericType>
  inline
  ChemKinParser<NumericType>::ChemKinParser(const std::string &filename)
  {
    if(!_doc.open(filename.c_str()))
      {
        std::cerr << "ERROR: unable to load ChemKin file " << filename << std::endl;
        antioch_error();
      }

      _map[ParsingKey::REACTION_DATA] = "REACTIONS";
      _map[ParsingKey::FALLOFF_LOW]   = "LOW";

   // typically chemkin files list 
   //      pre-exponential parameters in (m3/kmol)^(m-1)/s
   //      activation energy in cal/mol, but we want it in K.
   //      power parameter without unit
   // if falloff, we need to know who's k0 and kinfty
   // if photochemistry, we have a cross-section on a lambda grid
   //      cross-section typically in cm2/nm (cross-section on a resolution bin, 
   //                                          if bin unit not given, it is lambda unit (supposed to anyway), and a warning message)
   //      lambda typically in nm, sometimes in ang, default considered here is nm
   //                         you can also have cm-1, conversion is done with
   //                         formulae nm = cm-1 * / * adapted factor
      _default_unit[ParsingKey::PREEXP]                = "m3/kmol";
      _default_unit[ParsingKey::POWER]                 = "";
      _default_unit[ParsingKey::ACTIVATION_ENERGY]    = "cal/mol";
      _default_unit[ParsingKey::BERTHELOT_COEFFICIENT] = "K-1";
      _default_unit[ParsingKey::TREF]                  = "K";
      _default_unit[ParsingKey::HV_LAMBDA]             = "nm";
      _default_unit[ParsingKey::HV_CROSS_SECTION]      = "cm2/nm";
      _default_unit[ParsingKey::EFFICIENCY]            = "";
      _default_unit[ParsingKey::TROE_F_ALPHA]          = "";
      _default_unit[ParsingKey::TROE_F_TS]             = "K";
      _default_unit[ParsingKey::TROE_F_TSS]            = "K";
      _default_unit[ParsingKey::TROE_F_TSSS]           = "K";
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::initialize()
  {
        //we start here
      _reaction_block = _doc.FirstChildElement("ctml");
      if (!_reaction_block) 
        {
          std::cerr << "ERROR:  no <ctml> tag found in input file"
                    << std::endl;
          antioch_error();
        }

      _reaction_block = _reaction_block->FirstChildElement(_map.at(ParsingKey::REACTION_DATA).c_str());

      _reaction = NULL;
      _rate_constant = NULL;

      return _reaction_block;
  }

  template <typename NumericType>
  inline
  ChemKinParser<NumericType>::~ChemKinParser()
  {
     return;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reaction()
  {
      antioch_assert(_reaction_block);
      _reaction = (_reaction)?
                        _reaction->NextSiblingElement(_map.at(ParsingKey::REACTION).c_str()):
                        _reaction_block->FirstChildElement(_map.at(ParsingKey::REACTION).c_str());

      _rate_constant = NULL;

      return _reaction;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant(const std::string & kinetics_model)
  {
      if(_reaction)
      {
        if(_rate_constant)
        {
           _rate_constant = _rate_constant->NextSiblingElement(kinetics_model.c_str());
        }else
        {
            antioch_assert(_reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str()));
            _rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str())->FirstChildElement(kinetics_model.c_str());
        }
      }else
      {
         _rate_constant = NULL;
      }

      return _rate_constant;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe()
  {
     if(_rate_constant)
     {
         _Troe = _rate_constant->FirstChildElement(_map[ParsingKey::TROE_FALLOFF].c_str());
     }else
     {
        _Troe = NULL;
     }

     return _Troe;
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_id() const
  {
    std::stringstream id;
    id << _reaction->IntAttribute(_map.at(ParsingKey::ID).c_str());
    return id.str();
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_equation() const
  {
     return _reaction->FirstChildElement(_map.at(ParsingKey::EQUATION).c_str())->GetText();
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_chemical_process() const
  {
    const char * chem_proc = _reaction->Attribute(_map.at(ParsingKey::CHEMICAL_PROCESS).c_str());

    return (chem_proc)?
                  std::string(chem_proc):
                  std::string();
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reaction_reversible() const
  {
     return (_reaction->Attribute(_map.at(ParsingKey::REVERSIBLE).c_str()))?
                      (std::string(_reaction->Attribute(_map.at(ParsingKey::REVERSIBLE).c_str())) == std::string("no"))?false:true //explicit
                      :
                      true; //default
  }

  template <typename NumericType>
  inline
  const std::string ChemKinParser<NumericType>::reaction_kinetics_model(const std::vector<std::string> &kinetics_models) const
  {
     unsigned int imod(0);
     tinyxml2::ChemKinElement * rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str())->FirstChildElement(kinetics_models[imod].c_str());
     while(!rate_constant)
     {
       if(imod == kinetics_models.size() - 1)
       {
                std::cerr << "Could not find a suitable kinetics model.\n"
                          << "Implemented kinetics models are:\n";

                for(unsigned int m = 0; m < kinetics_models.size(); m++)
                std::cerr << "  " << kinetics_models[m] << "\n";

                std::cerr << "See Antioch documentation for more details."
                          << std::endl;
                antioch_not_implemented();
       }
       imod++;
       rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str())->FirstChildElement(kinetics_models[imod].c_str());
     }

     return kinetics_models[imod];
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::reactants_pairs(std::vector<std::pair<std::string,int> > & reactants_pair) const
  {
     tinyxml2::ChemKinElement* reactants = _reaction->FirstChildElement(_map.at(ParsingKey::REACTANTS).c_str());
     return this->molecules_pairs(reactants, reactants_pair);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::products_pairs(std::vector<std::pair<std::string,int> > & products_pair) const
  {
     tinyxml2::ChemKinElement* products = _reaction->FirstChildElement(_map.at(ParsingKey::PRODUCTS).c_str());
     return this->molecules_pairs(products, products_pair);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::molecules_pairs(tinyxml2::ChemKinElement * molecules, std::vector<std::pair<std::string,int> > & molecules_pairs) const
  {
     bool out(true);
     if(molecules)
     {

       std::vector<std::string> mol_pairs;
     
       // Split the reactant string on whitespace. If no entries were found,
       // there is no whitespace - and assume then only one reactant is listed.
       if( !SplitString(std::string(molecules->GetText()),
                        " ",
                        mol_pairs,
                        /* include_empties = */ false)     )
         {
            mol_pairs.push_back(molecules->GetText());
         }

        for( unsigned int p=0; p < mol_pairs.size(); p++ )
          {
             std::pair<std::string,int> pair( split_string_int_on_colon(mol_pairs[p]) );

             molecules_pairs.push_back(pair);
          }
      }else
      {
         out = false;
      }

      return out;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::is_k0(unsigned int nrc, const std::string & kin_model) const
  {
      bool k0(false);
      if(_rate_constant->Attribute("name"))
      {
         if(std::string(_rate_constant->Attribute("name")) == "k0")k0 = true;
      }else if(nrc == 0) // if we're indeed at the first reading
      {
         antioch_assert(_rate_constant->NextSiblingElement(kin_model.c_str()));
         if(!_rate_constant->NextSiblingElement(kin_model.c_str())->Attribute("name")) // and the next doesn't have a name
         {
               k0 = true;
         }else
         {
             if(std::string(_rate_constant->NextSiblingElement(kin_model.c_str())->Attribute("name")) == "k0")k0 = false; 
         }
      }
      return k0;
  }

  template <typename NumericType>
  inline
  unsigned int ChemKinParser<NumericType>::where_is_k0(const std::string & kin_model) const
  {
      antioch_assert(!_rate_constant); //should be done exterior to rate constant block
      antioch_assert(_reaction);       //should be done interior to reaction block

      tinyxml2::ChemKinElement * rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str());
      antioch_assert(rate_constant);
      rate_constant = rate_constant->FirstChildElement(kin_model.c_str());
      unsigned int k0(0);
      if(rate_constant->NextSiblingElement()->Attribute("name"))
      {
        if(std::string(rate_constant->NextSiblingElement()->Attribute("name")) == "k0")k0=1;
      }

      return k0;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::get_parameter(const tinyxml2::ChemKinElement * ptr, const std::string & par, NumericType & par_value, std::string & par_unit) const
  {
     antioch_assert(ptr);

     bool out(false);
     par_unit.clear();
     if(ptr->FirstChildElement(par.c_str()))
     {
        par_value = std::atof(ptr->FirstChildElement(par.c_str())->GetText());
        if(ptr->FirstChildElement(par.c_str())->Attribute(_map.at(ParsingKey::UNIT).c_str()))
            par_unit = ptr->FirstChildElement(par.c_str())->Attribute(_map.at(ParsingKey::UNIT).c_str());
        out = true;
     }

     return out;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::get_parameter(const tinyxml2::ChemKinElement * ptr, const std::string & par, std::vector<NumericType> & par_values, std::string & par_unit) const
  {
     antioch_assert(ptr);

     bool out(false);
     par_unit.clear();
     if(ptr->FirstChildElement(par.c_str()))
     {
        std::vector<std::string> values;
        if(!SplitString(ptr->FirstChildElement(par.c_str())->GetText()," ",values,false))
            values.push_back(ptr->FirstChildElement(par.c_str())->GetText());

        par_values.resize(values.size());
        for(unsigned int i = 0; i < values.size(); i++)
        {
          par_values[i] =  std::atof(values[i].c_str());
        }
        if(ptr->FirstChildElement(par.c_str())->Attribute(_map.at(ParsingKey::UNIT).c_str()))
            par_unit = ptr->FirstChildElement(par.c_str())->Attribute(_map.at(ParsingKey::UNIT).c_str());
        out = true;
     }

     return out;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_preexponential_parameter(NumericType & A, std::string & A_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::PREEXP);
      return this->get_parameter(_rate_constant,_map.at(ParsingKey::PREEXP).c_str(),A,A_unit);
  }


  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_power_parameter(NumericType & b, std::string & b_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::POWER);
      return this->get_parameter(_rate_constant,_map.at(ParsingKey::POWER).c_str(),b,b_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_activation_energy_parameter(NumericType & Ea, std::string & Ea_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::ACTIVATION_ENERGY);
      return this->get_parameter(_rate_constant,_map.at(ParsingKey::ACTIVATION_ENERGY).c_str(),Ea,Ea_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_Berthelot_coefficient_parameter(NumericType & D, std::string & D_unit, std::string & def_unit) const
  {
      return false; //not a chemkin model
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_lambda_parameter(std::vector<NumericType> & lambda, std::string & lambda_unit, std::string & def_unit) const
  {
      return false; //not a chemkin model
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_cross_section_parameter(std::vector<NumericType> & sigma, std::string & sigma_unit, std::string & def_unit) const
  {
      return false; //not a chemkin model
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::rate_constant_Tref_parameter(NumericType & Tref, std::string & Tref_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::TREF);
      return this->get_parameter(_rate_constant,_map.at(ParsingKey::TREF).c_str(),Tref,Tref_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::verify_Kooij_in_place_of_Arrhenius() const
  {
     bool out(false);
     tinyxml2::ChemKinElement * rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str());
     antioch_assert(rate_constant->FirstChildElement("Arrhenius"));
     rate_constant = rate_constant->FirstChildElement("Arrhenius");
     if(rate_constant->FirstChildElement(_map.at(ParsingKey::POWER).c_str()))
     {
        if(std::atof(rate_constant->FirstChildElement(_map.at(ParsingKey::POWER).c_str())->GetText()) != 0.) //not a very good test
                out = true;
     }
     
     return out;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::efficiencies(std::vector<std::pair<std::string,double> > & par_values) const
  {
     bool out = false;
     if(_reaction)
     {
       tinyxml2::ChemKinElement * rate_constant = _reaction->FirstChildElement(_map.at(ParsingKey::KINETICS_MODEL).c_str());
       if(rate_constant)
       {
          if(rate_constant->FirstChildElement(_map.at(ParsingKey::EFFICIENCY).c_str()))
          {
            std::vector<std::string> values;
            SplitString(std::string((rate_constant->FirstChildElement(_map.at(ParsingKey::EFFICIENCY).c_str())->GetText())?rate_constant->FirstChildElement(_map.at(ParsingKey::EFFICIENCY).c_str())->GetText():""),
                                     " ",values,false);
            for(unsigned int i = 0; i < values.size(); i++)
            {
               par_values.push_back(split_string_double_on_colon (values[i]));
            }
            out = true;
          }
        }
     }
     return out;
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_alpha_parameter(NumericType & alpha, std::string & alpha_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::TROE_F_ALPHA);
      return this->get_parameter(_Troe,_map.at(ParsingKey::TROE_F_ALPHA),alpha,alpha_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T1_parameter(NumericType & T1, std::string & T1_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::TROE_F_TS);
      return this->get_parameter(_Troe,_map.at(ParsingKey::TROE_F_TS),T1,T1_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T2_parameter(NumericType & T2, std::string & T2_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::TROE_F_TSS);
      return this->get_parameter(_Troe,_map.at(ParsingKey::TROE_F_TSS),T2,T2_unit);
  }

  template <typename NumericType>
  inline
  bool ChemKinParser<NumericType>::Troe_T3_parameter(NumericType & T3, std::string & T3_unit, std::string & def_unit) const
  {
      def_unit = _default_unit.at(ParsingKey::TROE_F_TSSS);
      return this->get_parameter(_Troe,_map.at(ParsingKey::TROE_F_TSSS),T3,T3_unit);
  }


}//end namespace Antioch

#endif
