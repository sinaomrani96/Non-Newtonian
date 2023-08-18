//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowFluidPropertiesBase.h"


class SinglePhaseFluidProperties;

/**
 * General single component fluid material. Provides quadpoint density, viscosity,
 * internal energy, enthalpy and derivatives wrt pressure and temperature
 * for a fluid defined in the FluidProperties module
 */
template <bool is_ad>
class MyNewFluidTempl : public PorousFlowFluidPropertiesBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  MyNewFluidTempl(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;


  const MaterialProperty<std::vector<RealVectorValue>> & _darcy_velocity_old;

  const SinglePhaseFluidProperties & _fp;


  usingPorousFlowFluidPropertiesMembers;
};

typedef MyNewFluidTempl<false> MyNewFluid;
typedef MyNewFluidTempl<true> ADMyNewFluid;
