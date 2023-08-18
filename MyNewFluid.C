//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MyNewFluid.h"
#include "SinglePhaseFluidProperties.h"


registerMooseObject("PorousFlowApp", MyNewFluid);
registerMooseObject("PorousFlowApp", ADMyNewFluid);

template <bool is_ad>
InputParameters
MyNewFluidTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowFluidPropertiesBaseTempl<is_ad>::validParams();
  params.addRequiredParam<UserObjectName>("fp", "The name of the user object for fluid properties");
    params.addRequiredParam<RealVectorValue>("gravity",
                                           "Gravitational acceleration vector downwards (m/s^2)");
  params.addClassDescription("This Material calculates fluid properties at the quadpoints or nodes "
                             "for a single component fluid");
  return params;
}

template <bool is_ad>
MyNewFluidTempl<is_ad>::MyNewFluidTempl(
    const InputParameters & parameters)
  : PorousFlowFluidPropertiesBaseTempl<is_ad>(parameters),
    _darcy_velocity_old(this->template getMaterialPropertyOldByName<std::vector<RealVectorValue>>("PorousFlow_darcy_velocity_qp")),
    _fp(this->template getUserObject<SinglePhaseFluidProperties>("fp"))
{
}

template <bool is_ad>
void
MyNewFluidTempl<is_ad>::initQpStatefulProperties()
{

  if (_compute_rho_mu)
  {
    (*_density)[_qp] = _fp.rho_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                        _temperature[_qp] + _t_c2k);


    (*_viscosity)[_qp] = _fp.mu_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                         _temperature[_qp] + _t_c2k) /
                         _pressure_to_Pascals / _time_to_seconds;
  }

  if (_compute_internal_energy)
    (*_internal_energy)[_qp] = _fp.e_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                              _temperature[_qp] + _t_c2k);

  if (_compute_enthalpy)
    (*_enthalpy)[_qp] = _fp.h_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                       _temperature[_qp] + _t_c2k);
}


template <bool is_ad>
void
MyNewFluidTempl<is_ad>::computeQpProperties()
{
  Real velocity_abs;
  if (_compute_rho_mu)
  {
    if (is_ad)
    {
      GenericReal<is_ad> rho, mu;
      _fp.rho_mu_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                          _temperature[_qp] + _t_c2k,
                          rho,
                          mu);

      (*_density)[_qp] = rho;
      (*_viscosity)[_qp] = mu / _pressure_to_Pascals / _time_to_seconds;
    }
    else
    {
      // Density and viscosity, and derivatives wrt pressure and temperature
      Real rho, drho_dp, drho_dT, mu, dmu_dp, dmu_dT;
      _fp.rho_mu_from_p_T(MetaPhysicL::raw_value(_porepressure[_qp][_phase_num]) *
                              _pressure_to_Pascals,
                          MetaPhysicL::raw_value(_temperature[_qp]) + _t_c2k,
                          rho,
                          drho_dp,
                          drho_dT,
                          mu,
                          dmu_dp,
                          dmu_dT);
      (*_density)[_qp] = rho;
      (*_ddensity_dp)[_qp] = drho_dp * _pressure_to_Pascals;
      (*_ddensity_dT)[_qp] = drho_dT;
      velocity_abs = std::sqrt(MetaPhysicL::raw_value(_darcy_velocity_old[_qp][0]).norm() * MetaPhysicL::raw_value(_darcy_velocity_old[_qp][0]).norm());
      (*_viscosity)[_qp] = velocity_abs;
      (*_dviscosity_dp)[_qp] = 0 / _time_to_seconds;
      (*_dviscosity_dT)[_qp] = 0 / _pressure_to_Pascals / _time_to_seconds;

    }
  }

  if (_compute_enthalpy)
  {
    if (is_ad)
      (*_enthalpy)[_qp] = _fp.h_from_p_T(_porepressure[_qp][_phase_num] * _pressure_to_Pascals,
                                         _temperature[_qp] + _t_c2k);
    else
    {
      // Enthalpy and derivatives wrt pressure and temperature at the qps
      Real h, dh_dp, dh_dT;
      _fp.h_from_p_T(MetaPhysicL::raw_value(_porepressure[_qp][_phase_num]) * _pressure_to_Pascals,
                     MetaPhysicL::raw_value(_temperature[_qp]) + _t_c2k,
                     h,
                     dh_dp,
                     dh_dT);
      (*_enthalpy)[_qp] = h;
      (*_denthalpy_dp)[_qp] = dh_dp * _pressure_to_Pascals;
      (*_denthalpy_dT)[_qp] = dh_dT;
    }
  }
}

template class MyNewFluidTempl<false>;
template class MyNewFluidTempl<true>;

