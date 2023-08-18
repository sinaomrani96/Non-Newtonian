
[Mesh]
  [gen]
  type = GeneratedMeshGenerator
  dim = 2
  nx = 10
  ny = 10
  xmax = 1
  ymax = 1
  []
   [efluent]
     input = gen
     type = SubdomainBoundingBoxGenerator
     block_id = 1
     bottom_left = '0.9 0 0'
     top_right = '1 0.1 0'
   []
[]

[GlobalParams]
  PorousFlowDictator = 'dictator'
  gravity = '0 -9.8 0'
[]

[AuxVariables]
  [saturation_gas]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_liq]
    order = CONSTANT
    family = MONOMIAL
  []
  [x00] 
    # order = CONSTANT
    # family = MONOMIAL
    initial_condition = 0
  []
  [x10]
    # order = CONSTANT
    # family = MONOMIAL
    initial_condition = 1
  []
  [x01]
    # order = CONSTANT
    # family = MONOMIAL
  []
  [x11]
    # order = CONSTANT
    # family = MONOMIAL
    initial_condition = 0
  []
  [x02]
    # order = CONSTANT
    # family = MONOMIAL
    initial_condition = 1
  []
  [x12]
    # order = CONSTANT
    # family = MONOMIAL
    initial_condition = 0
  []
  [mf1]
    order = CONSTANT
    family = MONOMIAL
  []
  [mf2]
    order = CONSTANT
    family = MONOMIAL
  []
  [mf3]
    order = CONSTANT
    family = MONOMIAL
  []
  [mf4]
    order = CONSTANT
    family = MONOMIAL
  []
  [mf5]
    order = CONSTANT
    family = MONOMIAL
  []
  [mf6]
    order = CONSTANT
    family = MONOMIAL
  []
  [d1]
    order = FIRST
    family = LAGRANGE
  []
  [d2]
    order = FIRST
    family = LAGRANGE
  []
  [d3]
    order = FIRST
    family = LAGRANGE
  []
  [d4]
    order = FIRST
    family = LAGRANGE
  []
  [d5]
    order = FIRST
    family = LAGRANGE
  []
  [d6]
    order = FIRST
    family = LAGRANGE
  []
  [d7]
    order = FIRST
    family = LAGRANGE
  []
  [d8]
    order = FIRST
    family = LAGRANGE
  []
  [d9]
    order = FIRST
    family = LAGRANGE
  []
  [aperture1]
      family = MONOMIAL
      order = CONSTANT
  []
  [aperture2]
    family = MONOMIAL
    order = CONSTANT
  []
  [density]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosityn]
    family = MONOMIAL
    order = CONSTANT
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
  [gas_velocityx]
    family = MONOMIAL
    order = CONSTANT
  []
  [gas_velocityy]
    family = MONOMIAL
    order = CONSTANT
  []
  [liq_velocityx]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0
  []
  [liq_velocityy]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0
  []
  [test]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'timestep_end'
  []
  [viscosity_gas]
    type = PorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
    execute_on = 'timestep_end'
  []
  [viscosity_liq]
    type = PorousFlowPropertyAux
    variable = viscosity_liq
    property = viscosity
    phase = 0
    execute_on = 'timestep_end'
  []
  [mf1]
    type = PorousFlowPropertyAux
    variable = mf1
    property = mass_fraction
    phase = 0
    fluid_component = 0
    execute_on = 'timestep_end'
  []
  [mf2]
    type = PorousFlowPropertyAux
    variable = mf2
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = 'timestep_end'
  []
  [mf3]
    type = PorousFlowPropertyAux
    variable = mf3
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'timestep_end'
  []
  [mf4]
    type = PorousFlowPropertyAux
    variable = mf4
    property = mass_fraction
    phase = 1
    fluid_component = 1
    execute_on = 'timestep_end'
  []
  [mf5]
    type = PorousFlowPropertyAux
    variable = mf5
    property = mass_fraction
    phase = 0
    fluid_component = 2
    execute_on = 'timestep_end'
  []
  [mf6]
    type = PorousFlowPropertyAux
    variable = mf6
    property = mass_fraction
    phase = 1
    fluid_component = 2
    execute_on = 'timestep_end'
  []
  [Density]
    type = PorousFlowPropertyAux
    variable = density
    property = density
    phase = 0
    #fluid_component = 2
    execute_on = 'timestep_end'
  []
  [Porosityn]
    type = PorousFlowPropertyAux
    variable = porosityn
    property = porosity
    execute_on = 'TIMESTEP_BEGIN'
  []  
  [permeability]
    type = PorousFlowPropertyAux
    variable = permeability
    property = permeability
    execute_on = 'TIMESTEP_BEGIN'
  []  
  [gas_velocityx]
    type = PorousFlowDarcyVelocityComponent
    variable = gas_velocityx
    component = x
    fluid_phase = 1
    # execute_on = 'LINEAR TIMESTEP_END'
    execute_on = 'INITIAL TIMESTEP_BEGIN'

  []   
  [gas_velocityy]
    type = PorousFlowDarcyVelocityComponent
    variable = gas_velocityy
    component = y
    fluid_phase = 1
    # execute_on = 'LINEAR TIMESTEP_END'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [liq_velocityx]
    type = PorousFlowDarcyVelocityComponent
    variable = liq_velocityx
    component = x
    fluid_phase = 0
    execute_on = 'LINEAR TIMESTEP_END'
  []   
  [liq_velocityy]
    type = PorousFlowDarcyVelocityComponent
    variable = liq_velocityy
    component = y
    fluid_phase = 0
    execute_on = 'LINEAR TIMESTEP_END'
  []
[]

[Variables]
  [pwater]
    initial_condition = 20e6
    #  scaling = 1e4
    # scaling = 1e-2
  []
  [satg]
    initial_condition = 0
  #  scaling = 1e-2
  []
  [tracer]
    initial_condition =  0
  #  scaling = 1e-2
  []
[]


[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 2
    variable = pwater
    save_in = d1
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 2
    variable = pwater
    save_in = d2
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = tracer
    save_in = d4
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = tracer
    save_in = d5
  []

  [mass2]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = satg
    save_in = d7
  []
  [flux2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = satg
    save_in = d8
  []
 
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater tracer satg'
    number_fluid_phases = 2
    number_fluid_components = 3
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [pc]
    type = PorousFlowCapillaryPressureBC
    lambda = 2.0
    pe = 1e4
    # type = PorousFlowCapillaryPressureConst
    # pc = 0
  []

[]

[FluidProperties]
  [true_water]
    type = Water97FluidProperties
  []
  # [tabulated_water]
  #   type = TabulatedFluidProperties
  #   fp = true_water
  #   temperature_min = 275
  #   pressure_max = 1E8
  #   interpolated_properties = 'density viscosity enthalpy internal_energy'
  #   fluid_property_file = water97_tabulated_11.csv
  # []
  [true_co2]
    type = CO2FluidProperties
  []
  [tabulated_co2]
    type = TabulatedFluidProperties
    fp = true_co2
    temperature_min = 275
    pressure_max = 1E8
    interpolated_properties = 'density viscosity enthalpy internal_energy'
    fluid_property_file = co2_tabulated_11.csv
  []
  [fp1]
    type = SimpleFluidProperties
    density0 = 1000
    bulk_modulus = 2e+109
    viscosity = 0.001
    execute_on = 'INITIAL TIMESTEP_END'
  []

[]


[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '323'
  []
  [water]
    type = MyNewFluid
    fp = fp1
    phase = 0
  []
  [co2]
    type = PorousFlowSingleComponentFluid
    fp = true_co2
    phase = 1
  []
  [velocity1]
    type = PorousFlowDarcyVelocityMaterial
    outputs = 'exodus console'
  []


  [massfrac1]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'x00 tracer x10 x11'
  []
  [saturation_calculator]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = satg
    capillary_pressure = pc
  []

  [por]
    type = PorousFlowPorosityConst
    porosity = 0.2
  []
  [per]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0 0 1e-13 0 0 0 1e-13'
  []

  [relperm_water]
    type = PorousFlowRelativePermeabilityBC
    nw_phase = false
    lambda = 1
    s_res = 0.1
    sum_s_res = 0.2
    phase = 0
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    nw_phase = true
    lambda = 2
    s_res = 0.1
    sum_s_res = 0.2
    phase = 1
  []

[]

[Functions]
  [tracerfunc]
    type = PiecewiseConstant
    x = '0       5000     10000    50000'
    y = '0  0    0    0'
    direction = LEFT_INCLUSIVE
  []
  [gasfunc]
    type = PiecewiseConstant
    x = '0       5000     10000    50000'
    y = '0    0     -1e-6     -1e-6'
    direction = LEFT_INCLUSIVE
  []
  [waterfunc]
      type = PiecewiseConstant
      x = '0       5000     10000    50000'
      y = '-1e-6    -1e-6    0     0'
      direction = LEFT_INCLUSIVE
  []

[]




[BCs]

  [injection2]
    type = PorousFlowSink
    variable = tracer
    flux_function = tracerfunc
    boundary = 'left'
  []
  [injection1]
    type = PorousFlowSink
    variable = satg
    flux_function = gasfunc
    boundary = 'left'
  []
  [injection3]
    type = PorousFlowSink
    variable = pwater
    flux_function = waterfunc
    boundary = 'left'
  []

  [outp]
    type = DirichletBC
    variable = pwater
    value = 20e6
    boundary = right
  []  

[]

 [Postprocessors]
   [tracer_mass]
     type = PorousFlowFluidMass
     phase = 0
     fluid_component = 1
     execute_on = TIMESTEP_END
     outputs = 'csv'
   []
   [sg]
     type = ElementAverageValue
     variable = satg
     execute_on = TIMESTEP_END
     outputs = 'csv'
   []
   [water_mass]
     type = PorousFlowFluidMass
     phase = 0
     fluid_component = 2
     execute_on = TIMESTEP_END
     outputs = 'csv'
   []
   [gas_mass]
     type = PorousFlowFluidMass
     phase = 1
     fluid_component = 0
     execute_on = TIMESTEP_END
     outputs = 'csv'
   []
   [density_water]
     type = ElementAverageValue
     variable = density
     execute_on = 'timestep_end'
     outputs = 'csv'
   []
   [porosity]
     type = ElementAverageValue
     variable = porosityn
     execute_on = 'INITIAL timestep_end'
     outputs = 'csv console'
   []
   [permeability]
    type = ElementAverageValue
    variable = permeability
    execute_on = 'INITIAL timestep_end'
    outputs = 'csv console'
   []
   [massfractiontracer]
     type = ElementAverageValue
     variable = mf3
     execute_on = 'timestep_end'
     outputs = 'csv'
   []
   [massfractionwater]
     type = ElementAverageValue
     variable = mf5
     execute_on = 'timestep_end'
     outputs = 'csv'
   []
   [vis_liq]
     type = ElementAverageValue
     variable = viscosity_liq
     execute_on = 'timestep_end'
     outputs = 'csv'
   []
   [vis_gas]
     type = ElementAverageValue
     variable = viscosity_gas
     execute_on = 'timestep_end'
     outputs = 'csv'
   []
   [efluent]
     type = PorousFlowFluidMass
     phase = 1
     fluid_component = 0
     execute_on = TIMESTEP_END
     block = '1'
     outputs = 'csv'
   []

  [pvx]
    type = ElementAverageValue
    variable = liq_velocityx
    execute_on = 'timestep_end'
    outputs = 'csv'
  []
  [pvy]
    type = ElementAverageValue
    variable = liq_velocityy
    execute_on = 'timestep_end'
    outputs = 'csv'
  []

[]


[Preconditioning]
  active = 'bounded'
  [basic]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres asm lu NONZERO 2'
  []
  [preferred]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'lu mumps'
  []
  [bounded]
    # must use --use-petsc-dm command line argument
      type = SMP
      full = true
      petsc_options = '-snes_linesearch_monitor'
      petsc_options_iname = '-snes_type -pc_factor_shift_type -snes_linesearch_type'
      petsc_options_value = 'vinewtonssls nonzero l2'
    []
[]


[Executioner]

  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 10
  dtmax = 10
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
  []

# controls for linear iterations
  l_max_its = 40
  l_abs_tol = 1e-12




# controls for nonlinear iterations
  nl_max_its = 40
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-11
[]

[Outputs]
  exodus = true
  csv = true
[]

[Debug]
  show_var_residual_norms = true
 # show_top_residuals = 2
  # show_material_props = true
  show_execution_order = ALWAYS
  # show_actions = true
  # show_action_dependencies = true
  # show_parser = true
  # show_functors = true
  # show_reporters = true
[]
