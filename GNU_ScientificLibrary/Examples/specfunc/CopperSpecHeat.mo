within GNU_ScientificLibrary.Examples.specfunc;

model CopperSpecHeat "Copper Specific Heat"
  // example by tommy.burch@physik.uni-r.de
  extends Modelica.Icons.Example;
  import const = Modelica.Constants;
  final constant Real m_e = 9.10938e-31 "mass of the electron in kg";
  final constant Real k_dens = 2*(2*const.pi*m_e*const.k/const.h^2)^1.5;
  final constant Real k_pres = k_dens*const.k;
  Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
    Placement(transformation(origin = {10, -10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T (displayUnit = "K")= 373)  annotation(
    Placement(transformation(origin = {50, -10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Constant Gc_const(k = 1)  annotation(
    Placement(transformation(origin = {-10, 30}, extent = {{-10, -10}, {10, 10}})));
  Blocks.Physics.DebyeFermiDiracSolid copperLump(T(start = 1, fixed = true))  annotation(
    Placement(transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(fixedTemperature.port, convection.fluid) annotation(
    Line(points = {{40, -10}, {20, -10}}, color = {191, 0, 0}));
  connect(Gc_const.y, convection.Gc) annotation(
    Line(points = {{2, 30}, {10, 30}, {10, 0}}, color = {0, 0, 127}));
  connect(copperLump.port, convection.solid) annotation(
    Line(points = {{-30, -10}, {0, -10}}, color = {191, 0, 0}));
  annotation(
    __OpenModelica_simulationFlags(lv = "LOG_STATS", s = "dassl", variableFilter = ".*"),
    experiment(StartTime = 0, StopTime = 120, Tolerance = 1e-06, Interval = 0.001),
    Documentation(info = "<html><head></head><body>Warming a lump of copper via convection. The heat-transfer modeling is rather simplistic...&nbsp;<br><div><br></div><div><img src=\"modelica://GNU_ScientificLibrary/Examples/specfunc/CopperSpecHeat_T_vs_t.png\"></div><div><br></div><div>But the main feature here -- namely, the temperature-dependent heat capacity -- is still thereby revealed (see below). Note the c<sub>v</sub>~T<sup>3</sup>-behavior for low temperatures, accurately predicted in the Debye model (for high-T, this correctly turns over to an approach to a constant). Note also the relatively small contribution from the conduction electrons (well approximated by a gas of degenerate fermions), which, despite not holding much heat, contribute significantly to its quick transport. (See also the Physics.DebyeEoS and Physics.FermiDiracEoS blocks.)</div><div><br><img src=\"modelica://GNU_ScientificLibrary/Examples/specfunc/CopperSpecHeat_cv_vs_T.png\"></div><div><br></div></body></html>"),
    Diagram(coordinateSystem(extent = {{-40, 40}, {60, -20}})));
end CopperSpecHeat;