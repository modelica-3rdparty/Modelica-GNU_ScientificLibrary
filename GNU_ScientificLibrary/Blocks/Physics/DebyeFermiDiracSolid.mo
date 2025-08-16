within GNU_ScientificLibrary.Blocks.Physics;

model DebyeFermiDiracSolid "Lumped thermal element with variable heat capacity"
  // contact: tommy.burch@physik.uni-r.de
  extends GNU_ScientificLibrary.Icons.Block;
  parameter Modelica.Units.SI.Mass mass = 1 "kg";
  parameter Real T_D = 343 "K ; Debye Temperature";
  parameter Real eff_am = 63.546 "Da ; Effective atomic mass";
  parameter Real mdens = 8935 "kg/m^3 ; Mass density";
  parameter Integer N_e = 1 "# conduction e-'s per atom";
  Real cp;
  Modelica.Units.SI.HeatCapacity C "Heat capacity of element (= cp*m)";
  Modelica.Units.SI.Temperature T(start = 293.15, displayUnit = "K") "Temperature of element";
  Modelica.Units.SI.TemperatureSlope der_T(start = 0) "Time derivative of temperature (= der(T))";
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(
    Placement(transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  DebyeEoS debyeEoS(T_D = T_D, eff_am = eff_am, mdens = mdens)  annotation(
    Placement(transformation(origin = {-50, 20}, extent = {{-10, -10}, {10, 10}})));
  FermiDiracEoS fermiDiracEoS(eff_am = eff_am, N_e = N_e, mdens = mdens)  annotation(
    Placement(transformation(origin = {-50, -20}, extent = {{-10, 10}, {10, -10}}, rotation = -0)));
  Modelica.Blocks.Math.Add add annotation(
    Placement(transformation(origin = {-10, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput eD annotation(
    Placement(transformation(origin = {10, 26}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput eFD_minus_e0 annotation(
    Placement(transformation(origin = {10, -26}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput c_v annotation(
    Placement(transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  cp = add.y;
  C = cp*mass;
  T = port.T;
  debyeEoS.T = T;
  fermiDiracEoS.T = T;
  der_T = der(T);
  C*der(T) = port.Q_flow;
  connect(debyeEoS.c_v, add.u1) annotation(
    Line(points = {{-38, 14}, {-30, 14}, {-30, 6}, {-22, 6}}, color = {0, 0, 127}));
  connect(fermiDiracEoS.c_v, add.u2) annotation(
    Line(points = {{-39, -14}, {-30, -14}, {-30, -6}, {-22, -6}}, color = {0, 0, 127}));
  connect(debyeEoS.edens, eD) annotation(
    Line(points = {{-38, 26}, {10, 26}}, color = {0, 0, 127}));
  connect(fermiDiracEoS.e_minus_e0, eFD_minus_e0) annotation(
    Line(points = {{-39, -26}, {10, -26}}, color = {0, 0, 127}));
  connect(add.y, c_v) annotation(
    Line(points = {{2, 0}, {30, 0}}, color = {0, 0, 127}));
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Text(extent = {{-150, 110}, {150, 70}}, textString = "%name", textColor = {0, 0, 255}), Polygon(points = {{0, 67}, {-20, 63}, {-40, 57}, {-52, 43}, {-58, 35}, {-68, 25}, {-72, 13}, {-76, -1}, {-78, -15}, {-76, -31}, {-76, -43}, {-76, -53}, {-70, -65}, {-64, -73}, {-48, -77}, {-30, -83}, {-18, -83}, {-2, -85}, {8, -89}, {22, -89}, {32, -87}, {42, -81}, {54, -75}, {56, -73}, {66, -61}, {68, -53}, {70, -51}, {72, -35}, {76, -21}, {78, -13}, {78, 3}, {74, 15}, {66, 25}, {54, 33}, {44, 41}, {36, 57}, {26, 65}, {0, 67}}, lineColor = {160, 160, 164}, fillColor = {192, 192, 192}, fillPattern = FillPattern.Solid), Polygon(points = {{-58, 35}, {-68, 25}, {-72, 13}, {-76, -1}, {-78, -15}, {-76, -31}, {-76, -43}, {-76, -53}, {-70, -65}, {-64, -73}, {-48, -77}, {-30, -83}, {-18, -83}, {-2, -85}, {8, -89}, {22, -89}, {32, -87}, {42, -81}, {54, -75}, {42, -77}, {40, -77}, {30, -79}, {20, -81}, {18, -81}, {10, -81}, {2, -77}, {-12, -73}, {-22, -73}, {-30, -71}, {-40, -65}, {-50, -55}, {-56, -43}, {-58, -35}, {-58, -25}, {-60, -13}, {-60, -5}, {-60, 7}, {-58, 17}, {-56, 19}, {-52, 27}, {-48, 35}, {-44, 45}, {-40, 57}, {-58, 35}}, fillColor = {160, 160, 164}, fillPattern = FillPattern.Solid), Text(extent = {{-69, 7}, {71, -24}}, textString = "%C")}),
    Documentation(info = "<html><head></head><body><blockquote><p>Similar to the HeatCapacitor block in the MSL, except that the specific heat capacity is allowed to vary according to the temperature: c_v(T)=c_p(T) (good approximation for solids). The heat capacity is the sum of those determined via the Debye model (lattice phonons) and the Fermi-Dirac distribution (conduction electrons). This approach is fitting for many solids, both electrically-conducting (N<sub>e</sub>&gt;0) and electrically-insulating (N<sub>e</sub>=0).</p><p>Numerical outputs include the Debye-model energy density, the total specific heat capacity, and the Fermi-Dirac energy density above the zero-point (T=0) energy density.</p></blockquote>
</body></html>"),
  Diagram(coordinateSystem(extent = {{-100, 40}, {40, -40}})));
end DebyeFermiDiracSolid;