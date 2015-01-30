within IDEAS.Fluid.HeatExchangers.RadiantSlab.Examples;
model TabsValidation
  extends Modelica.Icons.Example;
  replaceable package Medium = IDEAS.Media.Water.Simple;
  parameter Integer nCap = 5;
  parameter Real Rlayer = Rtot/nCap;
  parameter Real Rtot = 0.1/1.8/A/2;
  parameter Real A = 20;
  parameter Real rConv = 1/11/A/2;
  parameter Real Clayer = 2*0.1*840*2100/nCap*A;

  IDEAS.Fluid.HeatExchangers.RadiantSlab.EmbeddedPipe embeddedPipe(
    redeclare package Medium = Medium,
    m_flow_nominal=12/3600*A,
    A_floor=A,
    nDiscr=1,
    m_flowMin=12/3600*A,
    T_start=293.15,
    redeclare
      IDEAS.Fluid.HeatExchangers.Examples.BaseClasses.RadSlaCha_ValidationEmpa
      RadSlaCha,
    R_c=Rtot + rConv,
    nParCir=1)
    annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
  Sources.MassFlowSource_T boundary(
    nPorts=1,
    redeclare package Medium = Medium,
    use_T_in=true,
    m_flow=12/3600*A)
    annotation (Placement(transformation(extent={{-60,-50},{-40,-30}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=273.15
         + 20)
    annotation (Placement(transformation(extent={{-40,80},{-20,100}})));

  Sources.Boundary_pT bou(nPorts=1, redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{80,-50},{60,-30}})));
  Sensors.TemperatureTwoPort senTem(redeclare package Medium = Medium,
      m_flow_nominal=1)
    annotation (Placement(transformation(extent={{26,-50},{46,-30}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor[nCap] thermalResistor(R=cat(
        1,
        fill(Rlayer/2, 1),
        fill(Rlayer, nCap - 1)))                                   annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,-10})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor[nCap] heatCapacitor(
                           each T(fixed=true, start=293.15), each C=Clayer)
    "2 parallel layers"
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-20,10})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor thermalResistor2(R=Rlayer/
        2)                                                         annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,30})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor conv(R=rConv)
    "h = 11 W/m2K" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={0,66})));
  Modelica.Blocks.Sources.CombiTimeTable Tsurface(table={{-99.18633116,293.1319823},
        {3590.687902,293.4706311},{10139.48326,294.1847965},{15971.8536,294.7194988},
        {22990.53528,295.1629236},{30481.60078,295.4617313},{37614.81056,295.6888252},
        {44270.99831,295.8263065},{52352.1868,295.9263253},{59482.18551,295.9912602},
        {69938.85095,296.0528652}}, smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative)
    annotation (Placement(transformation(extent={{-80,20},{-60,40}})));
  Modelica.Blocks.Sources.CombiTimeTable Tcore(smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative,
      table={{0,293.150},{3614.592521,294.678},{8026.957116,295.502},{15051.34736,
        296.234},{22071.81297,296.767},{30872.28075,297.191},{40144.77557,297.452},
        {47989.41537,297.606},{56189.77017,297.724},{65815.12572,297.805}})
    annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
  Modelica.Thermal.FluidHeatFlow.Examples.Utilities.DoubleRamp doubleRamp(
    startTime=0,
    interval=1e5,
    duration_1=1,
    duration_2=1,
    offset=293.15,
    height_1=10,
    height_2=-10)
    annotation (Placement(transformation(extent={{-100,-50},{-80,-30}})));
equation

  for i in 1:(nCap-1) loop
    connect(thermalResistor[i].port_b,thermalResistor[i+1].port_a);
      connect(heatCapacitor[i].port, thermalResistor[i].port_b) annotation (Line(
      points={{-10,10},{0,10},{0,0}},
      color={191,0,0},
      smooth=Smooth.None));
  end for;
  for i in 1:nCap loop
      connect(heatCapacitor[i].port, thermalResistor[i].port_b) annotation (Line(
      points={{-10,10},{0,10},{0,0}},
      color={191,0,0},
      smooth=Smooth.None));
  end for;

  connect(boundary.ports[1], embeddedPipe.port_a) annotation (Line(
      points={{-40,-40},{-10,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(embeddedPipe.port_b, senTem.port_a) annotation (Line(
      points={{10,-40},{26,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(senTem.port_b, bou.ports[1]) annotation (Line(
      points={{46,-40},{60,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(embeddedPipe.heatPortEmb[1], thermalResistor[1].port_a) annotation (
      Line(
      points={{0,-30},{0,-20}},
      color={191,0,0},
      smooth=Smooth.None));

  connect(thermalResistor[nCap].port_b, thermalResistor2.port_a) annotation (Line(
      points={{0,0},{0,20}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(fixedTemperature.port, conv.port_b) annotation (Line(
      points={{-20,90},{0,90},{5.55112e-16,76}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(conv.port_a, thermalResistor2.port_b) annotation (Line(
      points={{-5.55112e-16,56},{-5.55112e-16,48},{0,48},{0,40}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(doubleRamp.y, boundary.T_in) annotation (Line(
      points={{-79,-40},{-70,-40},{-70,-36},{-62,-36}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
    experiment(StopTime=72000),
    __Dymola_experimentSetupOutput,
    __Dymola_Commands(file=
          "modelica://IDEAS/Resources/Scripts/Dymola/Fluid/HeatExchangers/RadiantSlab/Examples/TabsValidation.mos"
        "Simulate and plot"),
    Documentation(info="<html>
<p>Validation of concrete core activation (TABS) model based on Koschenz.</p>
</html>", revisions="<html>
<ul>
<li>
January 30, 2015 by Filip Jorissen:<br/>
First implementation.
</li>
</ul>
</html>"));
end TabsValidation;
