#!/bin/bash
set -e

###############################################################################
# physicalProperties — steel for all solid regions
###############################################################################
mkdir -p constant/solid_pipes constant/solid_shell constant/solid_baffles

cat > constant/solid_pipes/physicalProperties << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      physicalProperties;
}
thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIsoSolid;
    thermo          eConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleInternalEnergy;
}
mixture
{
    specie          { molWeight 55.845; }
    transport       { kappa 50; }
    thermodynamics  { Cv 500; hf 0; }
    equationOfState { rho 8000; }
}
EOF

for r in solid_shell solid_baffles; do
    mkdir -p constant/$r
    ln -sf $(pwd)/constant/solid_pipes/physicalProperties constant/$r/physicalProperties
done

###############################################################################
# physicalProperties — water (rhoConst) for fluid regions
###############################################################################
for r in fluid_1_tubeside fluid_2_shellside; do
    mkdir -p constant/$r
    cat > constant/$r/physicalProperties << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      physicalProperties;
}
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie          { molWeight 18.0; }
    equationOfState { rho 1000; }
    thermodynamics  { Cp 4181; Hf 0; }
    transport       { mu 1e-3; Pr 7.0; }
}
EOF
done

###############################################################################
# momentumTransport — laminar for fluid regions
###############################################################################
for r in fluid_1_tubeside fluid_2_shellside; do
    cat > constant/$r/momentumTransport << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      momentumTransport;
}
simulationType laminar;
EOF
done

###############################################################################
# gravity — disabled
###############################################################################
for r in fluid_1_tubeside fluid_2_shellside; do
    cat > constant/$r/g << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    object      g;
}
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 0);
EOF
done

###############################################################################
# regionProperties
###############################################################################
cat > constant/regionProperties << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      regionProperties;
}
regions
(
    fluid
    (
        fluid_1_tubeside
        fluid_2_shellside
    )
    solid
    (
        solid_shell
        solid_baffles
        solid_pipes
    )
);
EOF

###############################################################################
# system/controlDict — steady state uses large pseudo-timestep
###############################################################################
cat > system/controlDict << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
application     foamMultiRegion;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         2000;
deltaT          1;
adjustTimeStep  no;
writeControl    timeStep;
writeInterval   100;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable yes;
EOF

###############################################################################
# top-level fvSchemes and fvSolution
###############################################################################
cat > system/fvSchemes << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes          { default steadyState; }
gradSchemes         { default Gauss linear; }
divSchemes          { default none; }
laplacianSchemes    { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes       { default corrected; }
EOF

cat > system/fvSolution << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers {}
PIMPLE
{
    nOuterCorrectors        20;
    nCorrectors             2;
    nNonOrthogonalCorrectors 1;
    residualControl
    {
        U       1e-4;
        p_rgh   1e-6;
        h       1e-4;
    }
}
EOF

###############################################################################
# system/fvSchemes and fvSolution — solids
###############################################################################
mkdir -p system/solid_pipes

cat > system/solid_pipes/fvSchemes << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes          { default steadyState; }
gradSchemes         { default Gauss linear; }
divSchemes          { default none; }
laplacianSchemes    { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes       { default corrected; }
EOF

cat > system/solid_pipes/fvSolution << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    e
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0.1;
    }
    eFinal
    {
        $e;
        relTol          0;
    }
}
PIMPLE
{
    nNonOrthogonalCorrectors 1;
    residualControl { e 1e-4; }
}
relaxationFactors { equations { e 1; } }
EOF

for r in solid_shell solid_baffles; do
    mkdir -p system/$r
    ln -sf $(pwd)/system/solid_pipes/fvSchemes  system/$r/fvSchemes
    ln -sf $(pwd)/system/solid_pipes/fvSolution system/$r/fvSolution
done

###############################################################################
# system/fvSchemes and fvSolution — fluids (steady state)
###############################################################################
for r in fluid_1_tubeside fluid_2_shellside; do
    mkdir -p system/$r

    cat > system/$r/fvSchemes << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes          { default         steadyState; }
gradSchemes         { default         Gauss linear; }
divSchemes
{
    default                                     none;
    div(phi,U)                                  Gauss upwind;
    div(phi,h)                                  Gauss upwind;
    div(phi,K)                                  Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U)))))         Gauss linear;
}
laplacianSchemes    { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes       { default corrected; }
EOF

    cat > system/$r/fvSolution << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    p_rgh
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
    p_rghFinal
    {
        $p_rgh;
        relTol          0;
    }
    rho
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0;
    }
    rhoFinal { $rho; }
    "(U|h)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0.01;
    }
    "(U|h)Final"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-8;
        relTol          0;
    }
}
PIMPLE
{
    momentumPredictor   yes;
    nOuterCorrectors    20;
    nCorrectors         2;
    nNonOrthogonalCorrectors 1;
    residualControl
    {
        U       1e-4;
        p_rgh   1e-6;
        h       1e-4;
    }
}
relaxationFactors
{
    fields    { p_rgh 0.5; }
    equations { U 0.7; h 0.9; e 0.9; }
}
fieldBounds
{
    p_rgh   1e3  1e8;
    h       1e4  1e8;
    T       200  2000;
}
EOF
done

###############################################################################
# 0/ fields — solids
###############################################################################
mkdir -p 0/solid_pipes
cat > 0/solid_pipes/T << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 300;
boundaryField
{
    ".*"
    {
        type        coupledTemperature;
        Tnbr        T;
        value       uniform 300;
    }
}
EOF

mkdir -p 0/solid_shell
cat > 0/solid_shell/T << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 300;
boundaryField
{
    shell_outer_wall
    {
        type        zeroGradient;
    }
    ".*"
    {
        type        coupledTemperature;
        Tnbr        T;
        value       uniform 300;
    }
}
EOF

mkdir -p 0/solid_baffles
cp 0/solid_pipes/T 0/solid_baffles/T

###############################################################################
# 0/ fields — fluid_1_tubeside
###############################################################################
mkdir -p 0/fluid_1_tubeside

cat > 0/fluid_1_tubeside/T << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 300;
boundaryField
{
    bc_inner_inlet  { type fixedValue; value uniform 350; }
    bc_inner_outlet { type zeroGradient; }
    ".*"
    {
        type        coupledTemperature;
        Tnbr        T;
        value       uniform 300;
    }
}
EOF

cat > 0/fluid_1_tubeside/U << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 -0.001);
boundaryField
{
    bc_inner_inlet  { type fixedValue; value uniform (0 0 -0.001); }
    bc_inner_outlet { type zeroGradient; }
    ".*"            { type noSlip; }
}
EOF

cat > 0/fluid_1_tubeside/p_rgh << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{
    bc_inner_inlet  { type zeroGradient; }
    bc_inner_outlet { type fixedValue; value uniform 1e5; }
    ".*"            { type zeroGradient; }
}
EOF

cat > 0/fluid_1_tubeside/p << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{
    ".*"    { type zeroGradient; }
}
EOF

cat > 0/fluid_1_tubeside/h << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      h;
}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 1254300;
boundaryField
{
    bc_inner_inlet  { type fixedValue; value uniform 1463350; }
    bc_inner_outlet { type zeroGradient; }
    ".*"            { type zeroGradient; }
}
EOF

###############################################################################
# 0/ fields — fluid_2_shellside
###############################################################################
mkdir -p 0/fluid_2_shellside

cat > 0/fluid_2_shellside/T << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
dimensions      [0 0 0 1 0 0 0];
internalField   uniform 280;
boundaryField
{
    bc_outer_inlet  { type fixedValue; value uniform 280; }
    bc_outer_outlet { type zeroGradient; }
    ".*"
    {
        type        coupledTemperature;
        Tnbr        T;
        value       uniform 280;
    }
}
EOF

cat > 0/fluid_2_shellside/U << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 -0.001 0);
boundaryField
{
    bc_outer_inlet  { type fixedValue; value uniform (0 -0.001 0); }
    bc_outer_outlet { type zeroGradient; }
    ".*"            { type noSlip; }
}
EOF

cat > 0/fluid_2_shellside/p_rgh << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{
    bc_outer_inlet  { type zeroGradient; }
    bc_outer_outlet { type fixedValue; value uniform 1e5; }
    ".*"            { type zeroGradient; }
}
EOF

cat > 0/fluid_2_shellside/p << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{
    ".*"    { type zeroGradient; }
}
EOF

cat > 0/fluid_2_shellside/h << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      h;
}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 1170680;
boundaryField
{
    bc_outer_inlet  { type fixedValue; value uniform 1170680; }
    bc_outer_outlet { type zeroGradient; }
    ".*"            { type zeroGradient; }
}
EOF

echo ""
echo "Steady-state case setup complete."
echo "  - ddtSchemes: steadyState"
echo "  - deltaT=1 (pseudo-timestep), endTime=2000 iterations"
echo "  - PIMPLE nOuterCorrectors=50 with residualControl 1e-4"
echo "  - div schemes: upwind (stable for steady state)"
echo "  - relaxationFactors: p_rgh=0.3, U=0.7, h=0.9"
echo "  - U internalField initialised to inlet velocity"
echo "  - Gravity disabled"
echo ""
echo "Workflow:"
echo "  gmshToFoam hx_thick.msh"
echo "  splitMeshRegions -cellZonesOnly -overwrite"
echo "  transformPoints (scale all regions to metres)"
echo "  ./steady_state_setup.sh"
echo "  foamMultiRun 2>&1 | tee log.foamMultiRun"