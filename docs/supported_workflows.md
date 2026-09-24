# Supported Tulip Workflows

This document describes the paths implemented by the current code. Tulip analyzes 2D cross sections of multiconductor transmission lines and computes per-unit-length parameters or in-cell reconstruction data for open domains.

## Complete Flow

```mermaid
flowchart TD
    start([Run tulip -i file]) --> format{Valid extension?}
    format -->|.tulip.adapted.json| adapted[Read Gmsh mesh, materials, and options]
    format -->|.tulip.input.json| input[Read input JSON]
    format -->|Other| invalidExtension[[Error: unsupported extension]]

    input --> inputValid{materials and layers are arrays?}
    inputValid -->|No| invalidInput[[Error: required sections missing]]
    inputValid -->|Yes| materialIds{Does every layer reference an existing materialId?}
    materialIds -->|No| invalidMaterial[[Error: unknown materialId]]
    materialIds -->|Yes| step[Resolve configured STEP or CASE.step]
    step --> stepExists{Does the STEP exist?}
    stepExists -->|No| missingStep[[Error: STEP not found]]
    stepExists -->|Yes| import[Import STEP into Gmsh]
    import --> names{Do JSON layers and STEP match in both directions?}
    names -->|No| namesError[[Error: unmatched layer]]
    names -->|Yes| classify[Classify conductor, shield, dielectric, and open]

    classify --> geometry[Resolve overlaps and build vacuum]
    geometry --> mesh[Create physical groups and 2D mesh]
    mesh --> meshValid{Are duplicate nodes left?}
    meshValid -->|Yes| meshError[[Error: mesh has duplicate nodes]]
    meshValid -->|No| writeMesh[Write CASE.msh and generate adapted model]
    writeMesh --> adapted

    adapted --> parseModel[Validate and load model]
    parseModel --> conductors{Is there at least one conductor?}
    conductors -->|No| noConductors[[Error: model has no conductors]]
    conductors -->|Yes| electrostatic[Solve electrostatics for every conductor]
    electrostatic --> dielectric{Does the model contain dielectrics?}
    dielectric -->|Yes| magnetic[Solve magnetostatics with relative permittivity equal to 1]
    dielectric -->|No| reuse[Reuse electrostatic solution for magnetostatics]
    magnetic --> domains[Separate domains by connectivity]
    reuse --> domains
    domains --> output[Generate CASE.tulip.out.json]
```

## Geometry Adaptation

`shield` is supplied to the solver as a conductor, with additional transfer-impedance information when declared. A dielectric without `relativePermittivity` uses 1.0. A conductor without resistance or conductivity is treated as PEC.

```mermaid
flowchart TD
    begin([2D STEP entities]) --> types[Associate every layer with its material]
    types --> overlaps[Resolve dielectric-dielectric overlaps]
    overlaps --> priority{Which dielectric takes precedence?}
    priority -->|Higher relative permittivity| highEpsilon[Keep dielectric with higher epsilon_r]
    priority -->|Same permittivity| lexical[Keep lexicographically greater name]
    highEpsilon --> cutConductors[Cut conductors out of dielectrics and vacuum]
    lexical --> cutConductors

    cutConductors --> openCase{Is this an open case?}
    openCase -->|No| closed[Closed domain: root conductor minus interior elements]
    openCase -->|Yes, and an open layer exists| declaredOpen[Cut open layer with conductors and dielectrics]
    openCase -->|Yes, and no open layer exists| automaticOpen[Create rectangular InnerRegion and circular OuterRegion]
    automaticOpen --> openBoundary[Create OpenBoundary at the outer edge]
    closed --> interfaces[Fragment dielectric and vacuum for conformal interfaces]
    declaredOpen --> interfaces
    openBoundary --> interfaces
    interfaces --> boundaries[Convert conductors into boundaries]
    boundaries --> physical[Assign Gmsh physical groups]
```

Tulip classifies a case as open if there is a single `open`, the geometry has no single root conductor enclosing elements, the root is a dielectric, or in certain shield configurations with more than two non-intersecting conductors. Otherwise, it attempts to construct a closed domain from the root conductor.

## Materials and Properties

```mermaid
flowchart LR
    material{Material type} -->|conductor| conductor[Conductor for the solver]
    material -->|shield| shield[Conductor with isShield]
    material -->|dielectric| dielectric[Dielectric domain]
    material -->|open| open[Open boundary]

    conductor --> resistance{Defines resistancePerMeter and conductivity?}
    resistance -->|Both| invalidResistance[[Error: mutually exclusive properties]]
    resistance -->|Only resistancePerMeter| directR[Use declared resistance]
    resistance -->|Only conductivity| calculateR[Calculate R per meter = 1 divided by sigma and area]
    resistance -->|Neither| pec[Zero resistance: PEC]
    calculateR --> validSigma{Are conductivity and area positive?}
    validSigma -->|No| invalidSigma[[Error: invalid conductivity or area]]
    validSigma -->|Yes| directR

    shield --> transfer{Is transfer impedance present?}
    transfer -->|Declared object| transferObject[Use resistive, inductive, and direction terms]
    transfer -->|Individual terms| transferFields[Build transfer object]
    transfer -->|No| noTransfer[No transfer impedance]
    dielectric --> epsilon[Use declared epsilon_r or 1.0]
```

## Per-Domain Solution and Results

For each conductor, Tulip applies 1 V to the excited conductor and 0 V to the others. It builds the generalized capacitance matrix from the resulting charges. The magnetostatic solution ignores dielectric permittivities; `L` is derived from the inverse of the equivalent vacuum capacitance.

```mermaid
flowchart TD
    solved([Electric and magnetic fields solved]) --> open{Does the model have an open boundary?}
    open -->|No| closedDomain[Closed domain]
    open -->|Yes| inspectDomain[Inspect every connected domain]

    closedDomain --> closedPUL[Choose ground conductor and remove its row and column]
    closedPUL --> pul[shieldedMultiwire output]

    inspectDomain --> grounded{Does the domain have a ground conductor?}
    grounded -->|Yes| nestedPUL[Build local PUL matrices]
    grounded -->|No| inCell[Calculate floating potentials and multipolar expansion]
    nestedPUL --> transferShield{Is ground a shield?}
    transferShield -->|Yes| addTransfer[Attach transfer impedance]
    transferShield -->|No| localPUL[No transfer impedance]
    addTransfer --> pul
    localPUL --> pul
    inCell --> openOutput[unshieldedMultiwire output]

    pul --> json[Serialize FDTD JSON with material associations]
    openOutput --> json
    json --> optional[Optional: export Phi, E, and D for ParaView]
```

| Geometry | Expected result |
| --- | --- |
| Coaxial cable or several conductors inside an outer conductor | Closed domain and `R`, `L`, `C` matrices without the reference conductor |
| Conductors without an enclosure, with an `open` boundary | Open domain and `unshieldedMultiwire` reconstruction |
| Conductors without an enclosure or `open` boundary | Open domain with automatically generated inner/outer regions and boundary |
| Conductors inside a shield, with open exterior | Internal PUL domain plus external open domain |
| Geometry with several nested shields | One output per connected domain; every internal domain uses its ground |
| Touching or overlapping dielectrics | Conformal interfaces; overlapping dielectrics use the one with higher permittivity |

## Relevant Restrictions and Errors

```mermaid
flowchart LR
    request{Requested operation} -->|Simple PUL| oneDomain{Is there only one domain?}
    oneDomain -->|Yes| pulOk[Calculate PUL]
    oneDomain -->|No| pulError[[Error: simple PUL accepts one domain only]]

    request -->|In-cell potentials| isOpen{Is the model open?}
    isOpen -->|Yes| inCellOk[Calculate in-cell values]
    isOpen -->|No| inCellError[[Error: in-cell accepts open models only]]

    request -->|Closed-problem C| enough{Are there at least two conductors?}
    enough -->|Yes| capacitanceOk[Calculate C]
    enough -->|No| capacitanceError[[Error: a closed problem requires two conductors]]
```

The test cases in `test/adapter` and `test/driver` cover, among others, empty and partially filled coaxial cables, shielded and unshielded multiwire configurations, dielectrics, open boundaries, nested geometry, and in-cell parameters.
