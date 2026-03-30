# Tulip input data format

- [Tulip input data format](#tulip-input-data-format)
  - [`[analysis]`.](#analysis)
    - [adapter options](#adapter-options)
    - [driver options](#driver-options)
  - [`<materials>`](#materials)
    - [`conductor`](#conductor)
    - [`shield`](#shield)
    - [`dielectric`](#dielectric)
    - [`open`](#open)
  - [`<model>`](#model)

Tulip receives a JSON object as an input with the entries described below. Square brackets indicate that the entry is optional and a default value will be assumed, angle brackets indicate that the entry is mandatory. 

Unless specified otherwise all units are assumed to be in SI-MKS.

## `[analysis]`.
This object contains options to control the adapter and driver behaviors. 
Adapter is in 

### adapter options

### driver options
- `[exportParaviewSolution]` can be `true` or `false`. Exports visualization results for each simulation performed.
- 


## `<materials>`
These materials are associated with `model` `layers` to define regions with different material properties.
They are defined by an array of JSON objects with:
- `<id>` an integer identifier with a unique number.
- `<name>` a string with a human readable name.
- `<type>` a string with `conductor`, `shield`, `dielectric`, or `open`. Depending on their type  they may need additional entries, as described below.

### `conductor`
A conductor material can contain either `[resistancePerMeter]` or  `[conductivity]`, but not both. If `conductivity` is defined, it will be used to compute a resistance per meter for the conductor.
If none is specified, the conductor is considered to be a perfect electric conductor.

Conductors associated with layers cannot intersect any other layer which is also associated with a conductor. They can however intersect layers associated with `shield`, `dielectric`, or `open`, taking precedence over them.

### `shield`
A shield is a special kind of conductor which defines a transfer impedance model, specified with:
- `[resistancePerMeter]` defined by a real representing transfer impedance resistance. Defaults to `0.0`
- `[inductancePerMeter]` defined by a real representing transfer impedance inductance. Defaults to `0.0`.
- `[direction]` which can be `both`, `inwards`, or `outwards`. Indicating the type of coupling considered. Defaults to `both` meaning that fields can couple from the exterior to interior and the other way round.

Shield layers can never intersect any other `shield` layer. They can be used in two ways:
+ To represent a _closed problem_, or the interior of a shielded domain. In that case:
  - Shields can intersect conductor or dielectric layers.
  - The problem cannot contain any open layer. 
  - The shield is assumed to be the ground conductor. 
  - There can only be one shield.
+ In _open problems_ they act as the boundary between an internal and external domain. In that case:
  - They follow the same rules as conductors about intersections
  - Must have a non-null area.
  - The shield is assumed to be the ground of the domain which it encloses and is one conductor more for the domain to which it belongs.

### `dielectric`
A dielectric is defined with a `[relativePermittivity]` which defaults to `1.0`. If two dielectrics are intersecting, the one with the highest permittivity takes precedence.

### `open`
An `open` material serves to specify the computational boundary of the problem. It must intersect every other material layer. If no open boundary is specified for an open problem, one is computed automatically, together with _inner_ and _outer_ regions used to extract the unshielded multiwire coefficients.


## `<model>`  
This object can contain the following entries:
+ `[stepFilename]` defaults to the ```CASE_NAME.step```
  + `<layers>` which is an array which associates the layers present in the `.step` file with the different `materials`. Each layer is specified by:
  - `<name>` which must match exactly the name of the corresponding layer within the `.step` file. It must be unique.
  - `<id>` which is an integer unique identifier which will be used to order the results for the calculated PUL matrices.
  - `<materialId>` which must match an `id` from a material in the list of `materials`
  


 