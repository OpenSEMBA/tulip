# Tulip input data format

- [Tulip input data format](#tulip-input-data-format)
  - [`[options]`.](#options)
  - [`<materials>`](#materials)
    - [`conductor`](#conductor)
    - [`shield`](#shield)
    - [`dielectric`](#dielectric)
    - [`open`](#open)
  - [`<model>`](#model)


Tulip receives a JSON object as an input.

## `[options]`.
Defaults to the options described in ...

`exportParaviewSolution` is defined as `true` in `analysis`, `pulmtln` will also export visualization results for each simulation performed.


## `<materials>`

Is an array of objects. Containing:
- `<id>` an integer identifier with a unique number.
- `<name>` a string with a human readable name.
- `<type>` a string with `conductor`, `shield`, `dielectric`, or `open`. Which are described below.

### `conductor`

### `shield`

### `dielectric`

### `open`

## `<model>`  

`[stepFilename]` defaults to the ```CASE_NAME.step```

`<layers>` which is an array which associates the layers with the different materials.


This object must contain the following entries: 

 