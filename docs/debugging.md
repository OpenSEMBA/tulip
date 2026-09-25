# Debugging

Tulip is built and debugged inside Docker containers. The VS Code configuration
creates a temporary container with GDB and `ptrace` permissions, so no local
build directory is required.

## Requirements

- Docker must be available to the current user.
- VS Code must have the Microsoft C/C++ extension (`ms-vscode.cpptools`) installed.

## Debug Image

Build the image with debug symbols and GDB:

```shell
make build-debug
```

The `tulip-test-debug:latest` image contains the `tulip` and `driver_tests`
executables built in `Debug` mode. VS Code runs the same command before every
debug session; Docker reuses its layers when nothing has changed.

## Debugging a Test

1. Open the **Run and Debug** view in VS Code.
2. Select `LauncherTest dielectric (Docker gdb)`.
3. Set breakpoints and press F5.

The configuration runs only:

```text
LauncherTest.two_conductors_with_dielectric_material_associations
```

The GoogleTest filter is defined in `.vscode/launch.json`. Change
`--gtest_filter=SuiteName.TestName` to debug another test.

## Debugging an Input Case

1. Open `.vscode/launch.json`.
2. In `Tulip input case (Docker gdb)`, edit the `-i` and `-o` arguments.
3. Select that configuration in **Run and Debug** and press F5.

For example:

```json
"args": [
    "-i",
    "testData/two_conductors_with_dielectric/two_conductors_with_dielectric.tulip.input.json",
    "-o",
    "/tmp/tulip-debug-output"
]
```

The repository root is mounted at `/workspace` inside the container and is the
working directory. Relative input paths are resolved from that location. The
default output is stored inside the container; use a path under `/workspace` to
keep it on the host.

## Debugging from the Terminal

To run the same test without VS Code:

```shell
docker run --rm \
  --cap-add=SYS_PTRACE \
  --security-opt seccomp=unconfined \
  --mount "type=bind,src=$(pwd),dst=/workspace" \
  --workdir /workspace \
  tulip-test-debug:latest \
  gdb /build-dbg/driver_tests
```

In GDB, start the case with:

```gdb
run --gtest_filter=LauncherTest.two_conductors_with_dielectric_material_associations
```
