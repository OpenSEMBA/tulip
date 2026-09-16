.PHONY: help build test run

DOCKER ?= docker
ARTIFACT_DIR ?= dist
BUILD_JOBS ?= 2
IMAGE ?= tulip:latest
TEST_IMAGE ?= tulip-test:latest
RUN_FILE ?= $(or $(FILE),$(filter-out run,$(MAKECMDGOALS)))

help:
	@printf '%s\n' \
		'Available commands:' \
		'  make build                       Build the runtime and test images, and export dist/.' \
		'  make test                        Run all CTest tests from the prebuilt test image.' \
		'  make run FILE=path/to/input.json Run Tulip in a container and create a results folder beside the input.' \
		'  make help                        Show this help.'

# Builds the reusable container image and exports a runnable Ubuntu 26.04
# binary bundle to $(ARTIFACT_DIR).
build:
	$(DOCKER) build --build-arg BUILD_JOBS=$(BUILD_JOBS) --file Dockerfile --target runtime --tag $(IMAGE) .
	$(DOCKER) build --build-arg BUILD_JOBS=$(BUILD_JOBS) --file Dockerfile --target test --tag $(TEST_IMAGE) .
	$(DOCKER) build --build-arg BUILD_JOBS=$(BUILD_JOBS) --file Dockerfile --target artifact --output type=local,dest=$(ARTIFACT_DIR) .

# Runs every CTest test from the already-built image. Build it on first use.
test:
	@$(DOCKER) image inspect $(TEST_IMAGE) >/dev/null 2>&1 || $(MAKE) build
	$(DOCKER) run --rm --entrypoint ctest $(TEST_IMAGE) --test-dir /src/build --output-on-failure

# Usage: make run FILE=path/to/case.tulip.input.json
# `make run path/to/case.tulip.input.json` is also accepted for paths without spaces.
run:
	@test -n "$(RUN_FILE)" || { printf '%s\n' 'Usage: make run FILE=path/to/case.tulip.input.json'; exit 2; }
	@test -f "$(RUN_FILE)" || { printf 'Input file not found: %s\n' "$(RUN_FILE)"; exit 2; }
	@$(DOCKER) image inspect $(IMAGE) >/dev/null 2>&1 || $(MAKE) build
	@input_file="$$(realpath "$(RUN_FILE)")"; \
	input_dir="$$(dirname "$$input_file")"; \
	input_name="$$(basename "$$input_file")"; \
	case_name="$${input_name%.tulip.input.json}"; \
	case_name="$${case_name%.tulip.adapted.json}"; \
	output_dir="$$(mktemp -d "$$input_dir/$${case_name}.tulip-output.XXXXXX")"; \
	$(DOCKER) run --rm --user "$$(id -u):$$(id -g)" \
		--mount "type=bind,src=$$input_dir,dst=/source,readonly" \
		--mount "type=bind,src=$$output_dir,dst=/output" \
		--workdir /output \
		-e INPUT_NAME="$$input_name" \
		--entrypoint sh $(IMAGE) -ec 'mkdir -p /tmp/input && cp -a /source/. /tmp/input/ && exec /opt/tulip/bin/tulip -i "/tmp/input/$$INPUT_NAME" -o /output' && \
	printf 'Results written to: %s\n' "$$output_dir"

# Treat a positional input-file argument as data for the run target.
%:
	@:
