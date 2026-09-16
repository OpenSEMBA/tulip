.PHONY: build

DOCKER ?= docker
ARTIFACT_DIR ?= dist
BUILD_JOBS ?= 2

# Exports a runnable Ubuntu 26.04 binary bundle to $(ARTIFACT_DIR).
build:
	$(DOCKER) build --build-arg BUILD_JOBS=$(BUILD_JOBS) --file Dockerfile --target artifact --output type=local,dest=$(ARTIFACT_DIR) .
