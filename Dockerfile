# syntax=docker/dockerfile:1
FROM public.ecr.aws/docker/library/ubuntu:26.04 AS build-base

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_JOBS=2
ENV VCPKG_ROOT=/opt/vcpkg
ENV VCPKG_INSTALLED_DIR=/opt/vcpkg/installed
ENV _VCPKG_INSTALLED_DIR=/opt/vcpkg/installed
ENV VCPKG_MAX_CONCURRENCY=${BUILD_JOBS}

# These packages bootstrap vcpkg and provide the host tools required by gmsh[occ].
RUN apt-get update && apt-get install --yes --no-install-recommends \
    autoconf \
    autoconf-archive \
    automake \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    git \
    libalglib-dev \
    libann-dev \
    libblas-dev \
    libcgns-dev \
    libfontconfig1-dev \
    libfltk1.3-dev \
    libfreetype-dev \
    libgmsh-dev \
    libgl2ps-dev \
    libglu1-mesa-dev \
    libglut-dev \
    libgmp-dev \
    libjpeg-dev \
    libtool \
    libx11-dev \
    libxext-dev \
    libxi-dev \
    libxmu-dev \
    libxt-dev \
    libmetis-dev \
    libocct-data-exchange-dev \
    libocct-foundation-dev \
    libocct-modeling-algorithms-dev \
    libocct-modeling-data-dev \
    libocct-ocaf-dev \
    libpng-dev \
    mesa-common-dev \
    ninja-build \
    patchelf \
    pkg-config \
    python3 \
    tar \
    unzip \
    zip \
    && ln -s "$(ldconfig -p | perl -ne 'print $1 and exit if /libvoro\+\+\.so.* => (.*)/')" /usr/lib/x86_64-linux-gnu/libvoro++.so \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/microsoft/vcpkg.git "$VCPKG_ROOT" \
    && git -C "$VCPKG_ROOT" checkout eed289f6e06a5e7a5c9e6a729671b0a56af7dd69 \
    && "$VCPKG_ROOT/bootstrap-vcpkg.sh" -disableMetrics

WORKDIR /src
COPY vcpkg.json ./
COPY triplets ./triplets
RUN "$VCPKG_ROOT/vcpkg" install --triplet x64-linux --overlay-triplets=/src/triplets

COPY . .

FROM build-base AS build

RUN cmake --preset gnu -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build --parallel "$BUILD_JOBS"

# Bundle project dependencies, relying on the Ubuntu 26.04 runtime for its
# standard C/C++ and desktop libraries.
RUN set -eux; \
    mkdir -p /bundle /bundle/lib; \
    cp build/tulip /bundle/; \
    ldd build/tulip | awk '/=> \/[^ ]+/ { print $3 } /^\// { print $1 }' | sort -u | \
    while IFS= read -r library; do \
        case "$(basename "$library")" in \
            libc.so.*|libm.so.*|libresolv.so.*|libstdc++.so.*|libgcc_s.so.*|libgomp.so.*|\
            libGL.so.*|libGLdispatch.so.*|libGLX.so.*|libGLU.so.*|libOpenGL.so.*|libglut.so.*|\
            libX11.so.*|libXau.so.*|libxcb.so.*|libXcursor.so.*|libXdmcp.so.*|libXext.so.*|\
            libXfixes.so.*|libXft.so.*|libXi.so.*|libXinerama.so.*|libXrender.so.*|libXxf86vm.so.*|\
            libfontconfig.so.*|libfreetype.so.*|libexpat.so.*|libjpeg.so.*|libpng*.so.*|\
            libbz2.so.*|libz.so.*|libzstd.so.*) continue ;; \
        esac; \
        cp "$library" /bundle/lib/; \
    done; \
    patchelf --set-rpath '$ORIGIN/lib' /bundle/tulip; \
    for library in /bundle/lib/*; do patchelf --set-rpath '$ORIGIN' "$library" || true; done; \
    printf '%s\n' \
        'Runtime target: Ubuntu 26.04 with its standard C/C++ and desktop libraries installed.' \
        'The lib directory contains only bundled third-party dependencies.' \
        > /bundle/RUNTIME-REQUIREMENTS.txt

# This target retains the compiled executable and CTest suite so that `make
# test` can run without rebuilding the project.
FROM build AS test

# Debuggable test image. Keep artifacts outside /src so the workspace can be
# mounted at /workspace by VS Code without hiding the compiled binaries.
FROM build-base AS test-debug
RUN apt-get update && apt-get install --yes --no-install-recommends gdb \
    && rm -rf /var/lib/apt/lists/* \
    && cmake -S /src -B /build-debug -G Ninja \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" \
        -DVCPKG_OVERLAY_TRIPLETS=/src/triplets \
        -DTULIP_USE_MFEM_AS_SUBDIRECTORY=ON \
        -DMFEM_ENABLE_TESTING=OFF \
        -DMFEM_USE_OPENMP=OFF \
    && cmake --build /build-debug --parallel "$BUILD_JOBS"

FROM public.ecr.aws/docker/library/ubuntu:26.04 AS runtime
RUN apt-get update && apt-get install --yes --no-install-recommends libgmsh4.14 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /bundle/ /opt/tulip/
ENTRYPOINT ["/opt/tulip/tulip"]

FROM scratch AS artifact
COPY --from=build /bundle/ /
